#!/usr/bin/env python3
"""Summarize published MSBooster PIN or OpenMS PSMFeatureExtractor features."""

import argparse
import csv
import re
from collections import Counter, defaultdict

import pyopenms as oms

from QCMetrics import Metrics, as_text as text, finite

SIMILARITY_FEATURES = (
    "spectralangle",
    "spectralentropy",
    "spectralcorrelation",
    "spectralcosine",
)


def target_decoy(value):
    """Normalize a PIN label or idXML target/decoy annotation."""
    label = text(value).strip().lower()
    if label in ("1", "+1", "target"):
        return "target"
    if label in ("-1", "decoy"):
        return "decoy"
    return "unknown"


def selected_feature(name):
    """Select interpretable feature families, excluding identifiers and scores."""
    normalized = re.sub(r"[^a-z0-9]", "", name.lower())
    mass_error = normalized in ("deltamass", "dm", "absdm")
    retention_error = (
        ("rt" in normalized or "retention" in normalized)
        and ("delta" in normalized or "diff" in normalized)
    )
    spectral_similarity = any(word in normalized for word in SIMILARITY_FEATURES)
    return mass_error or retention_error or spectral_similarity


def metric_key(name):
    """Keep source column names recognizable in the metric key."""
    return "feature_" + re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def pin_rows(path):
    """Read tab-delimited edited PIN rows and their source columns."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if "Label" not in columns:
            raise ValueError("MSBooster PIN is missing Label")
        rows = [(target_decoy(row["Label"]), row) for row in reader]
    return columns, rows


def idxml_rows(path):
    """Read hit-level features listed by PSMFeatureExtractor."""
    proteins = []
    peptides = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(path, proteins, peptides)
    declared = []
    for protein in proteins:
        search = protein.getSearchParameters()
        if search.metaValueExists("extra_features"):
            declared.extend(text(search.getMetaValue("extra_features")).split(","))
    columns = list(dict.fromkeys(name.strip() for name in declared if name.strip()))
    rows = []
    for peptide in peptides:
        for hit in peptide.getHits():
            label = (target_decoy(hit.getMetaValue("target_decoy"))
                     if hit.metaValueExists("target_decoy") else "unknown")
            values = {name: text(hit.getMetaValue(name))
                      for name in columns if hit.metaValueExists(name)}
            rows.append((label, values))
    return columns, rows


def pin_category(row, fields, kind):
    """Decode a PIN enzyme or one-hot charge category, or its data issue."""
    if kind == "missed_cleavages":
        value = finite(row.get(fields[0]))
        if value is None:
            return "missing" if not row.get(fields[0]) else "invalid"
        return str(int(value)) if value >= 0 and value.is_integer() else "invalid"
    if kind == "cleavage_specificity":
        values = [finite(row.get(field)) for field in fields]
        if any(value is None for value in values):
            return "missing" if any(not row.get(field) for field in fields) else "invalid"
        if any(value not in (0, 1) for value in values):
            return "invalid"
        return ("full", "semi", "none")[2 - int(sum(values))]
    values = [finite(row.get(field)) for field in fields]
    if any(value is None for value in values):
        return "missing" if any(not row.get(field) for field in fields) else "invalid"
    if any(value not in (0, 1) for value in values) or sum(values) != 1:
        return "invalid"
    return str(int(re.fullmatch(r"charge([0-9]+)", fields[values.index(1)], re.I).group(1)))


def summarize_pin_categories(metrics, columns, by_group, groups):
    """Count PIN PSMs by missed cleavages, enzyme specificity, and charge."""
    charge_fields = [name for name in columns if re.fullmatch(r"charge[0-9]+", name, re.I)]
    specifications = (
        ("missed_cleavages", ["enzInt"] if "enzInt" in columns else [], ()),
        ("cleavage_specificity", ["enzN", "enzC"]
         if all(name in columns for name in ("enzN", "enzC")) else [],
         ("full", "semi", "none")),
        ("charge", charge_fields, ()),
    )
    for kind, fields, fixed_categories in specifications:
        prefix = f"feature_{kind}"
        if not fields:
            for group in groups:
                metrics.unavailable(
                    name=prefix + "_psm_count",
                    unit="count",
                    group=group,
                )
            continue
        counts = {}
        for group in groups:
            counts[group] = Counter(
                pin_category(row, fields, kind)
                for row in by_group[group]
            )
        categories = set(fixed_categories)
        for group in groups:
            categories.update(counts[group])
        categories.difference_update(("missing", "invalid"))
        for group in groups:
            for category in sorted(
                categories,
                key=lambda item: (not item.isdigit(),
                                  int(item) if item.isdigit() else item),
            ):
                metrics.add(
                    name=f"{prefix}_{category}_psm_count",
                    value=counts[group][category],
                    group=group,
                )
            for issue in ("missing", "invalid"):
                metrics.add(
                    name=f"{prefix}_{issue}_psm_count",
                    value=counts[group][issue],
                    group=group,
                )


class PSMFeaturesSummarizer:
    """Summarize published PSM features by target/decoy class."""

    GROUPS = ("all", "target", "decoy", "unknown")

    def __init__(
        self,
        input_path,
        dataset_id,
        sample_name,
        search_engine,
        route,
        fraction_id,
    ):
        self.input_path = input_path
        self.route = route
        self.metrics = Metrics(
            dataset_id=dataset_id,
            layer="psm_features",
            source=input_path,
            sample_name=sample_name,
            source_process=route,
            search_engine=search_engine,
            fraction_id=fraction_id,
        )
        self.columns = []
        self.by_group = defaultdict(list)

    def load_rows(self):
        """Read a PIN or idXML and partition feature rows by label."""
        self.columns, rows = (pin_rows(self.input_path) if self.route == "MSBooster"
                              else idxml_rows(self.input_path))
        for label, values in rows:
            self.by_group["all"].append(values)
            self.by_group[label].append(values)

    def summarize_counts_and_categories(self):
        """Count PSMs and MSBooster enzyme or charge categories."""
        for group in self.GROUPS:
            self.metrics.add(
                name="feature_psm_count",
                value=len(self.by_group[group]),
                group=group,
            )
        if self.route == "MSBooster":
            summarize_pin_categories(
                metrics=self.metrics,
                columns=self.columns,
                by_group=self.by_group,
                groups=self.GROUPS,
            )

    def summarize_source_features(self):
        """Keep each selected source feature and its absolute deviation where relevant."""
        feature_columns = [name for name in self.columns if selected_feature(name)]
        metric_names = [metric_key(name) for name in feature_columns]
        if len(metric_names) != len(set(metric_names)):
            raise ValueError("Feature column names collide after metric-name normalization")
        for name in feature_columns:
            unit = "m/z" if name.lower() == "deltamass" else ""
            for group in self.GROUPS:
                values = [row.get(name) for row in self.by_group[group]]
                self.metrics.distribution(
                    name=metric_key(name),
                    values=values,
                    unit=unit,
                    group=group,
                )
                if ("delta" in name.lower() or name.lower() == "dm" or
                        ("diff" in name.lower() and "rt" in name.lower())):
                    self.metrics.distribution(
                        name=metric_key(name) + "_absolute",
                        values=(abs(number) if (number := finite(value)) is not None
                                else None for value in values),
                        unit=unit,
                        group=group,
                    )

    def summarize_mass_error(self):
        """Derive ppm mass error when both MSBooster mass fields exist."""
        mass_column = next((name for name in self.columns if name.lower() == "deltamass"), None)
        calc_column = next((name for name in self.columns if name.lower() == "calcmass"), None)
        for group in self.GROUPS:
            if mass_column and calc_column:
                errors = []
                for row in self.by_group[group]:
                    delta = finite(row.get(mass_column))
                    theoretical = finite(row.get(calc_column))
                    if delta is not None and theoretical is not None and theoretical > 0:
                        error = delta / theoretical * 1_000_000
                    else:
                        error = None
                    errors.append(error)
                self.metrics.distribution(
                    name="feature_mass_error_ppm_signed",
                    values=errors,
                    unit="ppm",
                    group=group,
                )
                self.metrics.distribution(
                    name="feature_mass_error_ppm_absolute",
                    values=(abs(value) if value is not None else None for value in errors),
                    unit="ppm",
                    group=group,
                )
            else:
                self.metrics.unavailable_distribution(
                    name="feature_mass_error_ppm_signed",
                    unit="ppm",
                    group=group,
                )
                self.metrics.unavailable_distribution(
                    name="feature_mass_error_ppm_absolute",
                    unit="ppm",
                    group=group,
                )

    def run(self, output_path):
        """Read the features and write their QC summaries."""
        self.load_rows()
        self.summarize_counts_and_categories()
        self.summarize_source_features()
        self.summarize_mass_error()
        self.metrics.write(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("input", "output", "dataset-id", "sample-name", "search-engine",
                     "route"):
        parser.add_argument("--" + argument, required=True)
    parser.add_argument("--fraction-id", default="")
    args = parser.parse_args()
    summarizer = PSMFeaturesSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        sample_name=args.sample_name,
        search_engine=args.search_engine,
        route=args.route,
        fraction_id=args.fraction_id,
    )
    summarizer.run(args.output)
