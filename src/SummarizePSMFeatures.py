#!/usr/bin/env python3
"""Summarize published MSBooster PIN or OpenMS PSMFeatureExtractor features."""

import argparse
import csv
import re
from collections import defaultdict

import pyopenms as oms

from QCMetrics import Metrics, as_text, finite


def target_decoy(value):
    """Normalize a PIN label or idXML target/decoy annotation."""
    label = as_text(value).strip().lower()
    if label in ("1", "+1", "target"):
        return "target"
    if label in ("-1", "decoy"):
        return "decoy"
    return "unknown"


def selected_feature(name):
    """Select interpretable feature families, excluding identifiers and scores."""
    normalized = re.sub(r"[^a-z0-9]", "", name.lower())
    return (normalized in ("deltamass", "dm", "absdm") or
            ("rt" in normalized and ("delta" in normalized or "diff" in normalized)) or
            ("retention" in normalized and ("delta" in normalized or "diff" in normalized)) or
            any(word in normalized for word in ("spectralangle", "spectralentropy",
                                               "spectralcorrelation", "spectralcosine")))


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
            declared.extend(as_text(search.getMetaValue("extra_features")).split(","))
    columns = list(dict.fromkeys(name.strip() for name in declared if name.strip()))
    rows = []
    for peptide in peptides:
        for hit in peptide.getHits():
            label = (target_decoy(hit.getMetaValue("target_decoy"))
                     if hit.metaValueExists("target_decoy") else "unknown")
            values = {name: as_text(hit.getMetaValue(name))
                      for name in columns if hit.metaValueExists(name)}
            rows.append((label, values))
    return columns, rows


def summarize(args):
    """Write compact feature distributions by target/decoy class."""
    metrics = Metrics(args.dataset_id, args.sample_name, args.input, "psm_features",
                      source_process=args.route, search_engine=args.search_engine,
                      fraction_id=args.fraction_id)
    columns, rows = (pin_rows(args.input) if args.route == "MSBooster"
                     else idxml_rows(args.input))
    groups = ("all", "target", "decoy", "unknown")
    by_group = defaultdict(list)
    for label, values in rows:
        by_group["all"].append(values)
        by_group[label].append(values)
    for group in groups:
        metrics.add("feature_psm_count", len(by_group[group]), group=group)
    feature_columns = [name for name in columns if selected_feature(name)]
    metric_names = [metric_key(name) for name in feature_columns]
    if len(metric_names) != len(set(metric_names)):
        raise ValueError("Feature column names collide after metric-name normalization")
    for name in feature_columns:
        unit = "m/z" if name.lower() == "deltamass" else ""
        for group in groups:
            values = [row.get(name) for row in by_group[group]]
            metrics.distribution(metric_key(name), values, unit, group)
            if ("delta" in name.lower() or name.lower() == "dm" or
                    ("diff" in name.lower() and "rt" in name.lower())):
                metrics.distribution(metric_key(name) + "_absolute",
                                     (abs(number) if (number := finite(value)) is not None
                                      else None for value in values), unit, group)
    mass_column = next((name for name in columns if name.lower() == "deltamass"), None)
    calc_column = next((name for name in columns if name.lower() == "calcmass"), None)
    rt_columns = [name for name in feature_columns
                  if ("rt" in name.lower() or "retention" in name.lower()) and
                  ("delta" in name.lower() or "diff" in name.lower())]
    for group in groups:
        if mass_column and calc_column:
            errors = []
            for row in by_group[group]:
                delta, theoretical = finite(row.get(mass_column)), finite(row.get(calc_column))
                errors.append(delta / theoretical * 1_000_000
                              if delta is not None and theoretical is not None and theoretical > 0
                              else None)
            metrics.distribution("feature_mass_error_ppm_signed", errors, "ppm", group)
            metrics.distribution("feature_mass_error_ppm_absolute",
                                 (abs(value) if value is not None else None for value in errors),
                                 "ppm", group)
        else:
            metrics.unavailable("feature_mass_error_ppm_signed_p50", "ppm", group)
            metrics.unavailable("feature_mass_error_ppm_absolute_p50", "ppm", group)
        if not rt_columns:
            metrics.unavailable("feature_rt_delta_p50", group=group)
    metrics.write(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("input", "output", "dataset-id", "sample-name", "search-engine",
                     "route"):
        parser.add_argument("--" + argument, required=True)
    parser.add_argument("--fraction-id", default="")
    summarize(parser.parse_args())
