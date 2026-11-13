#!/usr/bin/env python3
"""Summarize the final aggregated DIA-NN report by sample."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from QCMetrics import Metrics, finite


def first(columns, candidates):
    """Return the first available column name from a list of candidates."""
    return next((name for name in candidates if name in columns), None)


def run_lookup(path):
    """Map unambiguous DIA-NN run labels to biological samples."""
    lookup = {}
    pairs = set()
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        if not {"sample_name", "fraction_id", "path"} <= columns:
            raise ValueError("mzML manifest needs sample_name,fraction_id,path")
        for row in reader:
            sample = row["sample_name"].strip()
            fraction = row["fraction_id"].strip()
            pair = (sample, fraction)
            if not sample or not fraction or pair in pairs:
                raise ValueError(f"Invalid or duplicate sample/fraction pair: {pair}")
            pairs.add(pair)
            paths = [Path(row["path"])]
            if row.get("acquisition_path"):
                paths.append(Path(row["acquisition_path"]))
            for path_value in paths:
                for key in (path_value.name, path_value.stem):
                    if key in lookup and lookup[key] != sample:
                        raise ValueError(f"Ambiguous DIA-NN run key: {key}")
                    lookup[key] = sample
    for sample, _ in pairs:
        if sample in lookup and lookup[sample] != sample:
            raise ValueError(f"Ambiguous DIA-NN run key: {sample}")
        lookup[sample] = sample
    return lookup


class DIANNReportSummarizer:
    """Summarize final aggregated DIA-NN quantities at the sample level."""

    def __init__(self, input_path, sample_csv, dataset_id):
        self.input_path = input_path
        self.lookup = run_lookup(sample_csv)
        self.metrics = Metrics(dataset_id, "", input_path, "quantification")
        self.precursors = defaultdict(set)
        self.proteins = defaultdict(set)
        self.precursor_values = defaultdict(dict)
        self.protein_values = defaultdict(dict)
        self.observed_samples = set()
        self.columns = {}

    def select_columns(self, columns):
        """Identify the current DIA-NN report columns used by QC."""
        self.columns = {
            "run": first(columns, ("Run", "File.Name")),
            "precursor": first(columns, ("Precursor.Id", "Precursor.ID")),
            "precursor_quantity": first(columns, (
                "Precursor.Quantity", "Precursor.Normalised", "Precursor.Translated")),
            "protein": first(columns, ("Protein.Group", "Protein.Ids")),
            "protein_quantity": first(columns, ("PG.Quantity", "PG.Normalised")),
        }
        if not self.columns["run"] or not self.columns["precursor"]:
            raise ValueError("DIA-NN report needs Run and Precursor.Id")

    def sample_for_row(self, row):
        """Map a DIA-NN run label to exactly one biological sample."""
        run_name = row[self.columns["run"]]
        candidates = (run_name, Path(run_name).name, Path(run_name).stem)
        matches = {self.lookup[key] for key in candidates if key in self.lookup}
        if len(matches) != 1:
            raise ValueError(f"Unknown or ambiguous DIA-NN run: {run_name}")
        return matches.pop()

    def record_row(self, row):
        """Collect distinct identities and positive quantities for one report row."""
        sample = self.sample_for_row(row)
        self.observed_samples.add(sample)
        precursor = row[self.columns["precursor"]]
        quantity_column = self.columns["precursor_quantity"]
        precursor_amount = finite(row[quantity_column]) if quantity_column else None
        if precursor and (not quantity_column or
                          (precursor_amount is not None and precursor_amount > 0)):
            self.precursors[sample].add(precursor)
        protein_column = self.columns["protein"]
        protein = row[protein_column] if protein_column else ""
        protein_quantity_column = self.columns["protein_quantity"]
        protein_amount = finite(row[protein_quantity_column]) if protein_quantity_column else None
        if protein and (not protein_quantity_column or
                        (protein_amount is not None and protein_amount > 0)):
            self.proteins[sample].add(protein)
        if precursor and precursor_amount is not None and precursor_amount > 0:
            self.precursor_values[sample][precursor] = max(
                self.precursor_values[sample].get(precursor, 0),
                precursor_amount,
            )
        if protein and protein_amount is not None and protein_amount > 0:
            self.protein_values[sample][protein] = max(
                self.protein_values[sample].get(protein, 0),
                protein_amount,
            )

    def load(self):
        """Read the aggregated report and collect per-sample values."""
        with open(self.input_path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.select_columns(set(reader.fieldnames or []))
            for row in reader:
                self.record_row(row)
        if not self.observed_samples:
            raise ValueError("DIA-NN report is empty")

    def summarize_sample(self, sample):
        """Write distinct identity counts and log10 quantity distributions."""
        self.metrics.sample_name = sample
        self.metrics.add("quantified_precursor_count", len(self.precursors[sample]))
        if self.columns["protein"]:
            self.metrics.add("quantified_protein_group_count", len(self.proteins[sample]))
        else:
            self.metrics.unavailable("quantified_protein_group_count")
        if self.columns["precursor_quantity"]:
            values = [math.log10(value) for value in self.precursor_values[sample].values()]
            self.metrics.distribution("precursor_quantity_log10", values, "log10 intensity")
        else:
            self.metrics.unavailable_distribution("precursor_quantity_log10", "log10 intensity")
        if self.columns["protein_quantity"]:
            values = [math.log10(value) for value in self.protein_values[sample].values()]
            self.metrics.distribution("protein_group_quantity_log10", values, "log10 intensity")
        else:
            self.metrics.unavailable_distribution(
                name="protein_group_quantity_log10",
                unit="log10 intensity",
            )

    def run(self, output_path):
        """Read the report and write one QC table for all observed samples."""
        self.load()
        for sample in sorted(self.observed_samples):
            self.summarize_sample(sample)
        self.metrics.write(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = DIANNReportSummarizer(
        input_path=args.input,
        sample_csv=args.sample_csv,
        dataset_id=args.dataset_id,
    )
    summarizer.run(args.output)
