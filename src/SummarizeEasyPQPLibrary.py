#!/usr/bin/env python3
"""Summarize EasyPQP library rows and the validation export summary."""

import argparse
import csv
from QCMetrics import Metrics, rows


class EasyPQPLibrarySummarizer:
    """Summarize cohort library size and its export diagnostics."""

    def __init__(self, library_path, export_summary_path, dataset_id):
        self.library_path = library_path
        self.export_summary_path = export_summary_path
        self.metrics = Metrics(dataset_id, "", library_path, "spectral_library")

    def summarize_library(self):
        """Count transitions and distinct peptide and precursor identifiers."""
        with open(self.library_path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = set(reader.fieldnames or [])
            peptide_names = ("ModifiedPeptideSequence", "PeptideSequence", "Peptide")
            peptide_column = next((name for name in peptide_names if name in columns), None)
            precursor_names = ("TransitionGroupId", "PrecursorId", "Precursor.Id")
            precursor_column = next((name for name in precursor_names if name in columns), None)
            if not peptide_column:
                raise ValueError("EasyPQP library is missing a peptide sequence column")
            if not precursor_column and "PrecursorCharge" not in columns:
                raise ValueError("EasyPQP library needs TransitionGroupId or PrecursorCharge")
            peptide_ids, precursor_ids = set(), set()
            transition_count = 0
            for row in reader:
                transition_count += 1
                peptide_ids.add(row[peptide_column])
                precursor_id = (row[precursor_column] if precursor_column else
                                (row[peptide_column], row["PrecursorCharge"]))
                precursor_ids.add(precursor_id)
        self.metrics.add("library_transition_count", transition_count)
        self.metrics.add("library_peptide_count", len(peptide_ids - {""}))
        self.metrics.add("library_precursor_count", len(precursor_ids - {""}))

    def summarize_export(self):
        """Record the source metric names from the export summary."""
        for row in rows(self.export_summary_path, ("metric", "count")):
            self.metrics.add("library_export_" + row["metric"], int(row["count"]))

    def run(self, output_path):
        """Read and write cohort spectral-library QC."""
        self.summarize_library()
        self.summarize_export()
        self.metrics.write(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True)
    parser.add_argument("--export-summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = EasyPQPLibrarySummarizer(
        library_path=args.library,
        export_summary_path=args.export_summary,
        dataset_id=args.dataset_id,
    )
    summarizer.run(args.output)
