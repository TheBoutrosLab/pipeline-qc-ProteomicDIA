#!/usr/bin/env python3
"""Summarize EasyPQP library rows and the validation export summary."""

import argparse
import csv
from QCMetrics import Metrics, rows


def run(args):
    metrics = Metrics(args.dataset_id, "", args.library)
    with open(args.library, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        peptide_column = next((name for name in ("ModifiedPeptideSequence", "PeptideSequence", "Peptide")
                               if name in columns), None)
        precursor_column = next((name for name in ("TransitionGroupId", "PrecursorId", "Precursor.Id")
                                 if name in columns), None)
        if not peptide_column:
            raise ValueError("EasyPQP library is missing a peptide sequence column")
        if not precursor_column and "PrecursorCharge" not in columns:
            raise ValueError("EasyPQP library needs TransitionGroupId or PrecursorCharge")
        peptide_ids, precursor_ids = set(), set()
        transition_count = 0
        for row in reader:
            transition_count += 1
            if peptide_column:
                peptide_ids.add(row[peptide_column])
            precursor_ids.add(row[precursor_column] if precursor_column else (row[peptide_column], row['PrecursorCharge']))
    metrics.add("library_transition_count", transition_count)
    metrics.add("library_peptide_count", len(peptide_ids - {""}))
    metrics.add("library_precursor_count", len(precursor_ids - {""}))
    for row in rows(args.export_summary, ("metric", "count")):
        metrics.add("library_export_" + row["metric"], int(row["count"]))
    metrics.write(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True)
    parser.add_argument("--export-summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    run(parser.parse_args())
