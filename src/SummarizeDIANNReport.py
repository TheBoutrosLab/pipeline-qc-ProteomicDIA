#!/usr/bin/env python3
"""Summarize per-run precursor and protein quantities in a DIA-NN report TSV."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from QCMetrics import Metrics, finite


def first(columns, candidates):
    """Return the first available column name from a list of candidates."""
    return next((name for name in candidates if name in columns), None)


def summarize(args):
    """Summarize per-sample quantities in the DIA-NN report."""
    sample_lookup = {}
    with open(args.sample_csv, newline="", encoding="utf-8") as handle:
        samples = csv.DictReader(handle)
        columns = set(samples.fieldnames or [])
        mzml_column = 'path' if 'path' in columns else 'mzML'
        if 'sample_name' not in columns or mzml_column not in columns:
            raise ValueError("mzML manifest needs sample_name and path (or mzML)")
        for row in samples:
            name = row['sample_name']
            mzml = Path(row[mzml_column])
            for key in (name, mzml.name, mzml.stem):
                if key in sample_lookup and sample_lookup[key] != name:
                    raise ValueError(f"Ambiguous DIA-NN run key: {key}")
                sample_lookup[key] = name
    with open(args.input, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        run_col = first(columns, ("Run", "File.Name"))
        precursor_col = first(columns, ("Precursor.Id", "Precursor.ID"))
        precursor_quantity = first(columns, (
            "Precursor.Quantity", "Precursor.Normalised", "Precursor.Translated"))
        protein_col = first(columns, ("Protein.Group", "Protein.Ids"))
        protein_quantity = first(columns, ("PG.Quantity", "PG.Normalised"))
        if not run_col or not precursor_col:
            raise ValueError("DIA-NN report needs Run and Precursor.Id columns")
        precursors = defaultdict(set)
        proteins = defaultdict(set)
        precursor_values = defaultdict(dict)
        protein_values = defaultdict(dict)
        run_rows = defaultdict(int)
        for row in reader:
            run_name = row[run_col]
            run = (sample_lookup.get(run_name)
                   or sample_lookup.get(Path(run_name).name)
                   or sample_lookup.get(Path(run_name).stem))
            if not run:
                raise ValueError(f"DIA-NN run is absent from the mzML manifest: {run_name}")
            run_rows[run] += 1
            precursor = row[precursor_col]
            precursor_amount = finite(row[precursor_quantity]) if precursor_quantity else None
            if precursor and (not precursor_quantity or
                              (precursor_amount is not None and precursor_amount > 0)):
                precursors[run].add(precursor)
            protein_amount = finite(row[protein_quantity]) if protein_quantity else None
            if (protein_col and row[protein_col] and
                    (not protein_quantity or
                     (protein_amount is not None and protein_amount > 0))):
                proteins[run].add(row[protein_col])
            if precursor and precursor_amount is not None and precursor_amount > 0:
                precursor_values[run][precursor] = max(
                    precursor_values[run].get(precursor, 0), precursor_amount)
            if (protein_col and row[protein_col] and protein_amount is not None
                    and protein_amount > 0):
                protein_values[run][row[protein_col]] = max(
                    protein_values[run].get(row[protein_col], 0), protein_amount)
    if not run_rows:
        raise ValueError("DIA-NN report is empty")
    cohort_precursors = set().union(*precursors.values())
    metrics = Metrics(args.dataset_id, "", args.input, "quantification")
    metrics.add("diann_run_count", len(run_rows))
    metrics.add("cohort_precursor_count", len(cohort_precursors))
    for run in sorted(run_rows):
        metrics.sample_name = run
        metrics.add("report_row_count", run_rows[run])
        metrics.add("quantified_precursor_count", len(precursors[run]))
        metrics.add("precursor_missing_from_cohort_union", len(cohort_precursors - precursors[run]))
        if cohort_precursors:
            missing_fraction = len(cohort_precursors - precursors[run]) / len(cohort_precursors)
            metrics.add("precursor_missing_fraction", missing_fraction, "fraction")
        if protein_col:
            metrics.add("quantified_protein_group_count", len(proteins[run]))
        if precursor_quantity:
            values = [math.log10(value) for value in precursor_values[run].values()]
            metrics.distribution("precursor_quantity_log10", values, "log10 intensity")
        if protein_quantity:
            values = [math.log10(value) for value in protein_values[run].values()]
            metrics.distribution("protein_quantity_log10", values, "log10 intensity")
    metrics.write(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    summarize(parser.parse_args())
