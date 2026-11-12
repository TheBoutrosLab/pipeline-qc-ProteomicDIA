#!/usr/bin/env python3
"""Summarize DIA-NN reports by acquisition and combined sample."""

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
    """Map unambiguous DIA-NN run labels to sample and fraction pairs."""
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
                    if key in lookup and lookup[key] != pair:
                        raise ValueError(f"Ambiguous DIA-NN run key: {key}")
                    lookup[key] = pair
    counts = defaultdict(int)
    for sample, _ in pairs:
        counts[sample] += 1
    for sample, fraction in pairs:
        target = (sample, fraction) if counts[sample] == 1 else (sample, "")
        if sample in lookup and lookup[sample] != target:
            raise ValueError(f"Ambiguous DIA-NN run key: {sample}")
        lookup[sample] = target
    return lookup, pairs


def summarize(args):
    """Write acquisition, sample, and cohort quantification metrics."""
    lookup, _ = run_lookup(args.sample_csv)
    precursors = defaultdict(set)
    proteins = defaultdict(set)
    precursor_values = defaultdict(dict)
    protein_values = defaultdict(dict)
    run_rows = defaultdict(int)
    run_names = {}
    has_protein = False
    has_precursor_quantity = False
    has_protein_quantity = False
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
            raise ValueError("DIA-NN report needs Run and Precursor.Id")
        has_protein |= bool(protein_col)
        has_precursor_quantity |= bool(precursor_quantity)
        has_protein_quantity |= bool(protein_quantity)
        for row in reader:
            run_name = row[run_col]
            candidates = [run_name, Path(run_name).name, Path(run_name).stem]
            matches = {lookup[key] for key in candidates if key in lookup}
            if len(matches) != 1:
                raise ValueError(f"Unknown or ambiguous DIA-NN run: {run_name}")
            pair = matches.pop()
            if run_name in run_names and run_names[run_name] != pair:
                raise ValueError(f"DIA-NN run maps to multiple acquisitions: {run_name}")
            run_names[run_name] = pair
            run_rows[pair] += 1
            precursor = row[precursor_col]
            precursor_amount = finite(row[precursor_quantity]) if precursor_quantity else None
            if precursor and (not precursor_quantity or
                              (precursor_amount is not None and precursor_amount > 0)):
                precursors[pair].add(precursor)
            protein = row[protein_col] if protein_col else ""
            protein_amount = finite(row[protein_quantity]) if protein_quantity else None
            if protein and (not protein_quantity or
                            (protein_amount is not None and protein_amount > 0)):
                proteins[pair].add(protein)
            if precursor and precursor_amount is not None and precursor_amount > 0:
                precursor_values[pair][precursor] = max(
                    precursor_values[pair].get(precursor, 0), precursor_amount)
            if protein and protein_amount is not None and protein_amount > 0:
                protein_values[pair][protein] = max(
                    protein_values[pair].get(protein, 0), protein_amount)
    if not run_rows:
        raise ValueError("DIA-NN report is empty")
    cohort_precursors = set().union(*precursors.values())
    metrics = Metrics(args.dataset_id, "", args.input, "quantification")
    metrics.add("diann_run_count", len(run_rows))
    metrics.add("cohort_precursor_count", len(cohort_precursors))

    def write_group(sample, fraction, row_count, precursor_set, protein_set,
                    precursor_amounts, protein_amounts):
        metrics.sample_name = sample
        metrics.context["fraction_id"] = fraction
        metrics.add("report_row_count", row_count)
        metrics.add("quantified_precursor_count", len(precursor_set))
        missing = len(cohort_precursors - precursor_set)
        metrics.add("precursor_missing_from_cohort_union", missing)
        if cohort_precursors:
            metrics.add("precursor_missing_fraction", missing / len(cohort_precursors), "fraction")
        if has_protein:
            metrics.add("quantified_protein_group_count", len(protein_set))
        if has_precursor_quantity:
            values = [math.log10(v) for v in precursor_amounts.values()]
            metrics.distribution("precursor_quantity_log10", values, "log10 intensity")
        if has_protein_quantity:
            values = [math.log10(v) for v in protein_amounts.values()]
            metrics.distribution("protein_quantity_log10", values, "log10 intensity")

    observed_pairs = sorted(run_rows)
    for sample, fraction in observed_pairs:
        pair = sample, fraction
        write_group(sample, fraction, run_rows[pair], precursors[pair], proteins[pair],
                    precursor_values[pair], protein_values[pair])
    samples = sorted({sample for sample, _ in observed_pairs})
    for sample in samples:
        sample_pairs = [pair for pair in observed_pairs if pair[0] == sample]
        if len(sample_pairs) == 1:
            continue
        combined_precursors = set().union(*(precursors[pair] for pair in sample_pairs))
        combined_proteins = set().union(*(proteins[pair] for pair in sample_pairs))
        combined_precursor_values = {}
        combined_protein_values = {}
        for pair in sample_pairs:
            for key, value in precursor_values[pair].items():
                combined_precursor_values[key] = combined_precursor_values.get(key, 0) + value
            for key, value in protein_values[pair].items():
                combined_protein_values[key] = combined_protein_values.get(key, 0) + value
        write_group(sample, "", sum(run_rows[pair] for pair in sample_pairs),
                    combined_precursors, combined_proteins,
                    combined_precursor_values, combined_protein_values)
    metrics.write(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    summarize(parser.parse_args())
