#!/usr/bin/env python3
"""Combine stage metric TSVs into one deterministically ordered table."""

import argparse
import csv
from QCMetrics import FIELDS

LAYER_ORDER = (
    "acquisition", "raw_identification", "psm_features",
    "validated_identification", "protein_inference", "spectral_library",
    "quantification",
)
LAYER_RANK = {layer: rank for rank, layer in enumerate(LAYER_ORDER)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    combined = []
    for path in args.inputs:
        with open(path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise ValueError(f"Unexpected metric columns in {path}: {reader.fieldnames}")
            combined.extend(reader)
    unknown_layers = {row["layer"] for row in combined} - LAYER_RANK.keys()
    if unknown_layers:
        raise ValueError(f"Unknown metric layers: {', '.join(sorted(unknown_layers))}")
    combined.sort(key=lambda row: (
        row["dataset_id"], LAYER_RANK[row["layer"]], row["sample_name"],
        row["search_engine"], row["fraction_id"], row["metric_name"],
        row["group"], row["source_file"]))
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(combined)
