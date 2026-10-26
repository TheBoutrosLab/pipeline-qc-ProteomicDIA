#!/usr/bin/env python3
"""Summarize cohort Fido groups in the parsed group JSON."""

import argparse
import json
from QCMetrics import Metrics, finite


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    with open(args.input, encoding="utf-8") as handle:
        groups = json.load(handle)
    if not isinstance(groups, list):
        raise ValueError("Fido groups JSON must be a list")
    if any(not isinstance(group.get("accessions"), list) or finite(group.get("probability")) is None
           for group in groups):
        raise ValueError("Fido group is missing accessions or probability")
    metrics = Metrics(args.dataset_id, "", args.input)
    metrics.add("protein_group_count", len(groups), group="Fido")
    metrics.distribution("proteins_per_group", [len(group["accessions"]) for group in groups], group="Fido")
    metrics.distribution("protein_group_probability", [group["probability"] for group in groups], group="Fido")
    metrics.write(args.output)
