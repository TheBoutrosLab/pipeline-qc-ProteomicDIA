#!/usr/bin/env python3
"""Summarize unique cohort ProtLoom group IDs across sample group TSVs."""

import argparse
from QCMetrics import Metrics, rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-files", nargs="+", required=True, help="Staged per-sample groups TSVs")
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    groups = {}
    for path in args.input_files:
        for row in rows(path, ('group_id', 'proteins')):
            members = tuple(sorted(filter(None, (member.strip() for member in row['proteins'].split(';')))))
            group_id = row['group_id']
            if group_id in groups and groups[group_id] != members:
                raise ValueError(f"ProtLoom group {group_id} has inconsistent members across samples")
            groups[group_id] = members
    metrics = Metrics(args.dataset_id, "", ",".join(args.input_files))
    metrics.add("protein_group_count", len(groups), group="ProtLoom")
    metrics.distribution("proteins_per_group", [len(members) for members in groups.values()], group="ProtLoom")
    metrics.write(args.output)
