#!/usr/bin/env python3
"""Summarize unique cohort ProtLoom group IDs across sample group TSVs."""

import argparse
from QCMetrics import Metrics, rows


class ProtLoomGroupSummarizer:
    """Count unique cohort ProtLoom groups across sample tables."""

    def __init__(self, input_paths, dataset_id):
        self.input_paths = input_paths
        self.metrics = Metrics(
            dataset_id=dataset_id,
            sample_name="",
            source=",".join(input_paths),
            layer="protein_inference",
        )
        self.groups = {}

    def load(self):
        """Read group memberships and reject conflicting group IDs."""
        for path in self.input_paths:
            for row in rows(path, ("group_id", "proteins")):
                members = tuple(sorted(
                    member.strip() for member in row["proteins"].split(";") if member.strip()
                ))
                group_id = row["group_id"]
                if group_id in self.groups and self.groups[group_id] != members:
                    raise ValueError(
                        f"ProtLoom group {group_id} has inconsistent members across samples"
                    )
                self.groups[group_id] = members

    def summarize(self):
        """Report cohort group count and membership distribution."""
        self.metrics.add("protein_group_count", len(self.groups))
        self.metrics.distribution(
            name="proteins_per_group",
            values=[len(members) for members in self.groups.values()],
        )

    def run(self, output_path):
        """Read, summarize, and write cohort ProtLoom QC."""
        self.load()
        self.summarize()
        self.metrics.write(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-files", nargs="+", required=True,
                        help="Staged per-sample groups TSVs")
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = ProtLoomGroupSummarizer(
        input_paths=args.input_files,
        dataset_id=args.dataset_id,
    )
    summarizer.run(args.output)
