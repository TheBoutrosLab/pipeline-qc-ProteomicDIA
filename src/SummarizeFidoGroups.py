#!/usr/bin/env python3
"""Summarize cohort Fido groups in the parsed group JSON."""

import argparse
import json
from QCMetrics import Metrics, finite


class FidoGroupSummarizer:
    """Read cohort Fido groups and report group size and probability."""

    def __init__(self, input_path, dataset_id):
        self.input_path = input_path
        self.metrics = Metrics(dataset_id, "", input_path, "protein_inference")
        self.groups = []

    def load(self):
        """Read the JSON group list and validate required fields."""
        with open(self.input_path, encoding="utf-8") as handle:
            self.groups = json.load(handle)
        if not isinstance(self.groups, list):
            raise ValueError("Fido groups JSON must be a list")
        if any(not isinstance(group.get("accessions"), list) or
               finite(group.get("probability")) is None for group in self.groups):
            raise ValueError("Fido group is missing accessions or probability")

    def summarize(self):
        """Count groups and summarize their sizes and probabilities."""
        self.metrics.add("protein_group_count", len(self.groups))
        self.metrics.distribution(
            name="proteins_per_group",
            values=[len(group["accessions"]) for group in self.groups],
        )
        self.metrics.distribution(
            name="protein_group_probability",
            values=[group["probability"] for group in self.groups],
        )

    def run(self, output_path):
        """Read, summarize, and write cohort Fido QC."""
        self.load()
        self.summarize()
        self.metrics.write(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = FidoGroupSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
    )
    summarizer.run(args.output)
