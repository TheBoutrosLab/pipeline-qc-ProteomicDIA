#!/usr/bin/env python3
"""Summarize retained PSMs and peptides in a final filtered idXML."""

import argparse
from QCMetrics import Metrics, score_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    metrics = Metrics(args.dataset_id, args.sample_name, args.input, "validated_identification")
    score_summary(metrics, args.input, name="validated_score")
    metrics.write(args.output)
