#!/usr/bin/env python3
"""Summarize cohort Epifany groups in idXML."""

import argparse
from ProteinGroupMetrics import ProteinGroupSummarizer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = ProteinGroupSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        route="Epifany",
    )
    summarizer.run(args.output)
