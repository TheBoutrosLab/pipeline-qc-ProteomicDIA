#!/usr/bin/env python3
"""Summarize cohort ProteinProphet groups in converted idXML."""

import argparse
from ProteinGroupMetrics import summarize_idxml


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarize_idxml(args.input, args.dataset_id, args.output, "ProteinProphet")
