#!/usr/bin/env python3
"""Summarize retained PSMs and peptides in a final filtered idXML."""

import argparse
from SummarizeIdXML import IdXMLSummarizer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    summarizer = IdXMLSummarizer(
        input_path=args.input,
        dataset_id=args.dataset_id,
        sample_name=args.sample_name,
        search_engine="MSFragger",
        stage="filtered",
    )
    summarizer.run(args.output)
