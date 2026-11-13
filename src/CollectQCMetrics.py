#!/usr/bin/env python3
"""Combine DIA QC metric TSVs into one deterministic table."""
import argparse
import csv
from QCMetrics import FIELDS

LAYER_ORDER = ('acquisition', 'raw_identification', 'psm_features',
               'validated_identification', 'protein_inference', 'spectral_library',
               'quantification')
LAYER_RANK = {layer: rank for rank, layer in enumerate(LAYER_ORDER)}


class QCMetricsCollector:
    """Combine stage metrics into one validated, ordered QC table."""

    def __init__(self, input_paths):
        self.input_paths = input_paths
        self.rows = []

    def read_inputs(self):
        """Read each stage table and require the shared column schema."""
        for path in self.input_paths:
            with open(path, newline='', encoding='utf-8') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if reader.fieldnames != FIELDS:
                    raise ValueError(f'Unexpected QC metric columns in {path}: {reader.fieldnames}')
                self.rows.extend(reader)

    def validate_layers(self):
        """Reject layer names that cannot be placed in pipeline order."""
        unknown = {row['layer'] for row in self.rows} - LAYER_RANK.keys()
        if unknown:
            raise ValueError(f"Unknown metric layers: {', '.join(sorted(unknown))}")

    @staticmethod
    def sort_key(row):
        """Return the stable cohort, layer, assay, and metric order."""
        return (
            row['dataset_id'],
            LAYER_RANK[row['layer']],
            row['sample_name'],
            row['search_engine'],
            row['fraction_id'],
            row['metric_name'],
            row['group'],
            row['source_file'],
        )

    def sort_rows(self):
        """Order rows by cohort, pipeline layer, assay, and metric."""
        self.rows.sort(key=self.sort_key)

    def write(self, output_path):
        """Write the merged long-form TSV."""
        with open(output_path, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows)

    def run(self, output_path):
        """Collect, validate, sort, and write metrics."""
        self.read_inputs()
        self.validate_layers()
        self.sort_rows()
        self.write(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    QCMetricsCollector(args.inputs).run(args.output)
