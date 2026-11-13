#!/usr/bin/env python3
"""Extract a bounded, reproducible MS1 peak matrix for a 3D stick spectrum."""

import argparse
import csv
import heapq
import json
import math
from pathlib import Path

import pyopenms as oms

from AcquisitionReader import is_bruker, load_acquisition


def evenly_spaced_indices(count, limit):
    """Choose MS1 scans across the complete run, including both endpoints."""
    if count <= limit:
        return list(range(count))
    if limit == 1:
        return [count // 2]
    return [round(i * (count - 1) / (limit - 1)) for i in range(limit)]


def ms1_spectra(path, max_scans):
    """Select evenly spaced MS1 spectra from mzML or Bruker .d input."""
    if is_bruker(path):
        experiment = load_acquisition(path, ms1_only=True)
        spectra = [spectrum for spectrum in experiment.getSpectra()
                   if spectrum.getMSLevel() == 1]
        if not spectra:
            raise ValueError('Bruker .d acquisition has no MS1 spectra')
        selected = evenly_spaced_indices(len(spectra), max_scans)
        return (spectra[i] for i in selected), len(spectra)
    ondisk = oms.OnDiscMSExperiment()
    if ondisk.openFile(str(path)):
        metadata = ondisk.getMetaData()
        indices = [i for i in range(ondisk.getNrSpectra())
                   if metadata.getSpectrum(i).getMSLevel() == 1]
        if not indices:
            raise ValueError('mzML has no MS1 spectra')
        selected = [indices[i] for i in evenly_spaced_indices(len(indices), max_scans)]
        return (ondisk.getSpectrum(i) for i in selected), len(indices)

    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(path), experiment)
    spectra = [s for s in experiment.getSpectra() if s.getMSLevel() == 1]
    if not spectra:
        raise ValueError('mzML has no MS1 spectra')
    selected = evenly_spaced_indices(len(spectra), max_scans)
    return (spectra[i] for i in selected), len(spectra)


def extract_peaks(spectra, mz_min, mz_max, peaks_per_scan):
    """Retain the strongest valid peaks from each selected spectrum."""
    if not (math.isfinite(mz_min) and math.isfinite(mz_max) and mz_min < mz_max and
            peaks_per_scan > 0):
        raise ValueError('Invalid MS1 peak-selection limits')
    rows = []
    selected_scan_count = 0
    for spectrum in spectra:
        rt_sec = float(spectrum.getRT())
        if not math.isfinite(rt_sec):
            continue
        mz_values, intensities = spectrum.get_peaks()
        valid = ((float(intensity), float(mz)) for mz, intensity
                 in zip(mz_values, intensities)
                 if math.isfinite(float(mz)) and math.isfinite(float(intensity))
                 and mz_min <= float(mz) <= mz_max and float(intensity) > 0)
        strongest = heapq.nlargest(peaks_per_scan, valid)
        for intensity, mz in strongest:
            rows.append((selected_scan_count, rt_sec, mz, intensity))
        selected_scan_count += 1
    if not rows:
        raise ValueError('No positive MS1 peaks within the selected m/z range')
    return rows, selected_scan_count


class MS1MatrixExtractor:
    """Select MS1 peaks and publish a matrix with its selection settings."""

    def __init__(
        self,
        input_path,
        sample_name,
        fraction_id,
        mz_min,
        mz_max,
        max_scans,
        peaks_per_scan,
    ):
        self.input_path = input_path
        self.sample_name = sample_name
        self.fraction_id = fraction_id
        self.mz_min = mz_min
        self.mz_max = mz_max
        self.max_scans = max_scans
        self.peaks_per_scan = peaks_per_scan

    def write_matrix(self, path, rows):
        """Write scan index, retention time, m/z, and intensity."""
        with open(path, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['scan_index', 'rt_sec', 'mz', 'intensity'])
            for index, rt, mz, intensity in rows:
                writer.writerow([index, f'{rt:.6f}', f'{mz:.6f}', f'{intensity:.9g}'])

    def write_metadata(self, path, rows, total_scans, selected_scans):
        """Record acquisition coverage and peak-selection limits."""
        metadata = {
            'sample_name': self.sample_name,
            'fraction_id': self.fraction_id,
            'total_ms1_scans': total_scans,
            'selected_ms1_scans': selected_scans,
            'selected_peaks': len(rows),
            'rt_min_sec': min(row[1] for row in rows),
            'rt_max_sec': max(row[1] for row in rows),
            'mz_min': self.mz_min,
            'mz_max': self.mz_max,
            'max_scans': self.max_scans,
            'peaks_per_scan': self.peaks_per_scan,
            'selection': ('Evenly spaced MS1 scans across the run; strongest positive peaks '
                          'per selected scan within the m/z range'),
        }
        Path(path).write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')

    def run(self, matrix_path, metadata_path):
        """Load one acquisition, select its peaks, and write both files."""
        if self.max_scans < 1:
            raise ValueError('max-scans must be positive')
        spectra, total_scans = ms1_spectra(self.input_path, self.max_scans)
        rows, selected_scans = extract_peaks(
            spectra=spectra,
            mz_min=self.mz_min,
            mz_max=self.mz_max,
            peaks_per_scan=self.peaks_per_scan,
        )
        self.write_matrix(matrix_path, rows)
        self.write_metadata(
            path=metadata_path,
            rows=rows,
            total_scans=total_scans,
            selected_scans=selected_scans,
        )


def main():
    """Parse arguments and extract one MS1 peak matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--sample-name', required=True)
    parser.add_argument('--fraction-id', required=True)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--mz-min', type=float, default=300.0)
    parser.add_argument('--mz-max', type=float, default=1200.0)
    parser.add_argument('--max-scans', type=int, default=350)
    parser.add_argument('--peaks-per-scan', type=int, default=120)
    args = parser.parse_args()
    extractor = MS1MatrixExtractor(
        input_path=args.input,
        sample_name=args.sample_name,
        fraction_id=args.fraction_id,
        mz_min=args.mz_min,
        mz_max=args.mz_max,
        max_scans=args.max_scans,
        peaks_per_scan=args.peaks_per_scan,
    )
    extractor.run(args.matrix, args.metadata)


if __name__ == '__main__':
    main()
