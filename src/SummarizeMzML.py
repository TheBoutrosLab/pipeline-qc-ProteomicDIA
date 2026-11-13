#!/usr/bin/env python3
"""Summarize one DIA mzML or Bruker acquisition with pyOpenMS."""
import argparse
from collections import namedtuple
from AcquisitionReader import load_acquisition
from QCMetrics import Metrics, finite

ScanRecord = namedtuple('ScanRecord', ('rt', 'ms_level', 'tic'))


class MzMLSummarizer:
    """Collect acquisition QC from one mzML fraction."""

    def __init__(self, mzml_path, dataset_id, sample_name, fraction_id):
        self.mzml_path = mzml_path
        self.metrics = Metrics(
            dataset_id=dataset_id,
            layer='acquisition',
            source=mzml_path,
            sample_name=sample_name,
            fraction_id=fraction_id,
        )
        self.scans = []

    def load(self):
        """Read spectra with OpenMS and collect the fields needed for QC."""
        experiment = load_acquisition(self.mzml_path)
        for spectrum in experiment.getSpectra():
            self.record_spectrum(spectrum)
        if not self.scans:
            raise ValueError('acquisition has no spectra')

    def record_spectrum(self, spectrum):
        """Record retention time, MS level, and TIC for one spectrum."""
        rt = finite(spectrum.getRT())
        if rt is None or (self.scans and rt < self.scans[-1].rt):
            raise ValueError('acquisition has missing or nonmonotonic retention time')
        _, intensities = spectrum.get_peaks()
        tic = float(sum(intensities)) if len(intensities) else None
        self.scans.append(ScanRecord(rt, spectrum.getMSLevel(), tic))

    def summarize_duration_and_rates(self):
        """Report observed scan counts and rates over recorded duration."""
        duration = self.scans[-1].rt - self.scans[0].rt
        self.metrics.add(
            name='acquisition_duration',
            value=duration / 60,
            unit='min',
        )
        counts = (('scan', len(self.scans)),
                  ('ms1_scan', sum(scan.ms_level == 1 for scan in self.scans)),
                  ('ms2_scan', sum(scan.ms_level == 2 for scan in self.scans)))
        for name, count in counts:
            self.metrics.add(f'{name}_count', count)
            if duration > 0:
                self.metrics.add(
                    name=f'{name}_rate',
                    value=count / duration,
                    unit='scan/s',
                )
            else:
                self.metrics.unavailable(f'{name}_rate', 'scan/s')

    def summarize_cycles(self):
        """Report successive MS1 intervals and their MS2 scan counts."""
        ms1_times = []
        ms2_per_cycle = []
        current_ms2_count = 0
        for scan in self.scans:
            if scan.ms_level == 1:
                if ms1_times:
                    ms2_per_cycle.append(current_ms2_count)
                ms1_times.append(scan.rt)
                current_ms2_count = 0
            elif scan.ms_level == 2:
                current_ms2_count += 1
        intervals = (end - start for start, end in zip(ms1_times, ms1_times[1:]))
        self.metrics.distribution(
            name='ms1_cycle_time',
            values=intervals,
            unit='s',
        )
        self.metrics.distribution('ms2_scans_per_cycle', ms2_per_cycle)

    def summarize_tic(self):
        """Report per-spectrum total ion current."""
        self.metrics.distribution(
            name='total_ion_current',
            values=(scan.tic for scan in self.scans),
            unit='intensity',
        )

    def run(self, output_path):
        """Read the mzML and write its acquisition QC metrics."""
        self.load()
        self.summarize_duration_and_rates()
        self.summarize_cycles()
        self.summarize_tic()
        self.metrics.write(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('input', 'output', 'sample-name', 'fraction-id', 'dataset-id'):
        parser.add_argument('--' + argument, required=True)
    args = parser.parse_args()
    summarizer = MzMLSummarizer(
        mzml_path=args.input,
        dataset_id=args.dataset_id,
        sample_name=args.sample_name,
        fraction_id=args.fraction_id,
    )
    summarizer.run(args.output)
