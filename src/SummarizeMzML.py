#!/usr/bin/env python3
"""Summarize scan timing and DIA windows from one mzML."""

import argparse

import pyopenms as oms

from QCMetrics import Metrics, finite


def run(args):
    metrics = Metrics(args.dataset_id, args.sample_name, args.input)
    times, ms1_times, ms2_times, tic, windows = [], [], [], [], []
    window_lowers, window_uppers, window_widths, window_times = [], [], [], []
    ms2_per_cycle = []
    cycle_ms2 = 0
    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(args.input), experiment)
    for spectrum in experiment.getSpectra():
        level = spectrum.getMSLevel()
        rt = finite(spectrum.getRT())
        if rt is None:
            raise ValueError("Invalid scan start time")
        if times and rt < times[-1]:
            raise ValueError("Nonmonotonic mzML scan times")
        times.append(rt)
        if level == 1:
            if ms1_times:
                ms2_per_cycle.append(cycle_ms2)
            ms1_times.append(rt)
            cycle_ms2 = 0
        elif level == 2:
            ms2_times.append(rt)
            cycle_ms2 += 1
            for precursor in spectrum.getPrecursors():
                center = finite(precursor.getMZ())
                if center is None:
                    continue
                windows.append(center)
                window_times.append((rt, center))
                lo = finite(precursor.getIsolationWindowLowerOffset())
                hi = finite(precursor.getIsolationWindowUpperOffset())
                if lo is not None and hi is not None:
                    window_lowers.append(center - lo)
                    window_uppers.append(center + hi)
                    window_widths.append(lo + hi)
        _, intensities = spectrum.get_peaks()
        if len(intensities):
            tic.append(float(sum(intensities)))
    if not times:
        raise ValueError("mzML contains no spectra")
    duration = times[-1] - times[0]
    metrics.add("acquisition_duration", duration / 60, "min")
    metrics.add("scan_count", len(times))
    metrics.add("ms1_scan_count", len(ms1_times))
    metrics.add("ms2_scan_count", len(ms2_times))
    if duration > 0:
        metrics.add("scan_rate", len(times) / duration, "scan/s")
        metrics.add("ms1_scan_rate", len(ms1_times) / duration, "scan/s")
        metrics.add("ms2_scan_rate", len(ms2_times) / duration, "scan/s")
    metrics.distribution("ms1_cycle_time", [b - a for a, b in zip(ms1_times, ms1_times[1:])], "s")
    metrics.distribution("ms2_scans_per_cycle", ms2_per_cycle)
    if not ms1_times and windows:
        first = windows[0]
        repeated = [rt for rt, center in window_times if center == first]
        metrics.distribution("dia_cycle_time", [b - a for a, b in zip(repeated, repeated[1:])], "s")
    else:
        metrics.distribution("dia_cycle_time", [b - a for a, b in zip(ms1_times, ms1_times[1:])], "s")
    metrics.distribution("total_ion_current", tic)
    metrics.add("distinct_isolation_target_count", len(set(windows)))
    if window_lowers:
        metrics.add("isolation_mz_min", min(window_lowers), "m/z")
        metrics.add("isolation_mz_max", max(window_uppers), "m/z")
        metrics.distribution("isolation_width", window_widths, "m/z")
    metrics.write(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--dataset-id", required=True)
    run(parser.parse_args())
