#!/usr/bin/env python3
"""Render an offline Plotly 3D MS1 stick spectrum from a peak matrix."""

import argparse
import csv
import json
import math
from pathlib import Path

import plotly.graph_objects as go


COLORS = ('#F5C400', '#FF8900', '#F44336', '#B928DB', '#2935C4')


def quantile(sorted_values, fraction):
    position = (len(sorted_values) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (position - lower)


def make_figure(rows, sample_name):
    if not rows:
        raise ValueError('MS1 peak matrix contains no peaks')
    maximum = max(row[2] for row in rows)
    scale_power = math.floor(math.log10(maximum))
    scale = 10 ** scale_power
    log_intensities = sorted(math.log10(row[2]) for row in rows)
    cutoffs = [quantile(log_intensities, q) for q in (0.25, 0.5, 0.75, 0.95)]
    traces = [([], [], []) for _ in COLORS]
    for rt, mz, intensity in rows:
        log_intensity = math.log10(intensity)
        group = next((i for i, cutoff in enumerate(cutoffs) if log_intensity <= cutoff), 4)
        xs, ys, zs = traces[group]
        xs.extend((mz, mz, None))
        ys.extend((rt, rt, None))
        zs.extend((0, intensity / scale, None))
    figure = go.Figure()
    for color, (xs, ys, zs) in zip(COLORS, traces):
        if xs:
            figure.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs, mode='lines',
                line={'color': color, 'width': 2},
                showlegend=False, hoverinfo='skip'))
    figure.update_layout(
        title=f'{sample_name} MS1 peak map',
        scene={
            'xaxis': {'title': 'm/z', 'backgroundcolor': 'white', 'gridcolor': '#D9DDE2'},
            'yaxis': {'title': 'Retention time (s)', 'backgroundcolor': 'white', 'gridcolor': '#D9DDE2'},
            'zaxis': {'title': f'Intensity (×10^{scale_power})', 'backgroundcolor': 'white', 'gridcolor': '#D9DDE2'},
            'camera': {'eye': {'x': 1.8, 'y': -1.7, 'z': 0.75}},
            'aspectmode': 'manual', 'aspectratio': {'x': 2.3, 'y': 1.3, 'z': 0.85},
        },
        paper_bgcolor='white',
        margin={'l': 0, 'r': 0, 'b': 0, 't': 50},
    )
    return figure


def run(args):
    metadata = json.loads(Path(args.metadata).read_text(encoding='utf-8'))
    rows = []
    with open(args.matrix, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'scan_index', 'rt_sec', 'mz', 'intensity'} <= set(reader.fieldnames or []):
            raise ValueError('MS1 peak matrix is missing required columns')
        for row in reader:
            rt, mz, intensity = (float(row['rt_sec']), float(row['mz']), float(row['intensity']))
            if not all(math.isfinite(v) for v in (rt, mz, intensity)) or intensity <= 0:
                raise ValueError('MS1 peak matrix has invalid values')
            rows.append((rt, mz, intensity))
    if len(rows) != metadata['selected_peaks']:
        raise ValueError('MS1 peak matrix row count does not match metadata')
    figure = make_figure(rows, metadata['sample_name'])
    figure.write_html(args.output, include_plotlyjs=True, full_html=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--output', required=True)
    run(parser.parse_args())
