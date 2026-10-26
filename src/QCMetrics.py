"""Shared OpenMS readers and compact QC metric output."""

import csv
import math
from collections import Counter

import pyopenms as oms


FIELDS = ["dataset_id", "sample_name", "metric_name", "group", "value", "unit", "source_file"]
QUANTILES = (5, 25, 50, 75, 95)


def as_text(value):
    """Convert an OpenMS byte string or value to text."""
    return value.decode() if isinstance(value, bytes) else str(value)


def finite(value):
    """Return a finite numeric value, or None for invalid input."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def percentile(sorted_values, percent):
    """Interpolate a percentile from sorted numeric values."""
    if not sorted_values:
        return None
    position = (len(sorted_values) - 1) * percent / 100
    left = math.floor(position)
    right = math.ceil(position)
    return sorted_values[left] + (sorted_values[right] - sorted_values[left]) * (position - left)


class Metrics:
    """Collect and write long-form QC metric rows."""
    def __init__(self, dataset_id, sample_name, source):
        self.dataset_id = dataset_id
        self.sample_name = sample_name or ""
        self.source = str(source)
        self.rows = []

    def add(self, name, value, unit="count", group="all"):
        """Append one metric row when its value is available."""
        if value is None:
            return
        self.rows.append({
            "dataset_id": self.dataset_id, "sample_name": self.sample_name,
            "metric_name": name, "group": group, "value": value, "unit": unit,
            "source_file": self.source,
        })

    def distribution(self, name, values, unit="", group="all"):
        """Record the count and selected percentiles of a distribution."""
        numbers = sorted(v for value in values if (v := finite(value)) is not None)
        self.add(name + "_count", len(numbers), group=group)
        for q in QUANTILES:
            self.add(f"{name}_p{q:02d}", percentile(numbers, q), unit, group)

    def write(self, path):
        """Write metric rows as a tab-separated table."""
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)


def rows(path, required=()):
    """Yield TSV rows after checking required column names."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        yield from reader


def score_summary(metrics, path, name="search_score"):
    """Summarize the primary score of each best-scoring peptide hit."""
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(path), protein_ids, peptide_ids)
    groups = {"all": [], "target": [], "decoy": [], "unknown": []}
    candidates = 0
    peptides = set()
    charges = Counter()
    score_types = set()
    score_directions = set()
    for identification in peptide_ids:
        hits = identification.getHits()
        candidates += len(hits)
        if not hits:
            continue
        higher = identification.isHigherScoreBetter()
        score_types.add(as_text(identification.getScoreType()))
        score_directions.add(higher)
        scored = [(finite(hit.getScore()), hit) for hit in hits]
        valid = [(score, hit) for score, hit in scored if score is not None]
        if not valid:
            continue
        score, rank_one = (max(valid, key=lambda pair: pair[0]) if higher
                           else min(valid, key=lambda pair: pair[0]))
        groups["all"].append(score)
        label = (as_text(rank_one.getMetaValue("target_decoy")).lower()
                 if rank_one.metaValueExists("target_decoy") else "unknown")
        if label not in ("target", "decoy"):
            label = "unknown"
        groups[label].append(score)
        peptides.add(as_text(rank_one.getSequence().toString()))
        charges[rank_one.getCharge()] += 1
    metrics.add("peptide_identification_count", len(peptide_ids))
    metrics.add("candidate_psm_count", candidates)
    metrics.add("rank1_hit_count", len(groups["all"]))
    metrics.add("distinct_rank1_peptide_count", len(peptides - {""}))
    if len(score_types) == 1:
        metrics.add(name + "_type", next(iter(score_types)), unit="")
    if len(score_directions) == 1:
        metrics.add(name + "_higher_is_better", int(next(iter(score_directions))), unit="boolean")
    metrics.add("target_decoy_label_missing_count", len(groups["unknown"]))
    for label, values in groups.items():
        metrics.distribution(name, values, group=label)
    for charge, count in sorted(charges.items()):
        metrics.add("rank1_charge_count", count, group=charge)
