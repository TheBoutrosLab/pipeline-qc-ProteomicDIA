"""Shared OpenMS readers and compact QC metric output."""

import csv
import math


FIELDS = ["dataset_id", "layer", "source_process", "sample_name", "search_engine",
          "fraction_id", "metric_name", "group", "value", "unit", "status", "source_file"]
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
    def __init__(self, dataset_id, sample_name, source, layer, **dimensions):
        self.context = {
            "dataset_id": dataset_id, "layer": layer,
            "source_process": dimensions.get("source_process", ""),
            "search_engine": dimensions.get("search_engine", ""),
            "fraction_id": dimensions.get("fraction_id", ""),
        }
        self.sample_name = sample_name or ""
        self.source = str(source)
        self.rows = []

    def add(self, name, value, unit="count", group="all", status="available"):
        """Append a metric row, including an explicit unavailable status."""
        if value is None and status == "available":
            return
        self.rows.append({
            **self.context, "sample_name": self.sample_name,
            "metric_name": name, "group": group,
            "value": "" if value is None else value, "unit": unit,
            "status": status, "source_file": self.source,
        })

    def unavailable(self, name, unit="", group="all"):
        """Record a missing source feature without treating it as zero."""
        self.add(name, None, unit, group, "unavailable")

    def distribution(self, name, values, unit="", group="all"):
        """Record the count and selected percentiles of a distribution."""
        values = list(values)
        numbers = sorted(v for value in values if (v := finite(value)) is not None)
        self.add(name + "_count", len(numbers), group=group)
        self.add(name + "_missing_count", len(values) - len(numbers), group=group)
        for q in QUANTILES:
            if numbers:
                self.add(f"{name}_p{q:02d}", percentile(numbers, q), unit, group)
            else:
                self.unavailable(f"{name}_p{q:02d}", unit, group)

    def unavailable_distribution(self, name, unit="", group="all"):
        """Mark all seven rows unavailable when a required source field is absent."""
        for suffix in ("count", "missing_count", "p05", "p25", "p50", "p75", "p95"):
            self.unavailable(
                name=f"{name}_{suffix}",
                unit="" if suffix.endswith("count") else unit,
                group=group,
            )

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
