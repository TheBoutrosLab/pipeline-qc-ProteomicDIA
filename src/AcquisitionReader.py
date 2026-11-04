"""Load mzML files or Bruker timsTOF directories with pyOpenMS."""

from pathlib import Path

import pyopenms as oms


def is_bruker(path):
    """Identify a Bruker acquisition directory from its path."""
    source = Path(path)
    return source.is_dir() and source.name.lower().endswith(".d")


def load_acquisition(path, ms1_only=False):
    """Return an MSExperiment while retaining raw Bruker MS1 frames."""
    if is_bruker(path):
        if not hasattr(oms, "BrukerTimsFile"):
            raise RuntimeError("pyOpenMS needs BrukerTimsFile support for .d input")
        reader = oms.BrukerTimsFile()
        config = oms.BrukerTimsFile.Config()
        config.load_ms1 = True
        if ms1_only and hasattr(config, "load_ms2"):
            config.load_ms2 = False
        return reader.load(str(path), config)
    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(path), experiment)
    return experiment
