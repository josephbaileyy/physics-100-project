#!/usr/bin/env python3
"""Compatibility entrypoint for AM CVn sequence photometry."""
from pathlib import Path
import runpy
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

runpy.run_module("amcvn.sequence_analysis", run_name="__main__")
