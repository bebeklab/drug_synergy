"""Repository-relative paths used throughout the analysis."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "analysis.json"
DATA = ROOT / "data"
INPUTS = DATA / "inputs"
ANNOTATIONS = DATA / "annotations"
RESULTS = ROOT / "results" / "generated"
