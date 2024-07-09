"""Data test."""
import os
import glob

from linkml_runtime.loaders import yaml_loader
from gocam.datamodel import Model

ROOT = os.path.join(os.path.dirname(__file__), '..')
DATA_DIR = os.path.join(ROOT, "src", "data", "examples")

EXAMPLE_FILES = glob.glob(os.path.join(DATA_DIR, '*.yaml'))


def test_data():
    """Data test."""
    for path in EXAMPLE_FILES:
        obj = yaml_loader.load(path, target_class=Model)
        assert obj
