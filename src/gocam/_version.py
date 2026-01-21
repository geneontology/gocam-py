from importlib.metadata import PackageNotFoundError, version

try:
    assert __package__ is not None
    __version__ = version(__package__)
except (PackageNotFoundError, AssertionError):
    # package not installed
    __version__ = "0.0.0"
