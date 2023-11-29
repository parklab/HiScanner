from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("scanner")
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"

