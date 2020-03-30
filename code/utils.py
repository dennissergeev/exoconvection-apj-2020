"""Generic utilities."""
from datetime import datetime
import sys

from loguru import logger


def pprint_dict(d):
    """Print each dictionary key-value pair on a new line."""
    return "\n".join(f"{k} = {v}" for k, v in d.items())


def create_logger(script_path, subdir="logs"):
    """Create a logger using loguru."""
    logpath = script_path.parent / subdir
    logpath.mkdir(exist_ok=True)
    logger.configure(
        handlers=[
            {"sink": sys.stdout, "level": "INFO"},
            {
                "sink": logpath
                / f"log_{script_path.stem}_{datetime.now():%Y-%m-%d_%H%M%S}.log",
                "level": "DEBUG",
            },
        ],
    )
    return logger


def tex2cf_units(unit_str):
    """Convert a TeX string to a string that can be used in cf_units."""
    return (
        unit_str.replace("$", "").replace("{", "").replace("}", "").replace("^", "**")
    )


def lsdir(path):
    """List all files and directories in the given path, sorted."""
    return sorted(path.glob("*"))
