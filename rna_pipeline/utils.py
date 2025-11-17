"""
utils.py

General-purpose helpers: running shell commands, filesystem helpers,
configuration loading, simple logging.
"""

from pathlib import Path
import subprocess
import sys
from typing import Dict

import yaml


def run_cmd(cmd: list[str], log_file: Path | None = None) -> None:
    """
    Run a shell command and optionally log its output.

    Parameters
    ----------
    cmd : list of str
        Command and arguments.
    log_file : Path, optional
        If provided, stdout/stderr will be written here.

    Raises
    ------
    RuntimeError
        If the command exits with a non-zero status.
    """
    if log_file:
        with log_file.open("w") as lf:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    else:
        proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def ensure_dir(path: Path) -> None:
    """
    Create directory if it does not exist.

    Parameters
    ----------
    path : Path
    """
    path.mkdir(parents=True, exist_ok=True)


def load_reference_config(path: Path | None) -> Dict:
    """
    Load YAML reference configuration file.

    Parameters
    ----------
    path : Path or None

    Returns
    -------
    dict
    """
    if path is None:
        return {}
    with path.open() as f:
        return yaml.safe_load(f)


def subdir(base: Path, name: str) -> Path:
    """
    Create and return a subdirectory of base.

    Parameters
    ----------
    base : Path
    name : str

    Returns
    -------
    Path
    """
    p = base / name
    ensure_dir(p)
    return p
