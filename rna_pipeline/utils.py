"""
utils.py

General-purpose helpers: running shell commands, filesystem helpers,
configuration loading, simple logging.
"""

from pathlib import Path
import subprocess
import sys
import logging
import yaml

def get_logger(name: str) -> logging.Logger:
    """
    Returns a configured logger that prints to the console.
    """
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(levelname)s] %(name)s: %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger

def run_cmd(cmd: list[str], log_file: Path | None = None) -> None:
    """
    Run a shell command and optionally log its output.
    Raises RuntimeError if the command exits with a non-zero status.
    """
    print(f"Running command: {' '.join(cmd)}") # Helpful print to see what's happening
    
    if log_file:
        with log_file.open("w") as lf:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    else:
        proc = subprocess.run(cmd)
        
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")

def ensure_dir(path: Path) -> None:
    """Create directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)

def load_reference_config(path: Path | None) -> dict:
    """Load YAML reference configuration file."""
    if path is None:
        return {}
    
    # Handle string path input safely
    if isinstance(path, str):
        path = Path(path)
        
    with path.open() as f:
        return yaml.safe_load(f)

def subdir(base: Path, name: str) -> Path:
    """Create and return a subdirectory of base."""
    p = base / name
    ensure_dir(p)
    return p