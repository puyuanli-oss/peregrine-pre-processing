"""Environment and dependency checks for the PEREGRINE pre-processing pipeline."""

import shutil
import sys
from pathlib import Path


REQUIRED_COMMANDS = ["bedtools"]
REQUIRED_PACKAGES = ["pandas", "requests", "pybedtools"]


def check_python_version(min_major=3, min_minor=7):
    if sys.version_info < (min_major, min_minor):
        raise RuntimeError(
            f"Python >= {min_major}.{min_minor} is required; "
            f"found {sys.version_info.major}.{sys.version_info.minor}"
        )


def check_commands(commands=None):
    commands = commands or REQUIRED_COMMANDS
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if missing:
        raise RuntimeError(
            "Required command(s) not found on PATH: " + ", ".join(missing) +
            "\nInstall bedtools (e.g. conda install -c bioconda bedtools) and retry."
        )


def check_packages(packages=None):
    packages = packages or REQUIRED_PACKAGES
    missing = []
    for pkg in packages:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        raise RuntimeError(
            "Required Python package(s) missing: " + ", ".join(missing) +
            "\nInstall with: pip install -r requirements.txt (or the preprocess extras)."
        )


def check_static_inputs(repo_root: Path):
    """Verify biocurated static files are present."""
    required = [
        "CREbedDBenhancers_10092018",
        "pantherGeneList.txt",
        "tissuetable_10092018.txt",
    ]
    missing = [name for name in required if not (repo_root / name).is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing required static biocurated file(s): " + ", ".join(missing) +
            "\nThese are treated as static inputs and must be placed in the repository root."
        )


def run_all_checks(repo_root: Path):
    print("Running environment checks...")
    check_python_version()
    check_commands()
    check_packages()
    print("All environment checks passed.")
