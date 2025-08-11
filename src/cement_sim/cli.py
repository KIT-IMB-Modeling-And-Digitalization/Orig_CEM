"""
cement_sim.cli
==============
Thin CLI wrappers around the packaged native executables.

This expects the three compiled binaries to live at:
    cement_sim/bin/{genpartnew,distrib3d,disrealnew}

Usage examples
--------------
# read input from a file and write output to a file
genpartnew --input-file input.txt --output-file out.txt

# same for the others, with optional working directory
distrib3d --input-file d3.in --output-file d3.out --cwd /path/to/run
disrealnew --input-file dr.in --output-file dr.out --cwd .

You can also pipe input:
cat input.txt | genpartnew --output-file out.txt
"""

from __future__ import annotations

import argparse
import os
import sys
import stat
from typing import Optional

from importlib import resources as r

# We moved linux_func.py to src/cement_sim/linux_func.py
from .linux_func import (
    run_executable_genpartnew,
    run_executable_distrib3d,
    run_executable_disrealnew,
)


# ---------- helpers ----------

def _resolve_exe(exe_name: str) -> str:
    """
    Return an absolute path to the packaged executable inside cement_sim/bin/.
    Raises FileNotFoundError if not present.
    """
    try:
        bin_dir = r.files("cement_sim").joinpath("bin")
    except Exception as e:
        raise FileNotFoundError(f"Could not locate cement_sim/bin in the package: {e}")

    exe_path = os.fspath(bin_dir / exe_name)
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")

    # ensure it's executable (mainly for safety on some filesystems)
    try:
        st = os.stat(exe_path)
        os.chmod(exe_path, st.st_mode | stat.S_IXUSR)
    except Exception:
        # non-fatal; if chmod fails we still try to run it
        pass

    return exe_path


def _read_input_text(input_file: Optional[str]) -> str:
    """
    Read input text either from --input-file or stdin.
    """
    if input_file:
        with open(input_file, "r", encoding="utf-8") as f:
            return f.read()
    # read from stdin (works with pipes)
    return sys.stdin.read()


def _require_output_path(path: Optional[str]) -> str:
    if not path:
        raise SystemExit("ERROR: --output-file is required")
    # create parent directory if needed
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)
    return path


# ---------- command runners (entry points call these) ----------

def run_genpartnew(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        prog="genpartnew",
        description="Run the packaged genpartnew executable.",
    )
    parser.add_argument("--input-file", help="Read input from this file (default: stdin)")
    parser.add_argument("--output-file", required=True, help="Write stdout of the executable here")
    args = parser.parse_args(argv)

    exe = _resolve_exe("genpartnew")
    input_text = _read_input_text(args.input_file)
    output_file = _require_output_path(args.output_file)

    run_executable_genpartnew(executable=exe, input_data=input_text, output_file=output_file)


def run_distrib3d(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        prog="distrib3d",
        description="Run the packaged distrib3d executable.",
    )
    parser.add_argument("--input-file", help="Read input from this file (default: stdin)")
    parser.add_argument("--output-file", required=True, help="Write stdout of the executable here")
    parser.add_argument("--cwd", help="Working directory to run the executable in")
    args = parser.parse_args(argv)

    exe = _resolve_exe("distrib3d")
    input_text = _read_input_text(args.input_file)
    output_file = _require_output_path(args.output_file)

    run_executable_distrib3d(
        executable=exe,
        input_data=input_text,
        output_file=output_file,
        cwd=args.cwd,
    )


def run_disrealnew(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        prog="disrealnew",
        description="Run the packaged disrealnew executable.",
    )
    parser.add_argument("--input-file", help="Read input from this file (default: stdin)")
    parser.add_argument("--output-file", required=True, help="Write stdout of the executable here")
    parser.add_argument("--cwd", help="Working directory to run the executable in")
    args = parser.parse_args(argv)

    exe = _resolve_exe("disrealnew")
    input_text = _read_input_text(args.input_file)
    output_file = _require_output_path(args.output_file)

    run_executable_disrealnew(
        executable=exe,
        input_data=input_text,
        output_file=output_file,
        cwd=args.cwd,
    )


# Allow `python -m cement_sim.cli <subcommand>` style if you want later.
if __name__ == "__main__":
    # simple router if ever used directly
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    cmd, *rest = sys.argv[1:]
    if cmd == "genpartnew":
        run_genpartnew(rest)
    elif cmd == "distrib3d":
        run_distrib3d(rest)
    elif cmd == "disrealnew":
        run_disrealnew(rest)
    else:
        print(f"Unknown subcommand: {cmd}")
        sys.exit(2)
