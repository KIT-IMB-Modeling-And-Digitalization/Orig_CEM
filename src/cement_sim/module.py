import os
import platform
import subprocess
from pathlib import Path

# Path to the _bin folder inside the installed package
_BIN_DIR = Path(__file__).parent / "_bin"

def _is_windows():
    return platform.system() == "Windows"
    
def run_cemhyd3d(sim_id, inputs:Dict[str:list]):
    """run all three functions""""
    
def run_genpartnew(input_file):
    """Run genpartnew executable with the given input file."""
    if _is_windows():
        exe = _BIN_DIR / "genpartnew.exe"
        subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True, cwd=_BIN_DIR)
        return
    # --- Linux (unchanged) ---
    exe = _BIN_DIR / "genpartnew"
    subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True)

def run_distrib3d(input_file):
    """Run distrib3d executable with the given input file.
    Note: runs with cwd=_BIN_DIR so relative asset lookups work."""
    if _is_windows():
        exe = _BIN_DIR / "distrib3d.exe"
        subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True, cwd=_BIN_DIR)
        return
    # --- Linux (unchanged) ---
    exe = _BIN_DIR / "distrib3d"
    subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True, cwd=_BIN_DIR)

def run_disrealnew(input_file):
    """Run disrealnew executable with the given input file."""
    if _is_windows():
        exe = _BIN_DIR / "disrealnew.exe"
        subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True, cwd=_BIN_DIR)
        return
    # --- Linux (unchanged) ---
    exe = _BIN_DIR / "disrealnew"
    subprocess.run([str(exe)], stdin=open(input_file, "rb"), check=True)
