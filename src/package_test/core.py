import subprocess
import sys
import importlib.resources as res
from pathlib import Path
import platform

def run_genpartnew(args=None):
    """
    Run the packaged genpartnew executable with optional arguments.
    If args is None, runs without arguments (and the program may prompt for input).
    """
    exe_name = "genpartnew.exe" if platform.system() == "Windows" else "genpartnew"
    exe_path = res.files("package_test") / "_bin" / exe_name

    if not exe_path.exists():
        raise FileNotFoundError(f"Executable not found: {exe_path}")

    # Build the command: binary plus args list
    cmd = [str(exe_path)]
    if args:
        if isinstance(args, str):
            cmd.append(args)
        else:
            cmd.extend(args)

    # Run the executable, streaming output to the console
    result = subprocess.run(cmd, check=False)
    return result.returncode
