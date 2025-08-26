'''Process execution helpers for pycemhyd3d.

Author:
    Omid Jahromi <omid.esmaeelipoor@gmail.com>

Overview:
    These functions encapsulate low-level details of invoking the native
    executables (`genpartnew`, `distrib3d`, `disrealnew`), writing temporary
    stdin files, choosing correct working directories, and mirroring outputs.

Conventions:
    - Each function accepts already-rendered stdin payloads (strings).
    - Working directories:
        * genpartnew: runs in the per-run sandbox (`cwd=run_dir`)
        * distrib3d: runs in the executable folder (`cwd=bin_dir`)
        * disrealnew: runs in the executable folder (`cwd=bin_dir`)
    - Logs are written as `<tool>_<id>.out` inside the sandbox.
'''
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile


def run_genpartnew(id: str,
                   genpartnew_input: str,
                   exe: Path,
                   run_dir: Path,
                   results_dir: Path,
                   part_name: str):
    '''Run `genpartnew` with the given stdin payload.

    Args:
        id: Identifier used for log/IO filenames.
        genpartnew_input: Text to send to stdin.
        exe: Path to the genpartnew executable.
        run_dir: Per-run sandbox directory (working directory).
        results_dir: Destination where selected outputs are copied.
        part_name: Expected name of the particle-ID image produced.

    Effects:
        - Writes a temp stdin file, captures combined stdout+stderr into
          `genpartnew_<id>.out`, and copies key outputs into `results_dir`.
        - Removes the temp stdin file on success/failure.

    Returns:
        str: `part_name` (for downstream use).
    '''
    with NamedTemporaryFile('w+', delete=False) as f:
        f.write(genpartnew_input.rstrip() + "\n")
        temp_in = f.name
    log_path = run_dir / f"genpartnew_{id}.out"
    with open(temp_in, "rb") as fin, open(log_path, "wb") as flog:
        subprocess.run([str(exe)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                       check=True, cwd=run_dir)
    try:
        os.remove(temp_in)
    except OSError:
        pass
    for name in [f"cem140w04floc_{id}.img",
                 f"pcem140w04floc_{id}.img",
                 f"genpartnew_{id}.out"]:
        src = run_dir / name
        if src.exists():
            shutil.copy2(src, results_dir / name)
    print(f"[genpartnew] Per-run outputs in: {run_dir}")
    print(f"[genpartnew] Copied to results folder: {results_dir}")
    return part_name


def run_distrib3d(id: str,
                  distrib3d_input: str,
                  exe2: Path,
                  bin_dir: Path,
                  run_dir: Path,
                  results_dir: Path,
                  phase_name: str):
    '''Run `distrib3d` with the given stdin payload.

    Args:
        id: Identifier used for log/IO filenames.
        distrib3d_input: Text to send to stdin.
        exe2: Path to the distrib3d executable.
        bin_dir: Folder containing the executables (working directory).
        run_dir: Per-run sandbox directory (holds produced files).
        results_dir: Destination where outputs are copied.
        phase_name: Expected output phase image filename.

    Effects:
        - Writes a temp stdin file, logs to `distrib3d_<id>.out`, and copies
          the phase image + log to `results_dir`.

    Returns:
        str: `phase_name` (for downstream use).
    '''
    with NamedTemporaryFile('w+', delete=False) as f:
        f.write(distrib3d_input.rstrip() + "\n")
        temp_in2 = f.name
    log_path2 = run_dir / f"distrib3d_{id}.out"
    with open(temp_in2, "rb") as fin, open(log_path2, "wb") as flog:
        subprocess.run([str(exe2)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                       check=True, cwd=bin_dir)
    try:
        os.remove(temp_in2)
    except OSError:
        pass
    for name in [phase_name, f"distrib3d_{id}.out"]:
        src = run_dir / name
        if src.exists():
            shutil.copy2(src, results_dir / name)
    print(f"[distrib3d] Used sandbox: {run_dir}")
    print(f"[distrib3d] Copied outputs to: {results_dir}")
    return phase_name


def run_disrealnew(id: str,
                   disrealnew_input: str,
                   exe3: Path,
                   bin_dir: Path,
                   run_dir: Path,
                   results_dir: Path,
                   phase_name: str,
                   part_name: str):
    '''Run `disrealnew` with the given stdin payload and finalize outputs.

    Preconditions:
        - `{run_dir}/{phase_name}` and `{run_dir}/{part_name}` must exist.

    Args:
        id: Identifier used for log/IO filenames.
        disrealnew_input: Text to send to stdin.
        exe3: Path to the disrealnew executable.
        bin_dir: Folder containing the executables (working directory).
        run_dir: Per-run sandbox directory (source for final mirroring).
        results_dir: Final results directory (destination).
        phase_name: Name of the phase microstructure file.
        part_name: Name of the particle-ID microstructure file.

    Effects:
        - Runs disrealnew, moves any **new files** created in `bin_dir` into
          `run_dir`, mirrors the entire sandbox into `results_dir`, then
          removes the sandbox.

    Returns:
        pathlib.Path: The `results_dir` path.

    Raises:
        FileNotFoundError: if required inputs are missing in `run_dir`.
    '''
    phase_path = run_dir / phase_name
    part_path = run_dir / part_name
    if not phase_path.exists():
        raise FileNotFoundError(f"[disrealnew] Missing phase microstructure in sandbox: {phase_path}")
    if not part_path.exists():
        raise FileNotFoundError(f"[disrealnew] Missing particle-ID microstructure in sandbox: {part_path}")

    pre_files = {p.name for p in bin_dir.iterdir() if p.is_file()}
    with NamedTemporaryFile('w+', delete=False) as f:
        f.write(disrealnew_input)
        temp_in3 = f.name
    log_path3 = run_dir / f"disrealnew_{id}.out"
    with open(temp_in3, "rb") as fin, open(log_path3, "wb") as flog:
        subprocess.run([str(exe3)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                       check=True, cwd=bin_dir)
    try:
        os.remove(temp_in3)
    except OSError:
        pass

    post_files = {p.name for p in bin_dir.iterdir() if p.is_file()}
    for name in sorted(post_files - pre_files):
        if name in {phase_name, part_name}:
            continue
        shutil.move(str(bin_dir / name), str(run_dir / name))
    print(f"[disrealnew] Completed in sandbox: {run_dir}")

    # === finalize ===
    for item in run_dir.iterdir():
        dst = results_dir / item.name
        if item.is_dir():
            shutil.copytree(item, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(item, dst)
    shutil.rmtree(run_dir, ignore_errors=True)
    print(f"[finalize] Copied all sandbox files to: {results_dir}")
    print(f"[finalize] Removed sandbox: {run_dir}")

    return results_dir
