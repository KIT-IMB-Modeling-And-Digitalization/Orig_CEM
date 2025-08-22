# executors.py
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
    # === run genpartnew (cwd=run_dir) ===
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
    # === run distrib3d (cwd=_BIN_DIR) ===
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
    # === run disrealnew (cwd=_BIN_DIR) ===
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
