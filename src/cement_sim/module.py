import os
import platform
import subprocess
from pathlib import Path


# Path to the _bin folder inside the installed package
_BIN_DIR = Path(__file__).parent / "_bin"

def _is_windows():
    return platform.system() == "Windows"
    
# def run_cemhyd3d(sim_id, inputs:Dict[str:list]):
#     """run all three functions""""
    
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


# --- append to module.py ---

from tempfile import NamedTemporaryFile  # new import for this function
import uuid, shutil  # new imports used by the function

# BEFORE:
# def run_pipeline(id: str):

# AFTER:
def run_pipeline(id: str,
                 genpartnew_input: str,
                 distrib3d_input: str,
                 disrealnew_input: str):

    # NEW: put results next to the caller script (fallback to CWD)
    try:
        import inspect
        caller_file = Path(inspect.stack()[1].filename).resolve()
        script_dir = caller_file.parent
    except Exception:
        script_dir = Path.cwd()

    results_dir = script_dir / f"result_{id}"
    results_dir.mkdir(parents=True, exist_ok=True)

    bin_dir = Path(_BIN_DIR).resolve()

    run_id  = "gp_" + uuid.uuid4().hex[:6]
    run_dir = bin_dir / run_id
    run_dir.mkdir(parents=False, exist_ok=False)

    # pick exe paths (unchanged) ...
    if _is_windows():
        exe  = bin_dir / "genpartnew.exe"
        exe2 = bin_dir / "distrib3d.exe"
        exe3 = bin_dir / "disrealnew.exe"
    else:
        exe  = bin_dir / "genpartnew"
        exe2 = bin_dir / "distrib3d"
        exe3 = bin_dir / "disrealnew"

    # ---- NEW: substitute tokens {ID} and {RUN} in the provided inputs ----
    if isinstance(genpartnew_input, dict):
        genpartnew_input = _build_genpartnew_from_dict(id, genpartnew_input)
    else:
        genpartnew_input = genpartnew_input.format(ID=id, RUN=run_id)
    distrib3d_input  = distrib3d_input.format(ID=id, RUN=run_id)
    disrealnew_input = disrealnew_input.format(ID=id, RUN=run_id)
    # ---------------------------------------------------------------------



    phase_name = f"cement140w04flocf_{id}.img"
    part_name  = f"pcem140w04floc_{id}.img"
    

    # === run genpartnew (cwd=run_dir) ===
    with NamedTemporaryFile('w+', delete=False) as f:
        f.write(genpartnew_input.rstrip() + "\n")
        temp_in = f.name
    log_path = run_dir / f"genpartnew_{id}.out"
    with open(temp_in, "rb") as fin, open(log_path, "wb") as flog:
        subprocess.run([str(exe)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                       check=True, cwd=run_dir)
    try: os.remove(temp_in)
    except OSError: pass
    for name in [f"cem140w04floc_{id}.img", f"pcem140w04floc_{id}.img", f"genpartnew_{id}.out"]:
        src = run_dir / name
        if src.exists():
            shutil.copy2(src, results_dir / name)
    print(f"[genpartnew] Per-run outputs in: {run_dir}")
    print(f"[genpartnew] Copied to results folder: {results_dir}")

    # === run distrib3d (cwd=_BIN_DIR) ===
    with NamedTemporaryFile('w+', delete=False) as f:
        f.write(distrib3d_input.rstrip() + "\n")
        temp_in2 = f.name
    log_path2 = run_dir / f"distrib3d_{id}.out"
    with open(temp_in2, "rb") as fin, open(log_path2, "wb") as flog:
        subprocess.run([str(exe2)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                       check=True, cwd=bin_dir)
    try: os.remove(temp_in2)
    except OSError: pass
    for name in [f"cement140w04flocf_{id}.img", f"distrib3d_{id}.out"]:
        src = run_dir / name
        if src.exists():
            shutil.copy2(src, results_dir / name)
    print(f"[distrib3d] Used sandbox: {run_dir}")
    print(f"[distrib3d] Copied outputs to: {results_dir}")

    # === run disrealnew (cwd=_BIN_DIR) ===
    phase_path = run_dir / phase_name
    part_path  = run_dir / part_name
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
    try: os.remove(temp_in3)
    except OSError: pass

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
    # ====== end moved code ======


def build_genpartnew_input(id: str, cfg: dict) -> str:
    # Split can be provided as list ["0.515","0.041"] or a single "0.515 0.041"
    if isinstance(cfg["calcium_sulfate_split"], (list, tuple)):
        split = " ".join(map(str, cfg["calcium_sulfate_split"]))
    else:
        split = str(cfg["calcium_sulfate_split"])

    lines = [
        str(cfg["seed"]),
        str(cfg.get("place_menu", 2)),
        str(cfg["n_size_classes"]),
        str(cfg.get("dispersion_px", 0)),
        str(cfg["calcium_sulfate_vf"]),
        split,
    ]

    for count, radius, phase in cfg["size_classes"]:
        lines += [str(count), str(radius), str(phase)]

    lines += [
        str(cfg.get("report_phase_counts_menu", 4)),
        str(cfg.get("flocculate_menu", 3)),
        str(cfg.get("n_flocs", 1)),
        str(cfg.get("output_menu", 8)),
        cfg["out_image"].format(ID=id),
        cfg["out_particle_ids"].format(ID=id),
        str(cfg.get("exit_menu", 1)),
    ]

    return "\n".join(lines)

def _build_genpartnew_from_dict(id: str, cfg: dict) -> str:
    # cfg schema:
    # {
    #   "seed": -3034,
    #   "place_menu": 2,
    #   "n_size_classes": 16,
    #   "dispersion_px": 0,
    #   "calcium_sulfate_vf": "0.0604",
    #   "calcium_sulfate_split": ["0.515", "0.041"] or "0.515 0.041",
    #   "size_classes": [(count, radius, phase_id), ...]  # 16 tuples
    #   "report_phase_counts_menu": 4, "flocculate_menu": 3,
    #   "n_flocs": 1, "output_menu": 8,
    #   "out_image": "cem140w04floc_{ID}.img",
    #   "out_particle_ids": "pcem140w04floc_{ID}.img",
    #   "exit_menu": 1
    # }
    split = cfg["calcium_sulfate_split"]
    if isinstance(split, (list, tuple)):
        split = " ".join(map(str, split))
    lines = [
        str(cfg["seed"]),
        str(cfg.get("place_menu", 2)),
        str(cfg["n_size_classes"]),
        str(cfg.get("dispersion_px", 0)),
        str(cfg["calcium_sulfate_vf"]),
        str(split),
    ]
    for c, r, p in cfg["size_classes"]:
        lines += [str(c), str(r), str(p)]
    lines += [
        str(cfg.get("report_phase_counts_menu", 4)),
        str(cfg.get("flocculate_menu", 3)),
        str(cfg.get("n_flocs", 1)),
        str(cfg.get("output_menu", 8)),
        cfg.get("out_image", "cem140w04floc_{ID}.img").format(ID=id),
        cfg.get("out_particle_ids", "pcem140w04floc_{ID}.img").format(ID=id),
        str(cfg.get("exit_menu", 1)),
    ]
    return "\n".join(lines)
