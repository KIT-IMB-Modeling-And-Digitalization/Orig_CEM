import os
import platform
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile  # new import for this function
import uuid, shutil  # new imports used by the function


# Path to the original _bin inside the installed package (for initial copy)
_ORIG_BIN_DIR = Path(__file__).parent / "_bin"


def _is_windows():
    return platform.system() == "Windows"


def run_pipeline(id: str,
                 genpartnew_input: str,
                 distrib3d_input: str,
                 disrealnew_input: str):

    # results next to the caller script (fallback to CWD)
    try:
        import inspect
        caller_file = Path(inspect.stack()[1].filename).resolve()
        script_dir = caller_file.parent
    except Exception:
        script_dir = Path.cwd()

    results_dir = script_dir / f"result_{id}"
    results_dir.mkdir(parents=True, exist_ok=True)

    # --- NEW: ensure local _bin exists next to the caller script ---
    local_bin = script_dir / "_bin"
    if not local_bin.exists():
        shutil.copytree(_ORIG_BIN_DIR, local_bin)
    bin_dir = local_bin.resolve()
    # --------------------------------------------------------------

    run_id  = "gp_" + uuid.uuid4().hex[:6]
    run_dir = bin_dir / run_id
    run_dir.mkdir(parents=False, exist_ok=False)

    # pick exe paths
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
        part_name = Path(genpartnew_input["out_particle_ids"].format(ID=id)).name
        genpartnew_input = _build_genpartnew_from_dict(id, genpartnew_input)
    else:
        genpartnew_input = genpartnew_input.format(ID=id, RUN=run_id)
        _g_lines = [ln.strip() for ln in genpartnew_input.strip().splitlines()]
        # tail is: ..., "8", <out_image>, <out_particle_ids>, "1"
        part_name = Path(_g_lines[-2]).name
        
    if isinstance(distrib3d_input, dict):
    # remember final image basename for copies/checks
        phase_name = Path(distrib3d_input["out_name"].format(ID=id)).name
        distrib3d_input = _build_distrib3d_from_dict(id, run_id, distrib3d_input)
    else:
        distrib3d_input = distrib3d_input.format(ID=id, RUN=run_id)
        _d_lines = [ln.strip() for ln in distrib3d_input.strip().splitlines()]
        phase_name = Path(_d_lines[3]).name  # 4th line = out_name

    if isinstance(disrealnew_input, dict):
    # (optional) derive part_name from dict, used later in your skip set
        disrealnew_input = _build_disrealnew_from_dict(id, run_id, disrealnew_input)
    else:
        # keep string path behavior
        disrealnew_input = disrealnew_input.format(ID=id, RUN=run_id)
    # ---------------------------------------------------------------------



    # phase_name = f"cement140w04flocf_{id}.img"
    # part_name  = f"pcem140w04floc_{id}.img"
    

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
    for name in [phase_name, f"distrib3d_{id}.out"]:
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

def _build_distrib3d_from_dict(id: str, run_id: str, cfg: dict) -> str:
    """
    Build stdin for distrib3d from a config dict.
    Fields:
      seed: int
      in_name: str          (e.g., "cem140w04floc_{ID}.img")
      filters_root: str     (e.g., "cement140")
      out_name: str         (e.g., "cement140w04flocf_{ID}.img")
      targets: list[str|num] 8 values in order:
               C3S vf, C3S sa, C2S vf, C2S sa, C3A vf, C3A sa, C4AF vf, C4AF sa
    Adds "{RUN}/" prefix to in/out if they’re not absolute or already a path.
    """
    from pathlib import Path as _P

    def _fmt_name(name: str) -> str:
        s = name.format(ID=id)
        # Prefix with run_id unless absolute or already contains a path separator
        if _P(s).is_absolute() or ("/" in s or "\\" in s):
            return s
        return f"{run_id}/{s}"

    seed = str(cfg.get("seed", -99))
    in_name = _fmt_name(cfg["in_name"])
    filters_root = str(cfg.get("filters_root", "cement140"))
    out_name = _fmt_name(cfg["out_name"])

    targets = cfg["targets"]
    if len(targets) != 8:
        raise ValueError(
            "distrib3d.targets must have exactly 8 values in order: "
            "C3S vf, C3S sa, C2S vf, C2S sa, C3A vf, C3A sa, C4AF vf, C4AF sa."
        )
    targets = [str(x) for x in targets]

    lines = [seed, in_name, filters_root, out_name] + targets
    return "\n".join(lines)

def _build_disrealnew_from_dict(id: str, run_id: str, cfg: dict) -> str:
    """
    Build stdin for disrealnew from a config dict (exact order preserved).
    """
    from pathlib import Path as _P

    def _fmt_name(name: str) -> str:
        s = str(name).format(ID=id)
        # Prefix with run_id unless absolute or already contains a path separator
        if _P(s).is_absolute() or ("/" in s or "\\" in s):
            return s
        return f"{run_id}/{s}"

    # phase map: list -> "1 2 3 ..." ; or accept string as-is
    pm = cfg.get("phase_map")
    if isinstance(pm, (list, tuple)):
        phase_map = " ".join(map(str, pm))
    else:
        phase_map = str(pm)

    lines = [
        str(cfg.get("seed", -2794)),
        _fmt_name(cfg["phase_file"]),
        phase_map,
        str(cfg.get("c3a_fa", 35)),
        _fmt_name(cfg["part_file"]),
    ]

    # 7 pairs (count, phase_id)
    one_px_pairs = cfg.get("one_px_pairs", [])
    for count, phase_id in one_px_pairs:
        lines += [str(count), str(phase_id)]

    # trailing lone count
    lines.append(str(cfg.get("one_px_extra", "0")))

    # main controls
    lines += [
        str(cfg["cycles"]),
        str(cfg["sat_flag"]),
        str(cfg["max_diff"]),
    ]

    # 4 nucleation params
    nuc = cfg["nuc_params"]
    if len(nuc) != 4:
        raise ValueError("disrealnew.nuc_params must have 4 entries.")
    lines += list(map(str, nuc))

    # 4 freqs
    fr = cfg["freqs"]
    if len(fr) != 4:
        raise ValueError("disrealnew.freqs must have 4 entries.")
    lines += list(map(str, fr))

    # 4 thermal
    th = cfg["thermal"]
    if len(th) != 4:
        raise ValueError("disrealnew.thermal must have 4 entries.")
    lines += list(map(str, th))

    # 3 activation energies
    ea = cfg["Ea"]
    if len(ea) != 3:
        raise ValueError("disrealnew.Ea must have 3 entries.")
    lines += list(map(str, ea))

    # final scalars
    lines += [
        str(cfg["cycle_to_time"]),
        str(cfg["agg_vf"]),
    ]

    # 8 flags
    fl = cfg["flags"]
    if len(fl) != 8:
        raise ValueError("disrealnew.flags must have 8 entries.")
    lines += list(map(str, fl))

    text = "\n".join(lines)
    if not text.endswith("\n"):
        text += "\n"
    return text
