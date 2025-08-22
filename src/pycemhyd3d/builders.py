# builders.py
# This file is part of pycemhyd3d that provides utility functions for building input configurations.

from pathlib import Path



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
