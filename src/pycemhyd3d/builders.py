'''Builders for stdin payloads used by pycemhyd3d executables.

Author:
    Omid Jahromi <omid.esmaeelipoor@gmail.com>

Overview:
    This module converts Python dictionaries into the exact line-oriented text
    formats expected on stdin by the native C executables:
        - genpartnew
        - distrib3d
        - disrealnew

Notes:
    * Paths inside the generated text may be automatically prefixed with a run
      directory identifier when appropriate (see each builder).
    * All builders preserve the ordering required by the original tools.

'''
from pathlib import Path



def _build_genpartnew_from_dict(id: str, cfg: dict) -> str:
    '''Build the stdin string for **genpartnew** from a configuration dict.

    Args:
        id: A short identifier used to substitute into output filenames
            (e.g., `{ID}` tokens in `out_image`, `out_particle_ids`).
        cfg: Configuration mapping with keys (examples / defaults shown):
            seed (int)
            place_menu (int, default 2)
            n_size_classes (int)
            dispersion_px (int, default 0)
            calcium_sulfate_vf (str or number)
            calcium_sulfate_split (list[str|num] | str)
            size_classes (list[tuple]): sequence of (count, radius, phase_id)
            report_phase_counts_menu (int, default 4)
            flocculate_menu (int, default 3)
            n_flocs (int, default 1)
            output_menu (int, default 8)
            out_image (str, may contain {ID})
            out_particle_ids (str, may contain {ID})
            exit_menu (int, default 1)

    Returns:
        A newline-joined string ready to be piped to genpartnew's stdin.
    '''
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
    '''Build the stdin string for **distrib3d** from a configuration dict.

    Behavior:
        - `in_name` and `out_name` are formatted with `{ID}` and, unless given
          as absolute or already containing a path separator, are prefixed with
          `{run_id}/` so files live inside the per-run sandbox.

    Args:
        id: Identifier substituted into filenames (via `{ID}`).
        run_id: Name of the per-run directory used as a prefix for relative IO.
        cfg: Keys:
            seed (int, default -99)
            in_name (str)
            filters_root (str, default "cement140")
            out_name (str)
            targets (list[8]): C3S vf, C3S sa, C2S vf, C2S sa,
                               C3A vf, C3A sa, C4AF vf, C4AF sa

    Returns:
        A newline-joined string for distrib3d's stdin.

    Raises:
        ValueError: if `targets` does not contain exactly 8 values.
    '''
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
    '''Build the stdin string for **disrealnew** from a configuration dict.

    Behavior:
        - `phase_file` and `part_file` are formatted with `{ID}` and, unless
          absolute or already a path, are prefixed with `{run_id}/`.
        - Preserves the exact ordering required by disrealnew.

    Args:
        id: Identifier substituted into filenames (via `{ID}`).
        run_id: Name of the per-run directory used as a prefix for relative IO.
        cfg: Keys (abbreviated):
            seed (int, default -2794)
            phase_file (str)
            phase_map (list|tuple|str): list/tuple joined as space-sep string
            c3a_fa (int, default 35)
            part_file (str)
            one_px_pairs (list[tuple[int,int]], optional): 7 (count, phase_id)
            one_px_extra (int|str, default "0")
            cycles (int), sat_flag (int), max_diff (num)
            nuc_params (list[4]), freqs (list[4]),
            thermal (list[4]), Ea (list[3])
            cycle_to_time (num), agg_vf (num)
            flags (list[8])

    Returns:
        A newline-joined string for disrealnew's stdin (with trailing newline).

    Raises:
        ValueError: if any of the fixed-length arrays have wrong lengths.
    '''
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
