# pipeline.py
import os
import uuid
import shutil
from pathlib import Path

# local imports from your package
from .builders import (
    _build_genpartnew_from_dict,
    _build_distrib3d_from_dict,
    _build_disrealnew_from_dict,
)
from .executors import (
    run_genpartnew,
    run_distrib3d,
    run_disrealnew,
)
from .utils import is_windows

# Path to the original _bin inside the installed package (for initial copy)
_ORIG_BIN_DIR = Path(__file__).parent / "_bin"


def run_pipeline(id: str,
                 genpartnew_input: str | dict,
                 distrib3d_input: str | dict,
                 disrealnew_input: str | dict):

    # results next to the caller script (fallback to CWD)
    try:
        import inspect
        caller_file = Path(inspect.stack()[1].filename).resolve()
        script_dir = caller_file.parent
    except Exception:
        script_dir = Path.cwd()

    results_dir = script_dir / f"result_{id}"
    results_dir.mkdir(parents=True, exist_ok=True)

    # --- ensure local _bin exists next to the caller script ---
    local_bin = script_dir / "_bin"
    if not local_bin.exists():
        shutil.copytree(_ORIG_BIN_DIR, local_bin)
    bin_dir = local_bin.resolve()
    # ----------------------------------------------------------

    run_id  = "gp_" + uuid.uuid4().hex[:6]
    run_dir = bin_dir / run_id
    run_dir.mkdir(parents=False, exist_ok=False)

    # pick exe paths
    if is_windows():
        exe  = bin_dir / "genpartnew.exe"
        exe2 = bin_dir / "distrib3d.exe"
        exe3 = bin_dir / "disrealnew.exe"
    else:
        exe  = bin_dir / "genpartnew"
        exe2 = bin_dir / "distrib3d"
        exe3 = bin_dir / "disrealnew"

    # ---- substitute tokens {ID} and {RUN} in the provided inputs ----
    # genpartnew
    
    part_name = Path(genpartnew_input["out_particle_ids"].format(ID=id)).name
    genpartnew_text = _build_genpartnew_from_dict(id, genpartnew_input)

   
    phase_name = Path(distrib3d_input["out_name"].format(ID=id)).name
    distrib3d_text = _build_distrib3d_from_dict(id, run_id, distrib3d_input)

    # disrealnew
    disrealnew_text = _build_disrealnew_from_dict(id, run_id, disrealnew_input)
    # ----------------------------------------------------------------

    # === execute steps via executors ===

    # genpartnew (cwd=run_dir)
    run_genpartnew(
        id=id,
        genpartnew_input=genpartnew_text,
        exe=exe,
        run_dir=run_dir,
        results_dir=results_dir,
        part_name=part_name,
    )

    # distrib3d (cwd=bin_dir)
    run_distrib3d(
        id=id,
        distrib3d_input=distrib3d_text,
        exe2=exe2,
        bin_dir=bin_dir,
        run_dir=run_dir,
        results_dir=results_dir,
        phase_name=phase_name,
    )

    # disrealnew (cwd=bin_dir) + finalize (mirrors and cleanup)
    # NOTE: this function mirrors run_dir -> results_dir and removes run_dir,
    # just like your original long pipeline did at the end.
    run_disrealnew(
        id=id,
        disrealnew_input=disrealnew_text,
        exe3=exe3,
        bin_dir=bin_dir,
        run_dir=run_dir,
        results_dir=results_dir,
        phase_name=phase_name,
        part_name=part_name,
    )

    return results_dir
