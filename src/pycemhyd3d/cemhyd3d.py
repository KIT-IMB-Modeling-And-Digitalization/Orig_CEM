'''High-level pipeline orchestration for pycemhyd3d.

Author:
    Omid Jahromi <omid.esmaeelipoor@gmail.com> 

Overview:
    Provides a single entry point `run_cemhyd3d(...)` that:
      1) Creates a per-run sandbox directory under a local `cempy3d/` copy,
      2) Builds stdin payloads for `genpartnew`, `distrib3d`, `disrealnew`,
      3) Executes the three native tools in sequence,
      4) Mirrors all produced files into `./results/result_<ID>/`,
      5) Cleans up the per-run sandbox and temporary local `cempy3d/`.

Notes:
    - Executables are resolved from a package-shipped `cempy3d/` directory and
      copied next to the caller if needed for the run.
    - Windows vs. Linux executable names are handled (e.g., `.exe` suffix).

'''
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
_ORIG_BIN_DIR = Path(__file__).parent / "cempy3d"


def run_cemhyd3d(id: str,
                 genpartnew_input: str | dict,
                 distrib3d_input: str | dict,
                 disrealnew_input: str | dict):

    '''Run the full CEMHYD3D pipeline: genpartnew → distrib3d → disrealnew.

    Args:
        id: Short identifier used in filenames and result folder naming.
        genpartnew_input: Either a ready-made stdin string or a dict config
            consumable by `_build_genpartnew_from_dict`.
        distrib3d_input: Either a ready-made stdin string or a dict config
            consumable by `_build_distrib3d_from_dict`.
        disrealnew_input: Either a ready-made stdin string or a dict config
            consumable by `_build_disrealnew_from_dict`.

    Side effects:
        - Creates `results/result_<id>/` next to the caller script.
        - Copies a temporary `cempy3d/` folder next to the caller (removed at end).
        - Creates a per-run sandbox inside that folder and deletes it after mirroring.

    Returns:
        pathlib.Path: The path to the results directory.

    Raises:
        FileNotFoundError: if required intermediate files are missing before
            running `disrealnew` (checked inside executors).
        subprocess.CalledProcessError: if any executable returns non-zero.
    '''
    try:
        import inspect
        caller_file = Path(inspect.stack()[1].filename).resolve()
        script_dir = caller_file.parent
    except Exception:
        script_dir = Path.cwd()

    results_dir = script_dir / "results" / f"result_{id}"
    results_dir.mkdir(parents=True, exist_ok=True)

    # --- ensure local cempy3d exists next to the caller script (idempotent) ---
    local_bin = script_dir / "cempy3d"
    marker = local_bin / ".ready"

# If missing OR previous copy was interrupted (no marker), rebuild cleanly
    if (not local_bin.exists()) or (not marker.exists()):
        # remove any partial/old dir
        shutil.rmtree(local_bin, ignore_errors=True)

        # copy into a temp dir, then move into place to avoid half-copies
        tmp_dir = script_dir / f".cempy3d.tmp-{os.getpid()}-{uuid.uuid4().hex[:6]}"
        shutil.copytree(_ORIG_BIN_DIR, tmp_dir)

        # move temp -> final
        # (shutil.move is fine here since we just removed local_bin)
        shutil.move(str(tmp_dir), str(local_bin))

        # write a marker to signal a complete copy
        marker.write_text("ok", encoding="utf-8")

    bin_dir = local_bin.resolve()
# --------------------------------------------------------------------------

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
    shutil.rmtree(local_bin, ignore_errors=True)
    print(f"✅ Pipeline completed successfully. Results in: {results_dir}")
    return results_dir
