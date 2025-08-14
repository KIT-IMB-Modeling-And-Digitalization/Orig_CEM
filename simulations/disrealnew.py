import package_test.module as m
from tempfile import NamedTemporaryFile
from pathlib import Path
import os, uuid, subprocess, shutil

# ---- New: run ID and results folder next to this script ----
id = "01"  # choose your run label
script_dir = Path(__file__).parent.resolve()
results_dir = script_dir / f"result_{id}"
results_dir.mkdir(parents=True, exist_ok=True)

# Installed _bin (where the executable + assets live)
bin_dir = Path(getattr(m, "_BIN_DIR", Path(m.__file__).parent / "_bin")).resolve()

# Per-run folder inside _bin (short + unique)
run_id  = "gp_" + uuid.uuid4().hex[:6]
run_dir = bin_dir / run_id
run_dir.mkdir(parents=False, exist_ok=False)

# Use ONLY short filenames in stdin (we'll run with cwd=run_dir)
genpartnew_input = "\n".join([
    "-3034",
    "2",
    "16",
    "0",
    "0.0604",
    "0.515 0.041",
    "1","17","1",
    "1","15","1",
    "1","14","1",
    "1","13","1",
    "2","12","1",
    "2","11","1",
    "4","10","1",
    "5","9","1",
    "8","8","1",
    "13","7","1",
    "21","6","1",
    "38","5","1",
    "73","4","1",
    "174","3","1",
    "450","2","1",
    "2674","1","1",
    "4",
    "3",
    "1",
    "8",
    f"cem140w04floc_{id}.img",      # short names (no long absolute paths)
    f"pcem140w04floc_{id}.img",
    "1"
])

# Write input to a temp file
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(genpartnew_input.rstrip() + "\n")
    temp_in = f.name

# Run inside the per-run folder; capture all terminal output to genpartnew.out
exe = bin_dir / "genpartnew"
log_path = run_dir / f"genpartnew_{id}.out"
with open(temp_in, "rb") as fin, open(log_path, "wb") as flog:
    subprocess.run([str(exe)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                   check=True, cwd=run_dir)

# Temp file cleanup (optional)
try:
    os.remove(temp_in)
except OSError:
    pass

# Copy outputs + log into result_{id} (keep originals in _bin)
for name in [f"cem140w04floc_{id}.img", f"pcem140w04floc_{id}.img", f"genpartnew_{id}.out"]:
    src = run_dir / name
    if src.exists():
        shutil.copy2(src, results_dir / name)

print(f"[genpartnew] Per-run outputs in: {run_dir}")
print(f"[genpartnew] Copied to results folder: {results_dir}")

# =========================
# Append: run distrib3d now
# =========================

# Build distrib3d stdin to use the SAME per-run folder inside _bin
# Input is the genpartnew output in run_dir; output also goes to run_dir
in_name  = f"{run_id}/cem140w04floc_{id}.img"
out_name = f"{run_id}/cement140w04flocf_{id}.img"

distrib3d_input = "\n".join([
    "-99",
    in_name,            # input under the same sandbox
    "cement140",        # correlation filters root in _bin
    out_name,           # output under the same sandbox
    "0.7344",
    "0.6869",
    "0.0938",
    "0.1337",
    "0.1311",
    "0.1386",
    "0.0407",
    "0.0408"
])

# Write temp stdin and run distrib3d with cwd=_bin; log to distrib3d.out in run_dir
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(distrib3d_input.rstrip() + "\n")
    temp_in2 = f.name

exe2 = bin_dir / "distrib3d"
log_path2 = run_dir / f"distrib3d_{id}.out"
with open(temp_in2, "rb") as fin, open(log_path2, "wb") as flog:
    subprocess.run([str(exe2)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                   check=True, cwd=bin_dir)

# Temp file cleanup (optional)
try:
    os.remove(temp_in2)
except OSError:
    pass

# Optionally copy distrib3d outputs + log to results folder (keeping originals in _bin)
for name in [f"cement140w04flocf_{id}.img", f"distrib3d_{id}.out"]:
    src = run_dir / name
    if src.exists():
        shutil.copy2(src, results_dir / name)

print(f"[distrib3d] Used sandbox: {run_dir}")
print(f"[distrib3d] Copied outputs to: {results_dir}")

# =========================
# Append: run disrealnew now (fixed)
# =========================

# Inputs for disrealnew come from the SAME per-run folder:
#  - phase microstructure: cement140w04flocf.img (from distrib3d)
#  - particle ID microstructure: pcem140w04floc.img (from genpartnew)
phase_name = f"cement140w04flocf_{id}.img"
part_name  = f"pcem140w04floc_{id}.img"

phase_path = run_dir / phase_name
part_path  = run_dir / part_name
if not phase_path.exists():
    raise FileNotFoundError(f"[disrealnew] Missing phase microstructure in sandbox: {phase_path}")
if not part_path.exists():
    raise FileNotFoundError(f"[disrealnew] Missing particle-ID microstructure in sandbox: {part_path}")

# IMPORTANT: run disrealnew with cwd = _bin so it can find its runtime data there.
# Use paths prefixed with the per-run folder so reads/writes land in run_dir.
in_phase = f"{run_id}/{phase_name}"
in_part  = f"{run_id}/{part_name}"

disrealnew_input = "\n".join([
    "-2794",                    # seed
    in_phase,                   # 3-D phase ID microstructure under run_id/
    "1 2 3 4 5 6 7 28 26",      # phase assignments
    "35",                       # C3A in fly ash ID
    in_part,                    # 3-D particle ID microstructure under run_id/
    "44990","1",
    "5850","2",
    "8692","3",
    "2631","4",
    "1100","5",
    "2062","6",
    "839","7",
    "0",
    "1000","0","500",
    "0.0001 9000.",
    "0.01 9000.",
    "0.00002 10000.",
    "0.002 2500.",
    "50","5","5000","100",
    "0.00","20.0","20.0","0.0",
    "40.0","83.14","80.0",
    "0.00035","0.72",
    "0","0","1",
    "10","1.0","0","0","1"
]) + "\n"  # ensure trailing newline

# Snapshot files in _bin before running (to catch anything disrealnew drops into _bin root)
pre_files = {p.name for p in bin_dir.iterdir() if p.is_file()}

# Write temp stdin and run with cwd=_bin; log to disrealnew.out inside run_dir
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(disrealnew_input)
    temp_in3 = f.name

exe3 = bin_dir / "disrealnew"
log_path3 = run_dir / f"disrealnew_{id}.out"
with open(temp_in3, "rb") as fin, open(log_path3, "wb") as flog:
    subprocess.run([str(exe3)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                   check=True, cwd=bin_dir)

try:
    os.remove(temp_in3)
except OSError:
    pass

# Collect any newly created files that landed in _bin root and move them into the run_dir
post_files = {p.name for p in bin_dir.iterdir() if p.is_file()}
for name in sorted(post_files - pre_files):
    # Skip our inputs (they live in run_dir already)
    if name in {phase_name, part_name}:
        continue
    shutil.move(str(bin_dir / name), str(run_dir / name))

print(f"[disrealnew] Completed in sandbox: {run_dir}")

# ===========================================================
# Finalize: copy EVERYTHING from sandbox → results, then remove sandbox
# ===========================================================
for item in run_dir.iterdir():
    dst = results_dir / item.name
    if item.is_dir():
        shutil.copytree(item, dst, dirs_exist_ok=True)
    else:
        shutil.copy2(item, dst)

shutil.rmtree(run_dir, ignore_errors=True)
print(f"[finalize] Copied all sandbox files to: {results_dir}")
print(f"[finalize] Removed sandbox: {run_dir}")
