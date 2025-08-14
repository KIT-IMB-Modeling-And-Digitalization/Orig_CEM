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
    "cem140w04floc.img",      # short names (no long absolute paths)
    "pcem140w04floc.img",
    "1"
])

# Write input to a temp file
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(genpartnew_input.rstrip() + "\n")
    temp_in = f.name

# Run inside the per-run folder; capture all terminal output to genpartnew.out
exe = bin_dir / "genpartnew"
log_path = run_dir / "genpartnew.out"
with open(temp_in, "rb") as fin, open(log_path, "wb") as flog:
    subprocess.run([str(exe)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                   check=True, cwd=run_dir)

# Temp file cleanup (optional)
try:
    os.remove(temp_in)
except OSError:
    pass

# Copy outputs + log into result_{id} (keep originals in _bin)
for name in ["cem140w04floc.img", "pcem140w04floc.img", "genpartnew.out"]:
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
in_name  = f"{run_id}/cem140w04floc.img"
out_name = f"{run_id}/cement140w04flocf.img"

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
log_path2 = run_dir / "distrib3d.out"
with open(temp_in2, "rb") as fin, open(log_path2, "wb") as flog:
    subprocess.run([str(exe2)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT,
                   check=True, cwd=bin_dir)

# Temp file cleanup (optional)
try:
    os.remove(temp_in2)
except OSError:
    pass

# Optionally copy distrib3d outputs + log to results folder (keeping originals in _bin)
for name in ["cement140w04flocf.img", "distrib3d.out"]:
    src = run_dir / name
    if src.exists():
        shutil.copy2(src, results_dir / name)

print(f"[distrib3d] Used sandbox: {run_dir}")
print(f"[distrib3d] Copied outputs to: {results_dir}")
