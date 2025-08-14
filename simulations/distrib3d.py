import package_test.module as m
from pathlib import Path
from tempfile import NamedTemporaryFile
import os, shutil, sys, uuid

# Working directory (where you launched this script)
here = Path.cwd()

# Installed _bin (where the executable + assets live)
bin_dir = Path(getattr(m, "_BIN_DIR", Path(m.__file__).parent / "_bin")).resolve()

# Create a short, unique sandbox inside _bin (visible during the run)
run_id  = "r_" + uuid.uuid4().hex[:6]
run_dir = bin_dir / run_id
run_dir.mkdir(parents=False, exist_ok=False)

# Inputs in your working dir
inp = (here / "cem140w04floc.img").resolve()
if not inp.exists():
    sys.exit(f"[distrib3d] Missing input image: {inp}")
out_basename = "cement140w04flocf.img"  # change if you want a different output file name

# Link (or copy) the input into the sandbox under a short local name
in_sbx = run_dir / inp.name
if in_sbx.exists() or in_sbx.is_symlink():
    in_sbx.unlink()
try:
    os.symlink(inp, in_sbx)
except OSError:
    shutil.copy(inp, in_sbx)

# Build stdin using ONLY short paths inside _bin (relative to cwd=_bin)
in_name  = f"{run_id}/{inp.name}"
out_name = f"{run_id}/{out_basename}"

distrib3d_input = "\n".join([
    "-99",
    in_name,            # input under sandbox
    "cement140",        # correlation filters root in _bin
    out_name,           # output under sandbox
    "0.7344",
    "0.6869",
    "0.0938",
    "0.1337",
    "0.1311",
    "0.1386",
    "0.0407",
    "0.0408"
])

# Write temp stdin and run (module runs with cwd=_bin)
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(distrib3d_input.rstrip() + "\n")
    temp_in = f.name

try:
    m.run_distrib3d(temp_in)
finally:
    try:
        os.remove(temp_in)
    except OSError:
        pass

# Move everything produced in the sandbox back to the working dir (except the input link/copy)
for item in sorted(run_dir.iterdir()):
    if item.name == in_sbx.name:
        continue
    dst = here / item.name
    if dst.exists():
        dst.unlink()
    shutil.move(str(item), str(dst))

# Cleanup sandbox
try:
    if in_sbx.exists() or in_sbx.is_symlink():
        in_sbx.unlink()
    run_dir.rmdir()
except OSError:
    shutil.rmtree(run_dir, ignore_errors=True)
