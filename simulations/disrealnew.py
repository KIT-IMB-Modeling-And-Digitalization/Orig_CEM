import package_test.module as m
from pathlib import Path
from tempfile import NamedTemporaryFile
import os, shutil, sys, subprocess

here = Path.cwd()
bin_dir = Path(getattr(m, "_BIN_DIR", Path(m.__file__).parent / "_bin")).resolve()

# Inputs expected next to this script (your working dir)
inp_phase = (here / "cem140w04floc.img").resolve()
inp_part  = (here / "pcem140w04floc.img").resolve()
if not inp_phase.exists():
    sys.exit(f"[disrealnew] Missing input image (phase IDs): {inp_phase}")
if not inp_part.exists():
    sys.exit(f"[disrealnew] Missing input image (particle IDs): {inp_part}")

# Short names inside _bin
name_phase = inp_phase.name
name_part  = inp_part.name

# Snapshot files in _bin before we run
pre_files = {p.name for p in bin_dir.iterdir() if p.is_file()}

# Make inputs visible in _bin via short names
link_phase = bin_dir / name_phase
link_part  = bin_dir / name_part
for src, dst in [(inp_phase, link_phase), (inp_part, link_part)]:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy(src, dst)

# Build stdin using only basenames (short names)
disrealnew_input = "\n".join([
    "-2794",                    # seed
    name_phase,                 # 3-D phase ID microstructure
    "1 2 3 4 5 6 7 28 26",      # phase assignments
    "35",                       # C3A in fly ash ID
    name_part,                  # 3-D particle ID microstructure
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

# Write temp stdin
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(disrealnew_input)
    temp_in = f.name

# Run the executable inside _bin so it finds its assets
exe = bin_dir / "disrealnew"
try:
    subprocess.run([str(exe)], stdin=open(temp_in, "rb"), check=True, cwd=bin_dir)
finally:
    try: os.remove(temp_in)
    except OSError: pass

# Collect any newly created files in _bin and move them back to your working dir
post_files = {p.name for p in bin_dir.iterdir() if p.is_file()}
new_files = {n for n in (post_files - pre_files) if n not in {name_phase, name_part}}

for fname in sorted(new_files):
    src = bin_dir / fname
    dst = here / fname
    try:
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), str(dst))
    except Exception as e:
        print(f"[disrealnew] Warning: could not move {src} -> {dst}: {e}")

# Cleanup the linked/copied inputs from _bin to keep it tidy
for p in [link_phase, link_part]:
    try:
        p.unlink()
    except OSError:
        pass
