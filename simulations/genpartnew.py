from package_test.core import run_genpartnew  # keeps dependency, even though we launch via python -c
import subprocess
import os
from pathlib import Path

# === Input for genpartnew ===
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
    "cem140w04floc.img",
    "pcem140w04floc.img",
    "1"
])

if __name__ == "__main__":
    work_dir = os.getcwd()
    out_path = Path(work_dir) / "genpartnew.out"

    # run your existing function in a child python, pipe stdin, and capture output to file
    with open(out_path, "wb") as fo:
        proc = subprocess.Popen(
            ["python", "-c", "from package_test.core import run_genpartnew; run_genpartnew()"],
            stdin=subprocess.PIPE,
            stdout=fo,
            stderr=subprocess.STDOUT,
            cwd=work_dir
        )
        # ensure a trailing newline so the last menu input is flushed properly
        proc.communicate(input=(genpartnew_input + "\n").encode("utf-8"))

    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"✅ Wrote log to {out_path.resolve()}")
    else:
        print(f"⚠ genpartnew.out not found or empty in {work_dir}")
