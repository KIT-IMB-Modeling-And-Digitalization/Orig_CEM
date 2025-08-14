import package_test.module as m
from tempfile import NamedTemporaryFile
from pathlib import Path

here = Path.cwd()

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
    str(here / "cem140w04floc.img"),   # absolute output path
    str(here / "pcem140w04floc.img"),  # absolute output path
    "1"
])

# Write input to a temp file and run
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(genpartnew_input.rstrip() + "\n")  # ensure trailing newline
    temp_in = f.name

m.run_genpartnew(temp_in)
