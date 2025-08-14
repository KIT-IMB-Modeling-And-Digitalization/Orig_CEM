import package_test.module as m
from tempfile import NamedTemporaryFile

distrib3d_input = "\n".join([
    "-99",
    "cem140w04floc.img",      # original microstructure image
    "cement140",              # correlation filters root name
    "cement140w04flocf.img",  # output final microstructure
    "0.7344",
    "0.6869",
    "0.0938",
    "0.1337",
    "0.1311",
    "0.1386",
    "0.0407",
    "0.0408"
])

# Write to temp file
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(distrib3d_input.rstrip() + "\n")  # ensure trailing newline
    temp_in = f.name

# Run distrib3d
m.run_distrib3d(temp_in)
