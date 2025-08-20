import cement_sim.module as m
from tempfile import NamedTemporaryFile
from pathlib import Path
import platform
import os, uuid, subprocess, shutil


#1ST STEP OF THE USER: CHOOSE THE run ID TO LABEL OUTPUTS AND RESULT FOLDER
id = "01"   

########### Creation of results folder ###########
script_dir = Path(__file__).parent.resolve()
results_dir = script_dir / f"result_{id}"
results_dir.mkdir(parents=True, exist_ok=True)

# Installed _bin (where the executable + assets live)
bin_dir = Path(getattr(m, "_BIN_DIR", Path(m.__file__).parent / "_bin")).resolve()

# Per-run folder inside _bin (short + unique)
run_id  = "gp_" + uuid.uuid4().hex[:6]
run_dir = bin_dir / run_id
run_dir.mkdir(parents=False, exist_ok=False)

# Find the OS and correspondingly the executables
def _is_windows():
    return platform.system() == "Windows"

if _is_windows():
    exe = bin_dir / "genpartnew.exe"
    exe2 = bin_dir / "distrib3d.exe"
    exe3 = bin_dir / "disrealnew.exe"
else:
    exe = bin_dir / "genpartnew"
    exe2 = bin_dir / "distrib3d"
    exe3 = bin_dir / "disrealnew"


# 2ND STEP OF THE USER: DEFINE INPUTS FOR EACH EXECUTABLE

# Use ONLY short filenames in stdin (we'll run with cwd=run_dir)
genpartnew_input = "\n".join([
    "-3034",  # Random number seed
    "2",  # Menu choice to place particles
    "16",  # Number of size classes to place
    "0",  # Dispersion distance between particles in pixels
    "0.0604",  # Calcium sulfate (total) volume fraction
    "0.515 0.041",  # Fractions of calcium sulfate that are hemihydrate and anhydrite
    "1", "17", "1",  # Size class 1: Number, Radius, Phase ID (1 = cement)
    "1", "15", "1",  # Size class 2
    "1", "14", "1",  # Size class 3
    "1", "13", "1",  # Size class 4
    "2", "12", "1",  # Size class 5
    "2", "11", "1",  # Size class 6
    "4", "10", "1",  # Size class 7
    "5", "9", "1",  # Size class 8
    "8", "8", "1",  # Size class 9
    "13", "7", "1",  # Size class 10
    "21", "6", "1",  # Size class 11
    "38", "5", "1",  # Size class 12
    "73", "4", "1",  # Size class 13
    "174", "3", "1",  # Size class 14
    "450", "2", "1",  # Size class 15
    "2674", "1", "1",  # Size class 16
    "4",  # Menu selection to report phase counts
    "3",  # Menu selection to flocculate particles
    "1",  # Number of separate flocs (particle clusters) to create
    "8",  # Menu selection to output current microstructure to file
    f"cem140w04floc_{id}.img",      # short name to save image to
    f"pcem140w04floc_{id}.img",     # short name to save particle IDs to
    "1"   # Menu selection to exit program (end)
])

in_name  = f"{run_id}/cem140w04floc_{id}.img"
out_name = f"{run_id}/cement140w04flocf_{id}.img"

distrib3d_input = "\n".join([
    "-99",              # Random number seed
    in_name,            # Filename of original microstructure image
    "cement140",        # File root name for correlation filters (sil, c3s, c3a, etc.)
    out_name,           # Filename under which to save final microstructure
    "0.7344",  # Volume fraction of C3S in microstructure to be created
    "0.6869",  # Surface area fraction of C3S in microstructure to be created
    "0.0938",  # Volume fraction of C2S in microstructure to be created
    "0.1337",  # Surface area fraction of C2S in microstructure to be created
    "0.1311",  # Volume fraction of C3A in microstructure to be created
    "0.1386",  # Surface area fraction of C3A in microstructure to be created
    "0.0407",  # Volume fraction of C4AF in microstructure to be created
    "0.0408"   # Surface area fraction of C4AF in microstructure to be created
])

phase_name = f"cement140w04flocf_{id}.img"
part_name  = f"pcem140w04floc_{id}.img"

in_phase = f"{run_id}/{phase_name}"
in_part  = f"{run_id}/{part_name}"

disrealnew_input = "\n".join([
    "-2794",                    # Random seed number
    in_phase,                   # Filename containing input 3-D phase ID microstructure
    "1 2 3 4 5 6 7 28 26",      # Phase assignments for C3S, C2S, etc.
    "35",                       # Phase ID for C3A in fly ash particles
    in_part,                    # Filename containing input 3-D particle ID microstructure
    "44990","1",                # Number of one-pixel particles to add of phase C3S
    "5850","2",                 # Number of one-pixel particles to add of phase C2S
    "8692","3",                 # Number of one-pixel particles to add of phase C3A
    "2631","4",                 # Number of one-pixel particles to add of phase C4AF
    "1100","5",                 # Number of one-pixel particles to add of phase gypsum
    "2062","6",                 # Number of one-pixel particles to add of phase hemihydrate
    "839","7",                  # Number of one-pixel particles to add of phase anhydrite
    "0",                        # Number of one-pixel particles to add
    "1000",                     # Number of cycles of hydration model to execute
    "0",                        # Flag for executing model under saturated conditions
    "500",                      # Maximum number of diffusion steps per cycle
    "0.0001 9000.",            # Nucleation parameters for CH
    "0.01 9000.",              # Nucleation parameters for gypsum (dihydrate)
    "0.00002 10000.",          # Nucleation parameters for C3AH6
    "0.002 2500.",             # Nucleation parameters for FH3
    "50",                      # Cycle frequency for checking pore space percolation
    "5",                       # Cycle frequency for checking total solids percolation
    "5000",                    # Cycle frequency for outputting particle hydration stats
    "100",                     # Cycle frequency for outputting hydrated microstructures
    "0.00",                    # Induction time in hours
    "20.0",                    # Initial specimen temperature in degrees Celsius
    "20.0",                    # Ambient temperature in degrees Celsius
    "0.0",                     # Overall heat transfer coefficient
    "40.0",                    # Activation energy for cement hydration
    "83.14",                   # Activation energy for pozzolanic reactions
    "80.0",                    # Activation energy for slag hydration
    "0.00035",                 # Conversion factor to go from cycles to real time
    "0.72",                    # Aggregate volume fraction in actual concrete mixture
    "0",                       # Flag indicating that hydration is under isothermal conditions
    "0",                       # Flag indicating that C-S-H conversion is prohibited
    "1",                       # Flag indicating that CH/aggregate precipitation is allowed
    "10",                      # Number of slices to include in a hydration movie file
    "1.0",                     # One-voxel dissolution bias factor
    "0",                       # Number of cycles to execute before resaturation
    "0",                       # Flag indicating that C-S-H morphology is random
    "1"                        # Flag indicating that pH does influence hydration kinetics
]) + "\n"  # ensure trailing newline




# =========================
# Append: run genpartnew now
# =========================

# Write input to a temp file
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(genpartnew_input.rstrip() + "\n")
    temp_in = f.name

# Run inside the per-run folder; capture all terminal output to genpartnew.out


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

# Write temp stdin and run distrib3d with cwd=_bin; log to distrib3d.out in run_dir
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(distrib3d_input.rstrip() + "\n")
    temp_in2 = f.name


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
# Append: run disrealnew now 
# =========================

phase_path = run_dir / phase_name
part_path  = run_dir / part_name
if not phase_path.exists():
    raise FileNotFoundError(f"[disrealnew] Missing phase microstructure in sandbox: {phase_path}")
if not part_path.exists():
    raise FileNotFoundError(f"[disrealnew] Missing particle-ID microstructure in sandbox: {part_path}")

# IMPORTANT: run disrealnew with cwd = _bin so it can find its runtime data there.
# Use paths prefixed with the per-run folder so reads/writes land in run_dir.




# Snapshot files in _bin before running (to catch anything disrealnew drops into _bin root)
pre_files = {p.name for p in bin_dir.iterdir() if p.is_file()}

# Write temp stdin and run with cwd=_bin; log to disrealnew.out inside run_dir
with NamedTemporaryFile('w+', delete=False) as f:
    f.write(disrealnew_input)
    temp_in3 = f.name
    
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
