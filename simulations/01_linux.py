import subprocess
import os
import shutil
from linux_func import *

id = "05"

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))          # → simulations/
parent_dir = os.path.abspath(os.path.join(base_dir, ".."))    # → Smart03 or Github_version
src = os.path.join(parent_dir, "scripts")
dst = os.path.join(parent_dir, f"scripts_{id}")

 # === Copy the entire folder ===
if os.path.exists(dst):
    shutil.rmtree(dst)
shutil.copytree(src, dst)
src_path = os.path.abspath(os.path.join(base_dir,"..", "scripts_%s" % id))


#base_dir = os.path.dirname(os.path.abspath(__file__))  # simulations/
bin_path = os.path.join(base_dir, "..", "bin")
#src_path = os.path.join(base_dir, "..", "scripts")              # scripts/
genpartnew_path = os.path.join(bin_path, "genpartnew")  
distrib3d_path = os.path.join(bin_path, "distrib3d")
disrealnew_path = os.path.join(bin_path, "disrealnew")

# === Input for genpartnew ===
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
    os.path.join(src_path, "cem140w04floc.img"), # Filename to save image to
    os.path.join(src_path, "pcem140w04floc.img"), # Filename to save particle IDs to
    "1" # Menu selection to exit program (end)
])



# === Input for distrib3d ===
distrib3d_input = "\n".join([
    "-99", # Random number seed
    os.path.join(src_path, "cem140w04floc.img"), # Filename of original microstructure image
    "cement140",  # File root name for correlation filters (sil, c3s, c3a, etc.)
    os.path.join(src_path, "cement140w04flocf.img"), # Filename under which to save final microstructure
    "0.7344",  # Volume fraction of C3S in microstructure to be created
    "0.6869",  # Surface area fraction of C3S in microstructure to be created
    "0.0938",  # Volume fraction of C2S in microstructure to be created
    "0.1337",  # Surface area fraction of C2S in microstructure to be created
    "0.1311",  # Volume fraction of C3A in microstructure to be created
    "0.1386",  # Surface area fraction of C3A in microstructure to be created
    "0.0407",  # Volume fraction of C4AF in microstructure to be created
    "0.0408"   # Surface area fraction of C4AF in microstructure to be created
]) + "\n"



# === Input for disrealnew ===
# === Use simple filenames ===
# img1 = "cement140w04flocf.img"
# img2 = "pcem140w04floc.img"

img1 = os.path.join(src_path, "cem140w04floc.img")
img2 = os.path.join(src_path, "pcem140w04floc.img")


# === Prepare input data ===
disrealnew_input = "\n".join([
    "-2794", # Random number seed
    img1, 
    # Filename containing input 3-D phase ID microstructure
    "1 2 3 4 5 6 7 28 26",  # Phase assignments for C3S, C2S, etc..
    "35",  # Phase ID for C3A in fly ash particles
    img2, # Filename containing input 3-D particle ID microstructure
    "44990",  # Number of one-pixel particles to add
    "1",  # Add one-pixel particles of phase C3S
    "5850",  # Number of one-pixel particles to add
    "2",  # Add one-pixel particles of phase C2S
    "8692",  # Number of one-pixel particles to add
    "3",  # Add one-pixel particles of phase C3A
    "2631",  # Number of one-pixel particles to add
    "4",  # Add one-pixel particles of phase C4AF
    "1100",  # Number of one-pixel particles to add
    "5",  # Add one-pixel particles of gypsum
    "2062",  # Number of one-pixel particles to add
    "6",  # Add one-pixel particles of hemihydrate
    "839",  # Number of one-pixel particles to add
    "7",  # Add one-pixel particles of anhydrite
    "0",  # Number of one-pixel particles to add
    "1000",  # Number of cycles of hydration model to execute
    "0",  # Flag for executing model under saturated conditions
    "500",  # Maximum number of diffusion steps per cycle
    "0.0001 9000.",  # Nucleation parameters for CH
    "0.01 9000.",  # Nucleation parameters for gypsum (dihydrate)
    "0.00002 10000.",  # Nucleation parameters for C3AH6
    "0.002 2500.",  # Nucleation parameters for FH3
    "50",  # Cycle frequency for checking pore space percolation
    "5",  # Cycle frequency for checking total solids percolation
    "5000",  # Cycle frequency for outputting particle hydration stats
    "100",  # Cycle frequency for outputting hydrated microstructures
    "0.00",  # Induction time in hours
    "20.0",  # Initial specimen temperature in degrees Celsius
    "20.0",  # Ambient temperature in degrees Celsius
    "0.0",  # Overall heat transfer coefficient
    "40.0",  # Activation energy for cement hydration
    "83.14",  # Activation energy for pozzolanic reactions
    "80.0",  # Activation energy for slag hydration
    "0.00035",  # Conversion factor to go from cycles to real time
    "0.72",  # Aggregate volume fraction in actual concrete mixture
    "0",  # Flag indicating that hydration is under isothermal conditions
    "0",  # Flag indicating that C-S-H conversion is prohibited
    "1",  # Flag indicating that CH/aggregate precipitation is allowed
    "10",  # Number of slices to include in a hydration movie file
    "1.0",  # One-voxel dissolution bias factor
    "0",  # Number of cycles to execute before resaturation
    "0",  # Flag indicating that C-S-H morphology is random
    "1"   # Flag indicating that pH does influence hydration kinetics
]) + "\n"


def main():
    results_dir = os.path.join(base_dir, "results_%s" % id)
    os.makedirs(results_dir, exist_ok=True)

    # Output text file for stdout
    output_path = os.path.join(src_path, "genpartnew.out")
    output_log_distrib3d = os.path.join(src_path, "distrib3d.out")
    output_log_disrealnew = os.path.join(src_path, "disrealnew.out")


    before_files = set(os.listdir(src_path))

    run_executable_genpartnew(genpartnew_path, genpartnew_input, output_path)

    run_executable_distrib3d(distrib3d_path, distrib3d_input, output_log_distrib3d, cwd=src_path)
    
    run_executable_disrealnew(disrealnew_path, disrealnew_input, output_log_disrealnew, cwd=src_path)
    # #=== Detect and move new output files to results ===
    after_files = set(os.listdir(src_path))
    new_files = after_files - before_files

    for fname in new_files:
        src = os.path.join(src_path, fname)
        dst = os.path.join(results_dir, fname)
        try:
            shutil.move(src, dst)
            print(f"[DEBUG] Moved output: {fname} → results/")
        except Exception as e:
            print(f"⚠️ Failed to move {fname}: {e}")
    try:
        shutil.rmtree(src_path)
        print(f"[DEBUG] Removed temporary folder: {src_path}")
    except Exception as e:
        print(f"⚠️ Failed to remove {src_path}: {e}")

if __name__ == "__main__":
    main()
