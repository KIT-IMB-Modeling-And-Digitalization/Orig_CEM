import subprocess

def compile_c_file(c_file, output_executable, flags="-lm"):
    """Compiles a C file into an executable."""
    compile_cmd = f"gcc {c_file} -o {output_executable} {flags}"
    result = subprocess.run(compile_cmd, shell=True)
    if result.returncode != 0:
        print(f"Compilation failed for {c_file}")
        exit(1)

def run_executable(executable, input_data, output_file):
    """Runs an executable with input data and redirects output."""
    process = subprocess.run(
        f"./{executable}", input=input_data, text=True, capture_output=True, shell=True
    )
    if process.returncode != 0:
        print(f"Execution failed for {executable}")
        print("Error:", process.stderr)
        exit(1)
    with open(output_file, "w") as f:
        f.write(process.stdout)
    print(f"Execution completed for {executable}, output saved to {output_file}")

# Step 1: Compile genpartnew
compile_c_file("genpartnew.c", "genpartnew")

# Step 2: Define input values for genpartnew
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
    "cem140w04floc.img",  # Filename to save image to
    "pcem140w04floc.img",  # Filename to save particle IDs to
    "1"   # Menu selection to exit program (end)
])
run_executable("genpartnew", genpartnew_input, "genpartnew.out")

# Step 3: Compile distrib3d
compile_c_file("distrib3d.c", "distrib3d")

# Step 4: Define input values for distrib3d
distrib3d_input = "\n".join([
    "-99",  # Random number seed
    "cem140w04floc.img",  # Filename of original microstructure image
    "cement140",  # File root name for correlation filters (sil, c3s, c3a, etc.)
    "cement140w04flocf.img",  # Filename under which to save final microstructure
    "0.7344",  # Volume fraction of C3S in microstructure to be created
    "0.6869",  # Surface area fraction of C3S in microstructure to be created
    "0.0938",  # Volume fraction of C2S in microstructure to be created
    "0.1337",  # Surface area fraction of C2S in microstructure to be created
    "0.1311",  # Volume fraction of C3A in microstructure to be created
    "0.1386",  # Surface area fraction of C3A in microstructure to be created
    "0.0407",  # Volume fraction of C4AF in microstructure to be created
    "0.0408"   # Surface area fraction of C4AF in microstructure to be created
])
run_executable("distrib3d", distrib3d_input, "distrib3d.out")

# Step 5: Compile disrealnew
compile_c_file("disrealnew.c", "disrealnew", "-lm -O2")

# Step 6: Define input values for disrealnew
disrealnew_input = "\n".join([
    "-2794",  # Random number seed
    "cement140w04flocf.img",  # Filename containing input 3-D phase ID microstructure
    "1 2 3 4 5 6 7 28 26",  # Phase assignments for C3S, C2S, etc.
    "35",  # Phase ID for C3A in fly ash particles
    "pcem140w04floc.img",  # Filename containing input 3-D particle ID microstructure
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
])
run_executable("disrealnew", disrealnew_input, "disrealnew.out")

print("All executions completed successfully!")
