import subprocess
import os

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))  # scripts/
bin_path = os.path.join(base_dir, "..", "bin")
genpartnew_path = os.path.join(bin_path, "genpartnew")

# Create results directory inside scripts/
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

# Output text file for stdout
output_path = os.path.join(results_dir, "genpartnew.out")

def run_executable(executable, input_data, output_file):
    """Runs an executable with input data and redirects output."""
    process = subprocess.run(
        [executable],
        input=input_data,
        text=True,
        capture_output=True
    )
    if process.returncode != 0:
        print(f"❌ Execution failed for {executable}")
        print("Error:", process.stderr)
        exit(1)
    with open(output_file, "w") as f:
        f.write(process.stdout)
    print(f"✅ Execution completed for {executable}, output saved to {output_file}")


# === Input for genpartnew ===
genpartnew_input = "\n".join([
    "-3034", "2", "16", "0", "0.0604", "0.515 0.041",
    "1", "17", "1", "1", "15", "1", "1", "14", "1", "1", "13", "1", "2", "12", "1",
    "2", "11", "1", "4", "10", "1", "5", "9", "1", "8", "8", "1", "13", "7", "1",
    "21", "6", "1", "38", "5", "1", "73", "4", "1", "174", "3", "1", "450", "2", "1",
    "2674", "1", "1", "4", "3", "1", "8",
    os.path.join("results", "cem140w04floc.img"),
    os.path.join("results", "pcem140w04floc.img"),
    "1"
])

run_executable(genpartnew_path, genpartnew_input, output_path)
