import subprocess
import os

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))       # scripts/
bin_path = os.path.join(base_dir, "..", "bin")              # bin/
src_path = os.path.join(base_dir, "..", "src")              # src/
results_dir = os.path.join(base_dir, "results")             # results/
disrealnew_path = os.path.abspath(os.path.join(bin_path, "disrealnew"))
output_log = os.path.join(results_dir, "disrealnew.out")

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# === Symlink required input files to results directory ===
input_files = {
    "cement140w04flocf.img": os.path.join(results_dir, "cement140w04flocf.img"),  # already in results/
    "pcem140w04floc.img": os.path.join(src_path, "pcem140w04floc.img"),           # from src/
}

for filename, source_path in input_files.items():
    link_path = os.path.join(results_dir, filename)
    if not os.path.exists(link_path):
        try:
            os.symlink(source_path, link_path)
            print(f"[DEBUG] Symlinked: {filename} → {source_path}")
        except Exception as e:
            print(f"⚠️ Failed to symlink {filename}: {e}")

# === Prepare input data ===
disrealnew_input = "\n".join([
    "-2794",
    "cement140w04flocf.img",  # Input microstructure
    "1 2 3 4 5 6 7 28 26",
    "35",
    "pcem140w04floc.img",     # Input particle IDs
    "44990", "1",
    "5850", "2",
    "8692", "3",
    "2631", "4",
    "1100", "5",
    "2062", "6",
    "839", "7",
    "0",
    "1000", "0", "500",
    "0.0001 9000.",
    "0.01 9000.",
    "0.00002 10000.",
    "0.002 2500.",
    "50", "5", "5000", "100",
    "0.00", "20.0", "20.0", "0.0",
    "40.0", "83.14", "80.0", "0.00035",
    "0.72", "0", "0", "1", "10",
    "1.0", "0", "0", "1"
]) + "\n"

# === Function to run executable ===
def run_executable(executable, input_data, output_file, cwd=None):
    print(f"[DEBUG] Running: {executable}")
    print(f"[DEBUG] CWD: {cwd}")
    process = subprocess.run(
        [executable],
        input=input_data,
        text=True,
        capture_output=True,
        cwd=cwd
    )
    with open(output_file, "w") as f:
        f.write(process.stdout)
    if process.returncode != 0:
        print(f"⚠️ Non-zero exit code ({process.returncode}) for {executable}")
        print("STDOUT (trimmed):", process.stdout[:500])
        print("STDERR:", process.stderr)
    else:
        print(f"✅ Execution completed for {executable}, output saved to {output_file}")

# === Run the executable ===
for idx, line in enumerate(disrealnew_input.strip().split("\n"), 1):
    print(f"[INPUT {idx:02}] {line}")

run_executable(disrealnew_path, disrealnew_input, output_log, cwd=results_dir)


