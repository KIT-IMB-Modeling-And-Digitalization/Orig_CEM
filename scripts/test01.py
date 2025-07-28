import subprocess
import os
import shutil

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))       # scripts/
bin_path = os.path.join(base_dir, "..", "bin")              # bin/
src_path = os.path.join(base_dir, "..", "src")              # src/
scripts_results = os.path.join(base_dir, "results")         # scripts/results/
disrealnew_path = os.path.abspath(os.path.join(bin_path, "disrealnew"))
output_log = os.path.join(src_path, "disrealnew.out")       # save output in src/

# === Copy required input files to src/ ===
required_files = [
    "cement140w04flocf.img",
    "pcem140w04floc.img"
]

for filename in required_files:
    src_file = os.path.join(scripts_results, filename)  # from scripts/results/
    dest_file = os.path.join(src_path, filename)        # to src/
    if not os.path.exists(dest_file):
        try:
            shutil.copy2(src_file, dest_file)
            print(f"[DEBUG] Copied: {filename} → src/")
        except Exception as e:
            print(f"⚠️ Failed to copy {filename}: {e}")

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

# === Print and run ===
for idx, line in enumerate(disrealnew_input.strip().split("\n"), 1):
    print(f"[INPUT {idx:02}] {line}")

run_executable(disrealnew_path, disrealnew_input, output_log, cwd=src_path)
