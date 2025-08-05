import subprocess
import os

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))           # scripts/
bin_path = os.path.join(base_dir, "..", "bin")                  # bin/
src_path = os.path.join(base_dir, "..", "scripts")                  # src/
results_dir = os.path.join(base_dir, "results")                 # scripts/results/
distrib3d_path = os.path.join(bin_path, "distrib3d")
output_log = os.path.join(results_dir, "distrib3d.out")

os.makedirs(results_dir, exist_ok=True)

def run_executable_with_cwd(executable, input_data, output_file, cwd=None):
    print(f"[DEBUG] Running: {executable}")
    print(f"[DEBUG] From cwd: {cwd}")
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

# === Input for distrib3d ===
distrib3d_input = "\n".join([
    "-99",
    "../simulations/results/cem140w04floc.img",
    "cement140",
    "../simulations/results/cement140w04flocf.img",
    "0.7344", "0.6869", "0.0938", "0.1337",
    "0.1311", "0.1386", "0.0407", "0.0408"
]) + "\n"

run_executable_with_cwd(
    distrib3d_path,
    distrib3d_input,
    output_log,
    cwd=src_path   # Run from inside src/ to find cement140.*
)
