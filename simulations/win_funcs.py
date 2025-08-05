import subprocess
import os
import shutil

def run_executable_genpartnew(executable, input_data, output_file):
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

def run_executable_distrib3d(executable, input_data, output_file, cwd=None):
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

def run_executable_disrealnew(executable, input_data, output_file, cwd=None):
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