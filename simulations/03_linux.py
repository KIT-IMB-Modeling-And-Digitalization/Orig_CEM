import subprocess
import os
import shutil
from linux_func import *

id = "03"

# === Setup paths ===
base_dir = os.path.dirname(os.path.abspath(__file__))          # → simulations/
# parent_dir = os.path.abspath(os.path.join(base_dir, ".."))    # → Smart03 or Github_version
# src = os.path.join(parent_dir, "scripts")
# dst = os.path.join(parent_dir, f"scripts")

# # === Copy the entire folder ===
# if os.path.exists(dst):
#     shutil.rmtree(dst)
# shutil.copytree(src, dst)
src_path = os.path.abspath(os.path.join(base_dir,"..", "scripts_%s" % id))


#base_dir = os.path.dirname(os.path.abspath(__file__))  # simulations/
bin_path = os.path.join(base_dir, "..", "bin")
#src_path = os.path.join(base_dir, "..", "scripts")              # scripts/
genpartnew_path = os.path.join(bin_path, "genpartnew")  
distrib3d_path = os.path.join(bin_path, "distrib3d")
disrealnew_path = os.path.join(bin_path, "disrealnew")

# === Input for genpartnew ===
genpartnew_input = "\n".join([
    "-3034", "2", "16", "0", "0.0604", "0.515 0.041",
    "1", "17", "1", "1", "15", "1", "1", "14", "1", "1", "13", "1", "2", "12", "1",
    "2", "11", "1", "4", "10", "1", "5", "9", "1", "8", "8", "1", "13", "7", "1",
    "21", "6", "1", "38", "5", "1", "73", "4", "1", "174", "3", "1", "450", "2", "1",
    "2674", "1", "1", "4", "3", "1", "8",
    os.path.join(src_path, "cem140w04floc.img"),
    os.path.join(src_path, "pcem140w04floc.img"),
    "1"
])

# Run the genpartnew executable


# === Input for distrib3d ===
distrib3d_input = "\n".join([
    "-99",
    os.path.join(src_path, "cem140w04floc.img"),
    "cement140",
    os.path.join(src_path, "cement140w04flocf.img"),
    "0.7344", "0.6869", "0.0938", "0.1337",
    "0.1311", "0.1386", "0.0407", "0.0408"
]) + "\n"



# === Input for disrealnew ===
# === Use simple filenames ===
# img1 = "cement140w04flocf.img"
# img2 = "pcem140w04floc.img"

img1 = os.path.join(src_path, "cem140w04floc.img")
img2 = os.path.join(src_path, "pcem140w04floc.img")


# === Prepare input data ===
disrealnew_input = "\n".join([
    "-2794",
    img1,
    "1 2 3 4 5 6 7 28 26",
    "35",
    img2,
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


if __name__ == "__main__":
    main()