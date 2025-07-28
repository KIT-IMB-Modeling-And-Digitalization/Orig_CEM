from setuptools import setup
from setuptools.command.build_ext import build_ext
import os
import subprocess

class BuildExecutables(build_ext):
    def run(self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir = os.path.join(base_dir, 'src')
        bin_dir = os.path.join(base_dir, 'bin')
        os.makedirs(bin_dir, exist_ok=True)

        # Define file-specific flags
        flags = {
            "genpartnew": "-lm",
            "distrib3d": "-lm",
            "disrealnew": "-lm -O2"
        }

        executables = ["genpartnew", "distrib3d", "disrealnew"]
        all_success = True

        for exe in executables:
            src = os.path.join(src_dir, f"{exe}.c")
            out = os.path.join(bin_dir, exe)
            compile_flags = flags.get(exe, "")
            print(f"Compiling {src} -> {out} with flags: {compile_flags}")
            result = subprocess.run(
                f"gcc {src} -o {out} {compile_flags}",
                shell=True
            )
            if result.returncode != 0:
                print(f"❌ Failed to compile {src}")
                all_success = False
                break

        if all_success:
            print("✅ All executables compiled successfully.")

setup(
    name="cement_sim",
    version="0.2",
    description="Compile C files from src to bin",
    cmdclass={'build_ext': BuildExecutables},
    packages=[],
    py_modules=[],
)
