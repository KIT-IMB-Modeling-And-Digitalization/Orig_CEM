from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import os
import subprocess
import platform
import shutil

PKG_NAME = "cement_sim"
PKG_SRC_DIR = os.path.join(os.path.dirname(__file__), "src", PKG_NAME)
PKG_BIN_DIR = os.path.join(PKG_SRC_DIR, "bin")  # executables go here

class BuildExecutables(build_ext):
    def run(self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir  = os.path.join(base_dir, 'scripts')              # your C sources stay here
        bin_dir  = PKG_BIN_DIR                                    # <-- build into the package
        os.makedirs(bin_dir, exist_ok=True)

        is_windows = platform.system() == "Windows"
        is_linux   = platform.system() == "Linux"

        executables = ["genpartnew", "distrib3d", "disrealnew"]

        # ---------- L I N U X   B U I L D ----------
        if is_linux:
            flags = {
                "genpartnew":  "-lm",
                "distrib3d":   "-lm",
                "disrealnew":  "-lm -O2",
            }
            for exe in executables:
                src  = os.path.join(src_dir, f"{exe}.c")
                out  = os.path.join(bin_dir, exe)
                opts = flags.get(exe, "")
                print(f"[Linux] Compiling {src} -> {out} with flags: {opts}")
                if subprocess.run(f"gcc {src} -o {out} {opts}", shell=True).returncode:
                    raise RuntimeError(f"❌ Failed to compile {src}")
            print("✅ All Linux executables compiled successfully into package bin/.")

        # ---------- W I N D O W S   B U I L D ----------
        elif is_windows:
            # give only disrealnew a larger stack (32 MiB)
            stack_flags = {"disrealnew": "/STACK:33554432"}   # 0x02000000 bytes

            for exe in executables:
                src_file  = os.path.join(src_dir, f"{exe}.c").replace("\\", "/")
                build_dir = os.path.join(base_dir, f"build_{exe}")
                bin_cmake = PKG_BIN_DIR.replace("\\", "/")

                if os.path.exists(build_dir):
                    shutil.rmtree(build_dir)
                os.makedirs(build_dir)

                link_flag = stack_flags.get(exe, "")

                with open(os.path.join(build_dir, "CMakeLists.txt"), "w") as f:
                    f.write(f"""
cmake_minimum_required(VERSION 3.10)
project({exe} C)

add_executable({exe} "{src_file}")

# Ensure a larger stack when requested (compatible with old CMake)
set_target_properties({exe} PROPERTIES
    LINK_FLAGS "{link_flag}"
    RUNTIME_OUTPUT_DIRECTORY "{bin_cmake}"
)
""")

                print(f"[Windows] Configuring {exe} …")
                subprocess.run(
                    'cmake -G "Visual Studio 17 2022" -A x64 .',
                    cwd=build_dir, shell=True, check=True
                )
                print(f"[Windows] Building {exe} …")
                subprocess.run(
                    'cmake --build . --config Release',
                    cwd=build_dir, shell=True, check=True
                )

            print("✅ All Windows executables compiled into package bin/.")

        else:
            raise RuntimeError(f"⚠️ Unsupported platform: {platform.system()}")

        # Continue with the normal build_ext pipeline
        super().run()

setup(
    name        = "cement-sim",
    version     = "0.5.0",
    description = "Cement simulation tools packaged with three native executables",
    cmdclass    = {"build_ext": BuildExecutables},

    # --- packaging: src-layout with a real package ---
    package_dir = {"": "src"},
    packages    = find_packages("src"),
    include_package_data = True,

    # include the built executables inside the wheel/install
    package_data = {
        PKG_NAME: ["bin/*"],      # ship the 3 compiled binaries
    },

    # expose convenient CLI commands for users
    entry_points = {
        "console_scripts": [
            "genpartnew=cement_sim.cli:run_genpartnew",
            "distrib3d=cement_sim.cli:run_distrib3d",
            "disrealnew=cement_sim.cli:run_disrealnew",
        ]
    },

    python_requires = ">=3.8",
)
