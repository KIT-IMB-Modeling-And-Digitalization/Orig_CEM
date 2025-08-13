from setuptools import setup
from setuptools.command.build_ext import Command
import os
import subprocess
import platform
import shutil

class BuildExecutables(Command):
    def run(self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir  = os.path.join(base_dir, 'scripts')
        bin_dir  = os.path.join(base_dir, 'bin')
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
            print("✅ All Linux executables compiled successfully.")

        # ---------- W I N D O W S   B U I L D ----------
        elif is_windows:
            # give only disrealnew a larger stack (32 MiB)
            stack_flags = {"disrealnew": "/STACK:33554432"}   # 0x02000000 bytes

            for exe in executables:
                src_file  = os.path.join(src_dir, f"{exe}.c").replace("\\", "/")
                build_dir = os.path.join(base_dir, f"build_{exe}")
                bin_cmake = bin_dir.replace("\\", "/")

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
                        # Group all individual build folders into one 'build' directory
            final_build_dir = os.path.join(base_dir, "build")
            os.makedirs(final_build_dir, exist_ok=True)
            for exe in executables:
                src = os.path.join(base_dir, f"build_{exe}")
                dst = os.path.join(final_build_dir, f"build_{exe}")
                if os.path.exists(dst):
                    shutil.rmtree(dst)
                shutil.move(src, dst)

            print("✅ All Windows executables compiled with Visual Studio + CMake.")

        else:
            raise RuntimeError(f"⚠️ Unsupported platform: {platform.system()}")
    #TODO: create ./bin folder and copy files to folder ./bin,         shutil.copy(join(fpath,'lib','libiphreeqc.so'),join(current_path,'src','IPhreeqcPy','lib'))

setup(
    name        = "cement_sim",
    version     = "0.4",
    description = "Compile C files from scripts to bin (Linux & Windows)",
    cmdclass={
            'build_exe': BuildExecutables,
        },
    packages    = [],
    py_modules  = [],
    #TODO: aadd package data here ./bin file
    #mkdir(), rmdir/shutil.rm #delete folder #copy
)

#TODO: two step compilation 
#python setup.py build_exe
#pip install . (copies files to the site-packages folder) where other python packages are !!


