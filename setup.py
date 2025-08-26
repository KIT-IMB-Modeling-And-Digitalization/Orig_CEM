'''setup.py for pycemhyd3d

Author:
    Omid Jahromi <omid.esmaeelipoor@gmail.com>

Overview:
    Custom setuptools script for the pycemhyd3d package.

    Responsibilities:
        - Defines a `build_exe` command that compiles the bundled C sources
          (cemhyd3d executables: genpartnew, distrib3d, disrealnew).
        - Handles both Linux (gcc) and Windows (CMake + Visual Studio) builds.
        - Copies compiled executables into the package tree
          (`src/pycemhyd3d/cempy3d/`) so they are shipped with the wheel/sdist.
        - Cleans up intermediate `bin/` and `build/` folders after compilation.

    Note:
        This script extends setuptools with a platform-aware build step, so
        users can install pycemhyd3d and have ready-to-use native executables
        across supported platforms.
'''


from setuptools import setup
from setuptools import Command
import os
import subprocess
import platform
import shutil

class BuildExecutables(Command):
    user_options = []
    def initialize_options(self): pass
    def finalize_options(self): pass
    def run(self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir = os.path.join(base_dir, 'src', 'pycemhyd3d', 'cempy3d')
        if not os.path.isdir(src_dir):
            raise RuntimeError(f"Source dir not found: {src_dir}")

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
    # Create _bin directory inside src/pycemhyd3d
        if is_linux:
            pkg_dir = os.path.join(base_dir, 'src', 'pycemhyd3d')
            bin_dir_pkg = os.path.join(pkg_dir, 'cempy3d')
            os.makedirs(bin_dir_pkg, exist_ok=True)

            for exe in executables:
                compiled_path = os.path.join(bin_dir, exe)
                if os.path.isfile(compiled_path):
                    shutil.copy(compiled_path, bin_dir_pkg)

            if os.path.exists(bin_dir):
                shutil.rmtree(bin_dir)

            # Copy everything from scripts into _bin so executables can run in-place
            for root, _, files in os.walk(src_dir):
                rel_path = os.path.relpath(root, src_dir)
                dest_dir = os.path.join(bin_dir_pkg, rel_path) if rel_path != '.' else bin_dir_pkg
                os.makedirs(dest_dir, exist_ok=True)
                for f in files:
                    shutil.copy(os.path.join(root, f), os.path.join(dest_dir, f))

        # ---------- Copy compiled Windows executables into src/pycemhyd3d/cempy3d ----------
        if is_windows:
            pkg_dir = os.path.join(base_dir, 'src', 'pycemhyd3d')
            bin_dir_pkg = os.path.join(pkg_dir, 'cempy3d')
            os.makedirs(bin_dir_pkg, exist_ok=True)

            for exe in executables:
                # check both bin\Release\exe.exe and bin\exe.exe
                exe_name = exe + ".exe"
                release_path = os.path.join(bin_dir, "Release", exe_name)
                flat_path    = os.path.join(bin_dir, exe_name)

                if os.path.isfile(release_path):
                    shutil.copy(release_path, bin_dir_pkg)
                elif os.path.isfile(flat_path):
                    shutil.copy(flat_path, bin_dir_pkg)

            # Copy everything from scripts/ into _bin so executables can run in-place
            for root, _, files in os.walk(src_dir):
                rel_path = os.path.relpath(root, src_dir)
                dest_dir = os.path.join(bin_dir_pkg, rel_path) if rel_path != '.' else bin_dir_pkg
                os.makedirs(dest_dir, exist_ok=True)
                for f in files:
                    shutil.copy(os.path.join(root, f), os.path.join(dest_dir, f))


setup(
    name        = "pycemhyd3d",
    version     = "0.4",
    description = "Compile C files from scripts to bin (Linux & Windows)",
    cmdclass={
            'build_exe': BuildExecutables,
        },
    packages=['pycemhyd3d'],
    package_dir={'': 'src'},
    package_data={
        'pycemhyd3d': ['cempy3d/*']
    },
    include_package_data=True,
    py_modules  = [],
)



