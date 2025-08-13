from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from pathlib import Path
import subprocess, shutil, os, platform, sys

# ---- settings you might tweak ------------------------------------------------
SRC_DIR = Path("csrc")                 # C sources live here
TARGET_NAME = "genpartnew"             # the executable we want to build
PACKAGE = "package_test"               # your package name under src/
BIN_DIR_IN_PKG = "_bin"                # where we ship binaries inside the package

def compile_executable(build_temp: Path, exe_out: Path):
    """
    Compile all .c files in csrc/ into a single executable.
    Assumes your flat-build works: gcc *.c -o genpartnew (plus -lm on POSIX).
    """
    build_temp.mkdir(parents=True, exist_ok=True)
    sources = sorted(str(p) for p in SRC_DIR.glob("*.c"))
    if not sources:
        raise RuntimeError(f"No C sources found in {SRC_DIR}")

    system = platform.system()
    if system == "Windows":
        # Prefer Microsoft cl.exe if available; otherwise try gcc (MinGW).
        cl = shutil.which("cl.exe")
        if cl:
            # MSVC: compile & link in one step
            # /Fe sets exe name, /Ox optimize; adjust flags as needed
            cmd = ["cl.exe", "/nologo", "/Ox", "/Fe:" + str(exe_out)] + sources
            env = os.environ.copy()
            # Ensure output dir exists; cl will drop .obj in cwd => use build_temp
            subprocess.check_call(cmd, cwd=str(build_temp), env=env)
        else:
            gcc = shutil.which("gcc") or shutil.which("clang")
            if not gcc:
                raise RuntimeError("No compiler found (need cl.exe or gcc/clang).")
            cmd = [gcc, "-O3", "-o", str(exe_out)] + sources
            subprocess.check_call(cmd, cwd=str(build_temp))
    else:
        # POSIX: gcc/clang with -lm is common
        cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
        if not cc:
            raise RuntimeError("No C compiler found (need cc/gcc/clang).")
        cmd = [cc, "-O3", "-o", str(exe_out)] + sources + ["-lm"]
        subprocess.check_call(cmd, cwd=str(build_temp))

class build_py(_build_py):
    def run(self):
        # 1) run normal Python build first
        super().run()

        # 2) build the executable into the build_lib package dir
        bin_dir = Path(self.build_lib) / PACKAGE / BIN_DIR_IN_PKG
        bin_dir.mkdir(parents=True, exist_ok=True)

        exe_name = TARGET_NAME + (".exe" if platform.system() == "Windows" else "")
        exe_out = bin_dir / exe_name

        # Use a temporary build dir to hold intermediates/objects
        build_temp = Path(self.build_lib) / "_build_bins"
        compile_executable(build_temp, exe_out)

        # Make sure it’s executable on POSIX
        if os.name == "posix":
            exe_out.chmod(exe_out.stat().st_mode | 0o111)

        print(f"[build_py] built {exe_out}")

setup(
    name="package-test",
    version="0.0.1",
    package_dir={"": "src"},              # src/ layout
    packages=find_packages(where="src"),
    python_requires=">=3.9",
    include_package_data=True,
    package_data={"package_test": ["_bin/*"]},  # ship whatever ends up in _bin
    cmdclass={"build_py": build_py},            # <-- hook our custom build step
)
