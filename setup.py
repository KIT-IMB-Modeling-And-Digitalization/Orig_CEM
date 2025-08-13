from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
try:
    from setuptools.command.build_editable import build_editable as _build_editable
except ImportError:
    _build_editable = None

from pathlib import Path
import subprocess, shutil, os, platform

# ---- settings --------------------------------------------------------------
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
            cmd = ["cl.exe", "/nologo", "/Ox", "/Fe:" + str(exe_out)] + sources
            env = os.environ.copy()
            subprocess.check_call(cmd, cwd=str(build_temp), env=env)
        else:
            gcc = shutil.which("gcc") or shutil.which("clang")
            if not gcc:
                raise RuntimeError("No compiler found (need cl.exe or gcc/clang).")
            cmd = [gcc, "-O3", "-o", str(exe_out)] + sources
            subprocess.check_call(cmd, cwd=str(build_temp))
    else:
        cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
        if not cc:
            raise RuntimeError("No C compiler found (need cc/gcc/clang).")
        cmd = [cc, "-O3", "-o", str(exe_out)] + sources + ["-lm"]
        subprocess.check_call(cmd, cwd=str(build_temp))

def _build_into(build_dir: Path):
    """Shared logic for both build_py and build_editable"""
    bin_dir = build_dir / PACKAGE / BIN_DIR_IN_PKG
    bin_dir.mkdir(parents=True, exist_ok=True)

    exe_name = TARGET_NAME + (".exe" if platform.system() == "Windows" else "")
    exe_out = bin_dir / exe_name

    build_temp = build_dir / "_build_bins"
    compile_executable(build_temp, exe_out)

    if os.name == "posix":
        exe_out.chmod(exe_out.stat().st_mode | 0o111)

    # Also copy to source tree so editable installs work
    src_bin_dir = Path("src") / PACKAGE / BIN_DIR_IN_PKG
    src_bin_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(exe_out, src_bin_dir / exe_name)

    print(f"[build] built {exe_out} and copied to {src_bin_dir / exe_name}")

class build_py(_build_py):
    def run(self):
        super().run()
        _build_into(Path(self.build_lib))

if _build_editable:
    class build_editable(_build_editable):
        def run(self):
            super().run()
            _build_into(Path(self.build_lib))
    cmdclass = {"build_py": build_py, "build_editable": build_editable}
else:
    cmdclass = {"build_py": build_py}

# ---- main setup() ----------------------------------------------------------

setup(
    name="package-test",
    version="0.0.1",
    package_dir={"": "src"},              # src/ layout
    packages=find_packages(where="src"),
    python_requires=">=3.9",
    include_package_data=True,
    package_data={"package_test": ["_bin/*"]},
    cmdclass=cmdclass,
)
