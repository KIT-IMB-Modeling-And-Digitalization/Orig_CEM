# utils.py
# This file is part of pycemhyd3d to detect the OS
import platform

def is_windows() -> bool:
    """
    Return True if running on Windows, False otherwise.
    """
    return platform.system() == "Windows"
