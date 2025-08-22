# utils.py
import platform

def is_windows() -> bool:
    """
    Return True if running on Windows, False otherwise.
    """
    return platform.system() == "Windows"
