'''Utility helpers for pycemhyd3d.

Author:
    Omid Jahromi <omid.esmaeelipoor@gmail.com>

Overview:
    Cross-platform helpers used by the pipeline and executors.
'''
import platform

def is_windows() -> bool:
    '''Return True if the current platform is Windows, else False.

    Uses:
        Helps select executable names (with/without `.exe`) and choose
        platform-specific behaviors during execution.
    '''
    return platform.system() == "Windows"
