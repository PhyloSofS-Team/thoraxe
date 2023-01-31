"""
Windows: Functions useful for testing and Windows support.
"""

import platform

def is_windows():
    """Return True in Windows."""
    name = platform.system().lower()
    return name[0:3] == 'win'
