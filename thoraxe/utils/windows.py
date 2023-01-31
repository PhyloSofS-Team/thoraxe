"""
Windows: Functions useful for testing and Windows support.
"""

import platform

def is_windows():
    """Return True in Windows."""
    current_platform = platform.system()
    print(f"Current platform: {current_platform}")
    name = current_platform.lower()
    return name[0:3] == 'win'
