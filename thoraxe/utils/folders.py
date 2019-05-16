"""
folders: Util functions to create output folder structure.
"""

import os


def create_subfolder(folder, subfolder):
    """
    Return path to the subfolder. It creates the folder if it doesn't exist.
    """
    subfolder_path = os.path.join(folder, subfolder)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path, exist_ok=True)  # Python 3.2+
    return subfolder_path
