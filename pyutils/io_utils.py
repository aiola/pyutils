""" This module contains I/O utilities.
"""

import os

def find_file(path, file_name):
    """ Find file in a given path
    """
    for root, _, files in os.walk(path):
        for current_file in files:
            if current_file == file_name:
                yield os.path.join(root, current_file)
