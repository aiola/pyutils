#!/usr/bin/env python
# I/O utilities

import os

def find_file(path, file_name):
    """ Find file in a given path
    """
    for root, _, files in os.walk(path):
        for current_file in files:
            if file == file_name:
                yield os.path.join(root, current_file)
