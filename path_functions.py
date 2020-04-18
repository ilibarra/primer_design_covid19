'''
Created on May 23, 2016

@author: ignacio
'''

# define basic import

import os
from os.path import join, exists, isdir, basename, abspath, dirname, getmtime
from os import remove, listdir, mkdir, rmdir, makedirs, system, chdir
from shutil import copy2, move


# the size of a given path
def filesize(p):
    return os.path.getsize(p)
def get_total_size(paths):
    return sum([get_file_size(p) for p in paths])
def get_file_size(p):
    return os.path.getsize(p)

def get_parent_path(p):
    return os.path.abspath(os.path.join(os.path.dirname(p), ".."))