#! /usr/bin/env python

import sys
import os
import platform
import random
import multiprocessing

WORK_FORCE = multiprocessing.Queue()
GLOBAL_RNG = random.Random()
PLATFORM = platform.system().lower()

PACKAGE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
BASE_DIR = os.path.abspath(os.path.dirname(PACKAGE_DIR))
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
BIN_DIR = None
if PLATFORM == 'linux':
    BIN_DIR = os.path.join(BASE_DIR, "bin", "linux", "bin")
elif PLATFORM == 'darwin':
    BIN_DIR = os.path.join(base_dir, "bin", "mac", "bin")
elif PLATFORM == 'windows':
    BIN_DIR = os.path.join(base_dir, "bin", "win", "bin")

