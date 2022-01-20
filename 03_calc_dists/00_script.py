#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
00_script.py -- Title of the script

Brief description

Requirements:

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import sys
import os
from optparse import OptionParser

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *

# Script options definition ----
parser = OptionParser()
parser.add_option('-d', '--def', dest='default',
                  help='Default configurations to execute with our data.',
                  action='store_true')
parser.add_option('-f', '--file', dest='ifile',
                  help='In file',
                  metavar='<path/to/file.txt>')
(options, args) = parser.parse_args()

# Definitions ----
