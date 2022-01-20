#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
tree_funs.py -- Define the functions to manage trees

Here are defined the functions to root, calculate interest distances and
retrieve data from them.

Requirements: ete3

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import sys
import os
import ete3

from pprint import pprint
from rooted_phylomes import ROOTED_PHYLOMES


# Definitions ----
def get_species_tag (node):
    return node.split("_")[1]


file = open('03_calc_dists/data/0003_best_trees.txt').read()

# For line
tree = file.split('\n')[0].split('\t')
t = ete3.PhyloTree(tree[-1], sp_naming_function=get_species_tag)
