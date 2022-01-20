#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
calc_dists.py -- Run in parallel the distances calculations.

Requirements: ete3

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import sys
import os
from optparse import OptionParser
import glob
from tree_thread import thread

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *

# Script options definition ----
# parser = OptionParser()
# parser.add_option('-d', '--def', dest='default',
#                   help='Default configurations to execute with our data.',
#                   action='store_true')
# parser.add_option('-f', '--file', dest='ifile',
#                   help='In file',
#                   metavar='<path/to/file.txt>')
# (options, args) = parser.parse_args()


# Definitions ----
def yet_calculated(dir, name):
    '''
    Check the existence and size of the file, returns True or False
    '''
    fp = '/'.join([dir, name])

    # Check file existence
    if os.path.isfile(fp):
        # Check the file is not empty
        if os.stat(fp).st_size != 0:
            return True
    else:
        return False


def main():
    files = glob.glob('../data/*_best_trees.txt')
    protfiles = glob.glob('../data/*_all_protein_names.txt')
    workdir = '../outputs'
    threads = 2

    tasks = list()

    for i, file in enumerate(files):
        phylome_id = file.rsplit('/', 1)[1].split('_', 1)[0]
        distfilenm = phylome_id + '_dist.tsv'

        if file != '':
            if (len(tasks) >= threads):
                # Wait for a process to finish
                done = False
                while not done:
                    for task in tasks:
                        if not task.is_alive():
                            # With these conditions the thread is ended, print
                            # and free a slot
                            task.show()
                            tasks.remove(task)
                            done = True

            # Download trees
            if not yet_calculated(workdir, distfilenm):
                tree_thread = thread(file, workdir, protfiles[i])
                tree_thread.show()
                tree_thread.start()
                tasks.append(tree_thread)

    for task in tasks:
        task.join()
        task.show()


if __name__ == '__main__':
    main()
