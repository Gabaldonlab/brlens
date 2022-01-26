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
parser = OptionParser()
parser.add_option('-d', '--def', dest='default',
                  help='Default configurations to execute with our data.',
                  action='store_true')
parser.add_option('-f', '--fold', dest='ifolder',
                  help='In files (in txt) folders',
                  metavar='<path/to/file.txt>')
parser.add_option('-w', '--wd', dest='workdir',
                  help='working directory',
                  metavar='<path/to/workdir>')
parser.add_option('-t', '--threads', dest='threads',
                  help='working directory',
                  metavar='<path/to/workdir>')
(options, args) = parser.parse_args()


# Definitions ----
def yet_calculated(name):
    '''
    Check the existence and size of the file, returns True or False
    '''

    # Check file existence
    if os.path.isfile(name):
        # Check the file is not empty
        if os.stat(name).st_size != 0:
            return True
    else:
        return False


def main():
    if options.default:
        files = glob.glob('../data/*_best_trees.txt')
        protfiles = glob.glob('../data/*_all_protein_names.txt')
        workdir = '../outputs'
        threads = 4
    else:
        files = glob.glob(options.ifolder + '*_best_trees.txt')
        protfiles = glob.glob(options.ifolder + '*_all_protein_names.txt')
        workdir = options.workdir
        threads = options.threads

    tasks = list()

    create_folder(workdir)

    for i, file in enumerate(files):
        phylome_id = file.rsplit('/', 1)[1].split('_', 1)[0]
        distfile = open(workdir + '/' + phylome_id + '_dist.tsv', 'w')
        sumfile = open(workdir + '/' + phylome_id + '_sum.tsv', 'w')

        prot_dict = csv_to_dict(protfiles[i], sep='\t')

        if not yet_calculated(distfile.name):
            distfile.write('phylome_id\tseed\tprot_id\tdist_seq\tog_dist' +
                           '\tseed_dist\tog_ndist\tseed_ndist\tseed_sps' +
                           '\tseed_dupls\tog_sps\tog_dupls\n')
            sumfile.write('phylome_id\tprot\twidth\tmean\tmedian\tskew' +
                          '\tkurt\tsd\n')

            for tree in open(file).read().split('\n'):
                tree = tree.split('\t')
                if file != '':
                    if (len(tasks) >= threads):
                        # Wait for a process to finish
                        done = False
                        while not done:
                            for task in tasks:
                                if not task.is_alive():
                                    # With these conditions the thread is ended,
                                    # print and free a slot
                                    # task.show()
                                    tasks.remove(task)
                                    done = True

                    tree_thread = thread(tree, workdir, prot_dict,
                                         phylome_id, distfile, sumfile)
                    # tree_thread.show()
                    tree_thread.start()
                    tasks.append(tree_thread)

                for task in tasks:
                    task.join()
                    task.show()

            distfile.close()
            sumfile.close()


if __name__ == '__main__':
    main()
