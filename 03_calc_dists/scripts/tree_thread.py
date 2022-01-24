#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
tree_thread.py -- Construct the threads to perform distances calculations

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import os
import sys
import threading
from tree_funs import phylome_tree
from rooted_phylomes import ROOTED_PHYLOMES

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *


# Definitions ----
class thread(threading.Thread):
    '''
    Run tree calculations in a thread
    '''

    # Class constructor
    def __init__(self, ifile, workdir, prot_file):
        threading.Thread.__init__(self)
        self.ifile = ifile
        self.pdict = prot_file
        self.workdir = workdir

        self.phylome_id = self.ifile.rsplit('/', 1)[1].split('_', 1)[0]
        self.phylome_no = int(self.ifile.rsplit('/', 1)[1].split('_', 1)[0])

        self.message = 'thread for the item %s is created' % self.phylome_id

    # Show method
    def show(self):
        '''
        Prints thread output
        '''

        print('Calculating: %s' % self.message)

    # Run method
    def run(self):
        '''
        Run the thread
        '''

        create_folder(self.workdir)

        self.message = 'Thread for the item %s is running' % self.phylome_id

        prot_dict = csv_to_dict(self.pdict, sep='\t')

        dist_file = open('/'.join([self.workdir,
                                   (self.phylome_id + '_dist.tsv')]), 'w')
        dist_file.write(('phylome_id\tseed\tprot_id\tdist_seq\tog_dist' +
                         '\tseed_dist\tog_ndist\tseed_ndist\tseed_sps' +
                         '\tseed_dupls\tog_sps\tog_dupls\n'))
        gene_sum = open('/'.join([self.workdir,
                                  (self.phylome_id + '_sum.tsv')]), 'w')
        gene_sum.write(('phylome_id\tprot\twidth\tmean\tmedian\tskew' +
                        '\tkurt\tsd\n'))
        for tree in open(self.ifile):
            tree = tree.split('\t')
            tree1 = phylome_tree(self.phylome_id, tree[0], tree[1], tree[2],
                                 tree[3], ROOTED_PHYLOMES[self.phylome_no],
                                 prot_dict, dist_file, gene_sum)
            tree1.run()
        dist_file.close()
        gene_sum.close()
