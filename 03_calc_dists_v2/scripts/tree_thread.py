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
    def __init__(self, itree, workdir, prot_file, phylome_id, dist_ofile,
                 sum_ofile):
        threading.Thread.__init__(self)
        self.itree = itree
        self.pdict = prot_file

        self.phylome_id = phylome_id
        self.phylome_no = int(phylome_id)

        self.dist_ofile = dist_ofile
        self.sum_ofile = sum_ofile

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

        self.message = 'Threads for the item %s is running' % self.phylome_id

        tree = phylome_tree(self.phylome_id, self.itree[0], self.itree[1],
                            self.itree[2], self.itree[3],
                            ROOTED_PHYLOMES[self.phylome_no], self.pdict,
                            self.dist_ofile, self.sum_ofile)
        tree.run()
