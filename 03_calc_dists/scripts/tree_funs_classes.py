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
import os
import sys
import ete3
from rooted_phylomes import ROOTED_PHYLOMES

import time

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *


# Definitions ----
def get_species_tag(node):
    return node.split("_")[1]


class phylome_tree(object):
    '''
    Object that allows to store the tree and its data
    '''

    def __init__(self, phylome_id, seed_id, model, lnl, tree,
                 root_dict, prot_dict):
        self.phylome_id = phylome_id
        self.seed_id = seed_id
        self.model = model
        self.lnl = lnl
        self.tree = ete3.PhyloTree(tree, sp_naming_function=get_species_tag)
        self.root_dict = root_dict
        self.prot_dict = prot_dict

    def root(self):
        '''
        Root the tree according to an outgroup
        '''

        og = self.tree.get_farthest_oldest_leaf(self.root_dict)
        self.ogseq = og.get_leaf_names()[0]
        self.tree.set_outgroup(og)

    def get_dists(self, seqsl):
        '''
        Calculate the distances between the
        '''
        self.ostr = ''
        for seq in self.seqsl:
            if seq != self.seed_id:
                seed_dist = self.tree.get_distance(seq, self.seed_id)
            else:
                seed_dist = 'NA'

            if seq != self.ogseq:
                og_dist = self.tree.get_distance(seq, self.ogseq)
            else:
                og_dist = 'NA'

            if seed_dist != 'NA' and og_dist != 'NA':
                self.ostr += ('%s\t%s\t%s\t%s\t%s\t%s\n' %
                              (self.phylome_id, self.seed_id,
                               self.prot_dict[self.seed_id], seq,
                               og_dist, seed_dist))

    def run(self):
        '''
        Run all the pipe to obtain the distances
        '''

        self.seqsl = self.tree.get_leaf_names()
        if 'Phy' in self.seed_id and len(self.seqsl) == len(set(self.seqsl)):
            self.root()
            # self.get_prot()
            self.get_dists(self.seqsl)
        else:
            print(('Tree cannot be parsed, seed or more than one sequences ',
                   'equally named.'))


def main():
    # Read tree
    start = time.time()
    file_path = '03_calc_dists/data/0003_best_trees.txt'
    file = open(file_path).read()
    phylome_id = file_path.rsplit('/', 1)[1].split('_', 1)[0]
    phylome_no = int(file_path.rsplit('/', 1)[1].split('_', 1)[0])
    end = time.time()

    print('Time of reading the phylome: %s' % (end - start))

    # Read protein file and get its protein
    start = time.time()
    prot_dict = csv_to_dict('03_calc_dists/data/0003_all_protein_names.txt',
                            sep='\t')
    end = time.time()

    print('Time of reading the protein file: %s' % (end - start))

    tree = file.split('\n')[0].split('\t')
    tree1 = phylome_tree(phylome_id, tree[0], tree[1], tree[2], tree[3],
                         ROOTED_PHYLOMES[phylome_no], prot_dict)

    start = time.time()
    tree1.run()
    end = time.time()

    print('Time of running the tree functions: %s' % (end - start))
    print(tree1.ostr)


if __name__ == '__main__':
    main()
