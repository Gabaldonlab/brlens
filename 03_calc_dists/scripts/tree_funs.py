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
import numpy as np
from scipy import stats
import statistics as st
import ete3

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
                 root_dict, prot_dict, distfile, sumfile):
        self.phylome_id = phylome_id
        self.seed_id = seed_id
        self.model = model
        self.lnl = lnl
        self.tree = ete3.PhyloTree(tree, sp_naming_function=get_species_tag)
        self.root_dict = root_dict
        self.prot_dict = prot_dict
        self.distfile = distfile
        self.sumfile = sumfile

    def root(self):
        '''
        Root the tree according to an outgroup
        '''

        if any(sp in self.root_dict for sp in self.tree.get_species()):
            ogdval = max([self.root_dict.get(sp, 0) for sp in self.tree.get_species()])
            ogsps = [k for k, val in self.root_dict.items() if val == ogdval][0]
            self.ogseq = [seq for seq in self.tree.get_leaf_names() if ogsps in seq][0]
        else:
            self.ogseq = self.tree.get_farthest_leaf()[0].get_leaf_names()[0]

        self.tree.set_outgroup(self.ogseq)

    def get_refpars(self):
        mphsptrees = dict()
        mphgtrees = dict()
        for i, subtree in enumerate(self.tree.traverse()):
            sps = [sp.split('_')[1] for sp in subtree.get_leaf_names()]
            if (len(sps) == len(set(sps)) and len(sps) > 10 and
                    subtree.get_farthest_leaf()[1] != 0):
                mphsptrees[len(sps)] = subtree
            elif len(sps) > 10 and subtree.get_farthest_leaf()[1] != 0:
                mphgtrees[len(sps) - len(set(sps))] = subtree

        if len(mphsptrees) > 0:
            self.rtree = mphsptrees.get(max(mphsptrees.keys()))
        elif len(mphgtrees) > 0:
            self.rtree = mphgtrees.get(min(mphgtrees.keys()))
        else:
            self.rtree = self.tree

        root = self.rtree.get_common_ancestor(self.rtree)
        self.rwdth = self.rtree.get_farthest_leaf()[1]

        rtldist = list()
        for leaf in self.rtree.get_leaf_names():
            rtldist.append(self.rtree.get_distance(root, leaf))

        self.rmean = np.mean(rtldist)
        self.rmed = np.median(rtldist)
        self.rskew = stats.skew(rtldist)
        self.rkurt = stats.kurtosis(rtldist)
        self.sd = st.stdev(rtldist)

        self.rstats = [self.phylome_id, self.prot_dict.get(self.seed_id, 'NA'),
                       self.rwdth, self.rmean, self.rmed, self.rskew,
                       self.rkurt, self.sd]

    def get_dists(self, seqsl):
        '''
        Calculate the distances between the
        '''
        distl = list()
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
                seqdistl = [self.phylome_id, self.seed_id,
                            self.prot_dict.get(self.seed_id, 'NA'), seq,
                            og_dist, seed_dist, og_dist / self.rmed,
                            seed_dist / self.rmed]
                distl.append(seqdistl)

        self.distl = distl

    def write_ofiles(self):
        for item in self.distl:
            self.distfile.write('\t'.join(str(v) for v in item) + '\n')
        self.sumfile.write('\t'.join(str(v) for v in self.rstats) + '\n')

    def run(self):
        '''
        Run all the pipe to obtain the distances
        '''

        self.seqsl = self.tree.get_leaf_names()
        if ('Phy' in self.seed_id and
                len(self.seqsl) == len(set(self.seqsl)) and
                len(self.seqsl) > 10):
            self.root()
            self.get_refpars()
            self.get_dists(self.seqsl)
            self.write_ofiles()
        else:
            print(('Tree cannot be parsed, seed or more than one sequences ',
                   'equally named.'))


# def main():
    # Read tree
    # start = time.time()
#     file_path = '03_calc_dists/data/0003_best_trees.txt'
#     phylome_id = file_path.rsplit('/', 1)[1].split('_', 1)[0]
#     phylome_no = int(file_path.rsplit('/', 1)[1].split('_', 1)[0])
#     # end = time.time()
#
#     # print('Time of reading the phylome: %s' % (end - start))
#
#     # Read protein file and get its protein
#     # start = time.time()
#     prot_dict = csv_to_dict('03_calc_dists/data/0003_all_protein_names.txt',
#                             sep='\t')
#     # end = time.time()
#
#     # print('Time of reading the protein file: %s' % (end - start))
#
#     dist_file = open('03_calc_dists/outputs/dist_ofile.txt', 'w')
#     dist_file.write('phylome_id\tseed\tprot_id\tdist_seq\tog_dist\tseed_dist\tog_ndist\tseed_ndist\n')
#     gene_sum = open('03_calc_dists/outputs/sum_ofile.txt', 'w')
#     gene_sum.write('phylome_id\tprot\twidth\tmean\tmedian\tskew\tkurt\tsd\n')
#     for tree in open(file_path):
#         # print(tree)
#         tree = tree.split('\t')
#         tree1 = phylome_tree(phylome_id, tree[0], tree[1], tree[2], tree[3],
#                              ROOTED_PHYLOMES[phylome_no], prot_dict,
#                              dist_file, gene_sum)
#
#         # start = time.time()
#         tree1.run()
#         # end = time.time()
#
#         # print('Time of running the tree functions: %s' % (end - start))
#
#     dist_file.close()
#     gene_sum.close()
#
#
# if __name__ == '__main__':
#     main()
