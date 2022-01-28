#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
calc_dists.py -- Phylome distances calculations

The program gets the phylome and calculate the distances between
the tree seed to each leaf and the tree outgroup to each leaf.
It also gets the speciation and duplication events at the lineage
branches of the species between each pair of sequences. It
returns a csv table with all calculations.

Requirements: ETE3, pandas, multiprocessing

Written by Moisès Bernabeu moigil.bernabeu@gmail.com
Month 2022
'''

# Import libraries ----
import sys
import os
from optparse import OptionParser
from rooted_phylomes import ROOTED_PHYLOMES as root_dict
import ete3
import pandas as pd
from multiprocessing import Process, Manager

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *


# Definitions ----
def get_species_tag(node):
    return node.split("_")[1]


def root(tree, root_dict):
    '''
    Root the tree according to an outgroup
    '''

    if any(sp in root_dict for sp in tree.get_species()):
        ogdval = max([root_dict.get(sp, 0) for sp in tree.get_species()])
        ogsps = [k for k, val in root_dict.items() if val == ogdval and k in tree.get_species()][0]
        ogseq = [seq for seq in tree.get_leaf_names() if ogsps in seq][0]
    else:
        ogseq = tree.get_farthest_leaf()[0].get_leaf_names()[0]

    tree.set_outgroup(ogseq)

    return ogseq


def get_events(tree, leaf, seqfrom):
    events = dict()
    events['S'] = 0
    events['D'] = 0
    tree.get_descendant_evol_events()
    ltstr = tree.get_common_ancestor(seqfrom, leaf)
    ltstreeln = ltstr.get_leaf_names()

    for node in ltstr.get_leaves_by_name(seqfrom)[0].get_ancestors():
        nln = node.get_leaf_names()
        if len(nln) <= len(ltstreeln):
            events[node.evoltype] += 1
    for node in ltstr.get_leaves_by_name(leaf)[0].get_ancestors():
        nln = node.get_leaf_names()
        if len(nln) <= len(ltstreeln):
            events[node.evoltype] += 1
    events[ltstr.evoltype] -= 1

    return events


def get_dists(tree, leaf, seed_id, ogseq, phylome_id, prot_dict):
    '''
    Calculate the distances between the
    '''

    if leaf != seed_id:
        seed_dist = tree.get_distance(leaf, seed_id)
        seed_events = get_events(tree, leaf, seed_id)
    else:
        seed_dist = 'NA'

    if leaf != ogseq:
        og_dist = tree.get_distance(leaf, ogseq)
        og_events = get_events(tree, leaf, ogseq)
    else:
        og_dist = 'NA'

    if seed_dist != 'NA' and og_dist != 'NA':
        leafdistd = dict()
        leafdistd['id'] = phylome_id
        leafdistd['seed'] = seed_id
        leafdistd['prot'] = prot_dict.get(seed_id, 'NA')
        leafdistd['leaf'] = leaf
        leafdistd['og_dist'] = og_dist
        leafdistd['seed_dist'] = seed_dist
        leafdistd['seed_sp'] = seed_events['S']
        leafdistd['seed_dupl'] = seed_events['D']
        leafdistd['og_sp'] = og_events['S']
        leafdistd['og_dupl'] = og_events['D']

        return leafdistd
    else:
        return None


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, prot_dict, olist):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.prot_dict = prot_dict
        self.olist = olist

    def run(self):
        tree = self.tree_row.split('\t')
        t = ete3.PhyloTree(tree[3], sp_naming_function=get_species_tag)

        if (len(t.get_species()) > 10 and
                len(t.get_leaf_names()) < 3 * len(t.get_species())):
            print('Calculating: %s, species no.: %s, leaves no.: %s' %
                  (tree[0], len(t.get_species()), len(t.get_leaf_names())))
            og = root(t, root_dict[int(self.phylome_id)])

            for leaf in t.get_leaf_names():
                leaf_dist = get_dists(t, leaf, tree[0], og,
                                      self.phylome_id, self.prot_dict)
                if leaf_dist is not None:
                    self.olist.append(leaf_dist)


def main():
    # Script options definition ----
    parser = OptionParser()
    parser.add_option('-d', '--def', dest='default',
                      help='Default configurations to execute with our data.',
                      action='store_true')
    parser.add_option('-f', '--file', dest='ifile',
                      help='In file',
                      metavar='<path/to/file.txt>')
    parser.add_option('-o', '--output', dest='odir',
                      help='Output directory',
                      metavar='<path/to/output>')
    parser.add_option('-p', '--prot', dest='prot',
                      help='File with protein codes',
                      metavar='<path/to/file.txt>')
    parser.add_option('-c', '--cpu', dest='cpus',
                      help='File with protein codes',
                      metavar='<path/to/file.txt>')
    (options, args) = parser.parse_args()

    if options.default:
        ifile = '../data/0469_best_trees.txt'
        prots = '../data/0469_all_protein_names.txt'
        odir = '../outputs'
        cpus = 4
    else:
        ifile = options.ifile
        odir = options.odir
        prots = options.prot
        cpus = int(options.cpus)

    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]

    dist_fn = '/'.join([odir, (phylome_id + '_dist.tsv')])

    if not file_exists(dist_fn):
        print('Creating: ', dist_fn)

        create_folder(odir)

        trees = open(ifile, 'r').read().split('\n')
        prot_dict = csv_to_dict(prots, '\t')

        with Manager() as manager:
            olist = manager.list()

            processes = list()
            for tree_row in trees:
                if tree_row != '':
                    if len(processes) >= cpus:
                        done = False
                        while not done:
                            for process in processes:
                                if not process.is_alive():
                                    processes.remove(process)
                                    done = True

                    process = dist_process(tree_row, phylome_id,
                                           prot_dict, olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            odf = pd.DataFrame(list(olist))
            odf.to_csv(dist_fn, index=False)


if __name__ == '__main__':
    main()
