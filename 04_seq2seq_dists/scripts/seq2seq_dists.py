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

    The tree is rooted with a dictionary containing species-to-age information,
    the farthest sequence from an species in the tree which has maximum age is
    selected to be the outgroup of the tree.

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function
        root_dict (dictionary): a dictionary containing the species age,
        indexes refers to the phylome.
        Eg.: 3: {'SP1': 1, 'SP2': 2}
        Where 3 is the phylome number and the SP2 is older than SP1.

    Returns:
        string: the outgroup sequence and the tree object is rooted without
        being returned.

    Raises:
        Exception: description
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
    '''
    Analyse the speciation and duplication events

    The function goes through the lineage branches and gets the number of
    speciation and duplication events between both sequences (leaf and
    seqfrom).

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function
        leaf (str): sequence name of the leaf.
        seqfrom (str): name of the reference sequence.

    Returns:
        dict: dictionary containing the number events of speciation ('S') and
        duplication ('D') events.

    Raises:
        Exception: description
    '''

    events = dict()
    events['S'] = 0
    events['D'] = 0
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


def get_dists(tree, from_seq, to_seq, seed_id, phylome_id, prot_dict):
    '''
    Calculate the distances between the
    '''

    dist = tree.get_distance(from_seq, to_seq)
    events = get_events(tree, from_seq, to_seq)

    leafdistd = dict()
    leafdistd['id'] = phylome_id
    leafdistd['tree'] = seed_id
    leafdistd['prot'] = prot_dict.get(seed_id, 'NA')
    leafdistd['from'] = seed_id
    leafdistd['from_sp'] = get_species_tag(from_seq)
    leafdistd['to'] = to_seq
    leafdistd['to_sp'] = get_species_tag(to_seq)
    leafdistd['dist'] = dist
    leafdistd['sp'] = events['S']
    leafdistd['dupl'] = events['D']

    return leafdistd


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
            root(t, root_dict[int(self.phylome_id)])
            t.get_descendant_evol_events()

            tnames = t.get_leaf_names()

            for i, from_sp in enumerate(tnames):
                for to_sp in tnames[i + 1:]:
                    if from_sp != to_sp:
                        leaf_dist = get_dists(t, from_sp, to_sp, tree[0],
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
        ifile = '../data/0435_best_trees.txt'
        prots = '../data/0435_all_protein_names.txt'
        odir = '../outputs'
        cpus = 4
    else:
        ifile = options.ifile
        odir = options.odir
        prots = options.prot
        cpus = int(options.cpus)

    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]

    dist_fn = '/'.join([odir, (phylome_id + '_dist.csv')])

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
