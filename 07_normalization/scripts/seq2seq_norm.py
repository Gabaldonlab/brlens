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
from rooted_phylomes import ROOTED_PHYLOMES as root_dict
import ete3
import pandas as pd
from multiprocessing import Process, Manager
from normalisation import subtree_tt_ref, mrca_tt_ref, root_tt_ref

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *


# Definitions ----
def get_species_tag(node):
    if '_' in node:
        return node.split("_")[1]
    else:
        return node


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
    events['MRCA'] = ltstr.evoltype

    return events


def get_dists(tree, from_seq, to_seq, seed_id, phylome_id, prot_dict):
    '''
    Retrieves distances between pairs of sequences

    The function gets a set of two sequences and calculates the distances
    between them and associates the tree with some information. Retrieves
    a dictionary with the main features and distances of the tree.

    Args:
      tree (PhyloTree): phylogenetic tree imported with ete3
      from_seq (char): string with the from leaf name
      to_seq (char): string with the to leaf name
      seed_id (char): name of the seed sequence
      phylome_id (char): code of the phylome in PhylomeDB
      prot_dict (dictionary): contains the protein id linked to the tree

    Returns:
      dict: main features and distances of the tree

    Raises:
      Exception: description

    '''

    root(tree, root_dict[int(phylome_id)])
    tree.get_descendant_evol_events()

    st_ref = subtree_tt_ref(tree)
    # mrca_ref = mrca_tt_ref(tree)
    root_ref = root_tt_ref(tree)

    dist = tree.get_distance(from_seq, to_seq)
    events = get_events(tree, from_seq, to_seq)

    leafdistd = dict()
    leafdistd['id'] = phylome_id
    leafdistd['tree'] = seed_id
    leafdistd['prot'] = prot_dict.get(seed_id, 'NA')
    leafdistd['from'] = from_seq
    leafdistd['from_sp'] = get_species_tag(from_seq)
    leafdistd['to'] = to_seq
    leafdistd['to_sp'] = get_species_tag(to_seq)
    leafdistd['dist'] = dist
    leafdistd['dist_norm_st'] = dist / st_ref['med']
    leafdistd['st_median'] = st_ref['med']
    # leafdistd['dist_norm_mrca'] = dist / mrca_ref['med']
    # leafdistd['mrca_median'] = mrca_ref['med']
    # leafdistd['dist_norm_root'] = dist / root_ref['med']
    # leafdistd['root_median'] = root_ref['med']
    # leafdistd['dist_norm_width'] = dist / root_ref['twdth']
    # leafdistd['root_median'] = root_ref['twdth']
    # leafdistd['dist_norm_rbls'] = dist / root_ref['rwdth']
    # leafdistd['rbls'] = root_ref['rwdth']
    leafdistd['sp'] = events['S']
    leafdistd['dupl'] = events['D']
    leafdistd['mrca_type'] = events['MRCA']

    return leafdistd


def get_sp_dist(tree, from_seq, to_seq):
    tree.get_descendant_evol_events()

    mrca_ref = mrca_tt_ref(tree)
    root_ref = root_tt_ref(tree)

    dist = tree.get_distance(from_seq, to_seq)

    leafdistd = dict()
    leafdistd['from_sp'] = from_seq
    leafdistd['to_sp'] = to_seq
    leafdistd['dist'] = dist
    leafdistd['dist_norm_mrca'] = dist / mrca_ref['med']
    leafdistd['mrca_median'] = mrca_ref['med']
    leafdistd['dist_norm_root'] = dist / root_ref['med']
    leafdistd['root_median'] = root_ref['med']
    leafdistd['dist_norm_width'] = dist / root_ref['twdth']
    leafdistd['root_median'] = root_ref['twdth']
    leafdistd['dist_norm_rbls'] = dist / root_ref['rwdth']
    leafdistd['rbls'] = root_ref['rwdth']

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

            tnames = t.get_leaf_names()

            for i, from_seq in enumerate(tnames):
                for to_seq in tnames[i + 1:]:
                    if from_seq != to_seq:
                        leaf_dist = get_dists(t, from_seq, to_seq, tree[0],
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
                      metavar='<path/to/file.txt>',
                      default=None)
    parser.add_option('-c', '--cpu', dest='cpus',
                      help='File with protein codes',
                      metavar='<path/to/file.txt>')
    parser.add_option('-s', '--sp', dest='sp_tree',
                      help='This option only requires a tree, not ' +
                           'proteins file.',
                      action='store_true', default=False)
    (options, args) = parser.parse_args()

    if options.default:
        ifile = '../data/0435_best_trees.txt'
        prots = '../data/0435_all_protein_names.txt'
        odir = '../outputs'
        cpus = 6
        sp_tree = False
    else:
        ifile = options.ifile
        odir = options.odir
        prots = options.prot
        cpus = int(options.cpus)
        sp_tree = options.sp_tree

    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]
    file_id = ifile.rsplit('/', 1)[1].split('.', 1)[0]

    if sp_tree:
        dist_fn = '/'.join([odir, (file_id + '_sptree_dist.csv')])

        if not file_exists(dist_fn):
            create_folder(odir)
            t = ete3.PhyloTree(ifile)

            tnames = t.get_leaf_names()

            olist = list()
            for i, from_seq in enumerate(tnames):
                for to_seq in tnames[i + 1:]:
                    if from_seq != to_seq:
                        olist.append(get_sp_dist(t, from_seq, to_seq))

            odf = pd.DataFrame(olist)
            odf.to_csv(dist_fn, index=False)

    else:
        dist_fn = '/'.join([odir, (file_id + '_dist.csv')])
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

                # Writing output files
                odf = pd.DataFrame(list(olist))
                odf.to_csv(dist_fn, index=False)


if __name__ == '__main__':
    main()
