#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
seq2seq_norm.py -- Sequence to sequence distances and normalisation

The script gets the sequence to sequence distances and calculates the
normalisation factors. It does not calculate the normalised distance.

Requirements: ete3, pandas, normalisation

Written by Moisès Bernabeu <moigil.bernabeu@gmail.com>
February 2022
'''

# Import libraries ----
import sys
import os
from optparse import OptionParser
import ete3
import pandas as pd
from multiprocessing import Process, Manager
from normalisation import clade_norm

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import file_exists, create_folder


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

    Args:rbls_ref
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


def get_dists(tree, from_seq, to_seq, seed_id, phylome_id, normfdic):
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

    Returns:
      dict: main features and distances of the tree

    Raises:
      Exception: description
    '''

    dist = tree.get_distance(from_seq, to_seq)
    events = get_events(tree, from_seq, to_seq)

    leafdistd = dict()
    leafdistd['id'] = phylome_id
    leafdistd['tree'] = seed_id
    leafdistd['from'] = from_seq
    leafdistd['from_sp'] = get_species_tag(from_seq)
    leafdistd['to'] = to_seq
    leafdistd['to_sp'] = get_species_tag(to_seq)
    leafdistd['sp'] = events['S']
    leafdistd['dupl'] = events['D']
    leafdistd['mrca_type'] = events['MRCA']
    leafdistd['dist'] = dist

    for item in normfdic:
        leafdistd['ndist_%s' % item['Normalising group']] = dist / item['mrca_to_tip_median']

    return leafdistd


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, gnmdf, olist):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.gnmdf = gnmdf
        self.olist = olist

    def run(self):
        if '\t' in self.tree_row:
            tree = self.tree_row.split('\t')
            t = ete3.PhyloTree(tree[3], sp_naming_function=get_species_tag)
            tname = tree[0]
        else:
            t = ete3.PhyloTree(self.tree_row, sp_naming_function=get_species_tag)
            tname = 'sp'

        if (len(t.get_species()) > 10 and
                len(t.get_leaf_names()) < 3 * len(t.get_species())):
            print('Calculating: %s, species no.: %s, leaves no.: %s' %
                  (tname, len(t.get_species()), len(t.get_leaf_names())))

            tnames = t.get_leaf_names()

            t.set_outgroup(t.get_midpoint_outgroup())
            t.get_descendant_evol_events()

            normfdic = clade_norm(t, tname, self.gnmdf,
                                  'Normalising group')

            for i, from_seq in enumerate(tnames):
                for to_seq in tnames[i + 1:]:
                    if from_seq != to_seq:
                        leaf_dist = get_dists(t, from_seq, to_seq, tname,
                                              self.phylome_id, normfdic)
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
    parser.add_option('-p', '--pinfo', dest='pinfo',
                      help='Phylome information in tsv format.',
                      metavar='<path/to/file.tsv>')
    parser.add_option('-c', '--cpu', dest='cpus',
                      help='File with protein codes',
                      metavar='<path/to/file.txt>', type='int')
    (options, args) = parser.parse_args()

    if options.default:
        ifile = '../data/0739_0.txt'
        odir = '../outputs'
        gnmdf = '../data/README'
        cpus = 4
    else:
        ifile = options.ifile
        odir = options.odir
        gnmdf = options.pinfo
        cpus = options.cpus

    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]
    file_id = ifile.rsplit('/', 1)[1].split('.', 1)[0]

    dist_fn = '/'.join([odir, (file_id + '_dist.csv')])
    norm_fn = '/'.join([odir, (file_id + '_norm.csv')])
    if not file_exists(dist_fn):
        print('Creating: ', dist_fn)

        create_folder(odir)

        trees = open(ifile, 'r').read().split('\n')
        gnmdf = pd.read_csv(gnmdf)

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

                    process = dist_process(tree_row, phylome_id, gnmdf, olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(dist_fn, index=False)


if __name__ == '__main__':
    main()
