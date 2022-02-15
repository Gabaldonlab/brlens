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
import numpy as np
from scipy import stats
import statistics as st
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf


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
    events['MRCA'] = ltstr.evoltype

    return events


def get_dists(tree, from_seq, to_seq, seed_id, phylome_id,
              prot_dict, st_refs, mrca_refs):


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
    leafdistd['dist_norm_st'] = dist / st_refs['med']
    leafdistd['dist_norm_mrca'] = dist / mrca_refs['med']
    leafdistd['sp'] = events['S']
    leafdistd['dupl'] = events['D']
    leafdistd['mrca_type'] = events['MRCA']

    return leafdistd


def subtree_ref(phylome_id, prot_dict, seed_id, tree):
    '''
    Reference statistics from subtree

    Subtree is selected according to maximum number of tips and only speciation
    events.

    Args:
        phylome_id (string): id from phylome
        prot_dict (dict): species to age dictionary
        seed_id (string): seed sequence idea
        tree (Phylotree): phylome tree

    Returns:
        dict: set of basic descriptive statistics of the lengths in the subtree

    Raises:
        Exception: description
    '''

    mphsptrees = dict()
    mphgtrees = dict()
    for i, subtree in enumerate(tree.traverse()):
        sps = [sp.split('_')[1] for sp in subtree.get_leaf_names()]
        if (len(sps) == len(set(sps)) and len(sps) > 10 and
                subtree.get_farthest_leaf()[1] != 0):
            mphsptrees[len(sps)] = subtree
        elif len(sps) > 10 and subtree.get_farthest_leaf()[1] != 0:
            mphgtrees[len(sps) - len(set(sps))] = subtree

    if len(mphsptrees) > 0:
        rtree = mphsptrees.get(max(mphsptrees.keys()))
    elif len(mphgtrees) > 0:
        rtree = mphgtrees.get(min(mphgtrees.keys()))
    else:
        rtree = tree

    root = rtree.get_common_ancestor(rtree)
    wdth = rtree.get_farthest_leaf()[1]

    rtldist = list()
    for leaf in rtree.get_leaf_names():
        rtldist.append(rtree.get_distance(root, leaf))

    mean = np.mean(rtldist)
    med = np.median(rtldist)
    skew = stats.skew(rtldist)
    kurt = stats.kurtosis(rtldist)
    sd = st.stdev(rtldist)

    rstats = {'id': phylome_id, 'prot': prot_dict.get(seed_id, 'NA'),
              'wdth': wdth, 'mean': mean, 'med': med, 'skew': skew,
              'kurt': kurt, 'sd': sd}

    return rstats


def mrca_tt_ref(tree, seed_id):
    ln = tree.get_leaf_names()

    # odl = list()
    distl = list()
    for seq_from in ln:
        for seq_to in ln:
            if seq_from != seq_to:
                subt = tree.get_common_ancestor(seq_from, seq_to)
                fromd = subt.get_distance(seq_from)
                tod = subt.get_distance(seq_to)
                if (fromd != tree.get_distance(seq_from) and
                        tod != tree.get_distance(seq_to) and
                        subt.evoltype == 'S'):
                    distl.extend([fromd, tod])

    mean = np.mean(distl)
    med = np.median(distl)
    skew = stats.skew(distl)
    kurt = stats.kurtosis(distl)
    sd = st.stdev(distl)

    rstats = {'mean': mean, 'med': med, 'skew': skew,
              'kurt': kurt, 'sd': sd}

    return rstats


def subtree_plot(phylome_id, prot_dict, seed_id, tree):
    mphsptrees = dict()
    mphgtrees = dict()
    for i, subtree in enumerate(tree.traverse()):
        sps = [sp.split('_')[1] for sp in subtree.get_leaf_names()]
        if (len(sps) == len(set(sps)) and len(sps) > 10 and
                subtree.get_farthest_leaf()[1] != 0):
            mphsptrees[len(sps)] = subtree
        elif len(sps) > 10 and subtree.get_farthest_leaf()[1] != 0:
            mphgtrees[len(sps) - len(set(sps))] = subtree

    if len(mphsptrees) > 0:
        rtree = mphsptrees.get(max(mphsptrees.keys()))
    elif len(mphgtrees) > 0:
        rtree = mphgtrees.get(min(mphgtrees.keys()))
    else:
        rtree = tree

    root = rtree.get_common_ancestor(rtree)

    rtldist = list()
    for leaf in rtree.get_leaf_names():
        rtldist.append(rtree.get_distance(root, leaf))

    odict = {'st_%s' % seed_id: rtldist}

    return odict


def mrca_tt_plot(tree, seed_id):
    ln = tree.get_leaf_names()

    distl = list()
    for seq_from in ln:
        for seq_to in ln:
            if seq_from != seq_to:
                subt = tree.get_common_ancestor(seq_from, seq_to)
                fromd = subt.get_distance(seq_from)
                tod = subt.get_distance(seq_to)
                if (fromd != tree.get_distance(seq_from) and
                        tod != tree.get_distance(seq_to) and
                        subt.evoltype == 'S'):
                    distl.extend([fromd, tod])

    odict = {'mrca_%s' % seed_id: distl}

    return odict


def root_tt_ref():
    return 0


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, prot_dict, olist, plots):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.prot_dict = prot_dict
        self.olist = olist
        self.plots = plots

    def run(self):
        tree = self.tree_row.split('\t')
        t = ete3.PhyloTree(tree[3], sp_naming_function=get_species_tag)

        if (len(t.get_species()) > 10 and
                len(t.get_leaf_names()) < 3 * len(t.get_species())):
            print('Calculating: %s, species no.: %s, leaves no.: %s' %
                  (tree[0], len(t.get_species()), len(t.get_leaf_names())))
            root(t, root_dict[int(self.phylome_id)])
            t.get_descendant_evol_events()

            st_ref = subtree_ref(self.phylome_id, self.prot_dict, tree[0], t)
            mrca_ref = mrca_tt_ref(t, tree[0])

            tnames = t.get_leaf_names()

            self.plots.append(mrca_tt_plot(t, tree[0]))
            self.plots.append(subtree_plot(self.phylome_id, self.prot_dict,
                                           tree[0], t))

            for i, from_sp in enumerate(tnames):
                for to_sp in tnames[i + 1:]:
                    if from_sp != to_sp:
                        leaf_dist = get_dists(t, from_sp, to_sp, tree[0],
                                              self.phylome_id, self.prot_dict,
                                              st_ref, mrca_ref)
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
        cpus = 6
    else:
        ifile = options.ifile
        odir = options.odir
        prots = options.prot
        cpus = int(options.cpus)

    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]
    file_id = ifile.rsplit('/', 1)[1].split('.', 1)[0]

    dist_fn = '/'.join([odir, (file_id + '_dist.csv')])

    if not file_exists(dist_fn):
        print('Creating: ', dist_fn)

        create_folder(odir)
        create_folder('%s/dens_plots' % odir)

        trees = open(ifile, 'r').read().split('\n')
        prot_dict = csv_to_dict(prots, '\t')

        with Manager() as manager:
            olist = manager.list()
            plots = manager.list()

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
                                           prot_dict, olist, plots)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            odf = pd.DataFrame(list(olist))
            odf.to_csv(dist_fn, index=False)

            print('Plotting densities')

            mrca_pdf = backend_pdf.PdfPages('/'.join([odir,
                                                      (file_id + '_mrca.pdf')]))
            st_pdf = backend_pdf.PdfPages('/'.join([odir,
                                                    (file_id + '_st.pdf')]))
            for plotl in list(plots):
                plotdf = pd.DataFrame(plotl)
                tree_id = list(plotl.keys())[0]
                plotdf.plot.density(color='darkorange',
                                    legend=False)
                plt.title('Density plot for %s' % tree_id)
                plt.axvline(plotdf[tree_id].mean(), color='k',
                            linestyle='dashed', linewidth=1)
                plt.axvline(plotdf[tree_id].median(), color='blue',
                            linestyle='dashed', linewidth=1)

                if 'mrca_' in tree_id:
                    plt.xlabel('MRCA-to-tip distance')
                    mrca_pdf.savefig()
                elif 'st_' in tree_id:
                    plt.xlabel('Subtree MRCA to tip distance')
                    st_pdf.savefig()

            mrca_pdf.close()
            st_pdf.close()


if __name__ == '__main__':
    main()
