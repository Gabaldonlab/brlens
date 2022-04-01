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
import numpy as np

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
    if any(sp in root_dict for sp in tree.get_species()):
        ogdval = max([root_dict.get(sp, 0) for sp in tree.get_species()])
        ogsps = [k for k, val in root_dict.items() if val == ogdval and k in tree.get_species()][0]
        ogseq = [seq for seq in tree.get_leaf_names() if ogsps in seq][0]
    else:
        ogseq = tree.get_farthest_leaf()[0].get_leaf_names()[0]

    tree.set_outgroup(ogseq)

    return ogseq


def annotate_lineages(tree, gnmdf, cols):
    for leaf in tree.get_leaves():
        for col in cols:
            sp = list(leaf.get_species())[0]
            feat_val = list(gnmdf[col][gnmdf['Proteome'] == sp])[0]
            leaf.add_feature(col, feat_val)

    return 0


def norm_factor(nodedict):
    ndlf = nodedict['node'].get_leaves()
    distl = list()
    for leaf in ndlf:
        distl.append(nodedict['node'].get_distance(leaf))
    nodedict['norm_factor'] = np.median(distl)

    return nodedict


def clade_norm(tree, tree_id, col):
    tlno = len(tree.get_leaf_names())

    mphysets = list()
    for st in tree.traverse():
        leaves = st.get_leaves()

        stwdth = st.get_farthest_leaf()[1]

        feat_list = list()
        for stleag in st.get_leaves():
            feat_list.append(getattr(stleag, col))

        if (len(set(feat_list)) == 1 and len(leaves) > 1 and
                len(leaves) != tlno and stwdth != 0):
            mphy = dict()
            mphy['tree'] = tree_id
            mphy['node'] = st
            mphy['seq_no'] = len(leaves)
            mphy[col] = feat_list[0]
            mphysets.append(mphy)

    nodedf = pd.DataFrame(mphysets)
    indexes = nodedf.index

    mphylist = list()
    for group in set(nodedf[col]):
        if str(group) != 'nan':
            maxval = max(nodedf.loc[nodedf[col] == group]['seq_no'])
            dfindex = indexes[(nodedf[col] == group) &
                              (nodedf['seq_no'] == maxval)]
            mphylist.append(mphysets[dfindex.values[0]])

    nodel = list()
    for nodedict in mphylist:
        nodel.append(norm_factor(nodedict))

    return nodel[0]


def get_group_mrca(tree, tree_id, col, value, seed):
    mphysets = list()
    tlno = len(tree.get_leaf_names())

    for st in tree.traverse():
        leaves = st.get_leaves()
        if len(leaves) > 1:
            stwdth = st.get_farthest_leaf()[1]

            feat_list = list()
            for stleag in st.get_leaves():
                feat_list.append(getattr(stleag, col))

            if (seed in st.get_leaf_names() and
                    len(set(feat_list)) == 1 and len(leaves) > 1 and
                    len(leaves) != tlno and stwdth != 0):
                mphy = dict()
                mphy['tree'] = tree_id
                mphy['node'] = st
                mphy['seq_no'] = len(leaves)
                mphy[col] = feat_list[0]
                mphysets.append(mphy)

    nodedf = pd.DataFrame(mphysets)
    indexes = nodedf.index

    mphylist = list()
    for group in set(nodedf[col]):
        if str(group) != 'nan':
            maxval = max(nodedf.loc[nodedf[col] == group]['seq_no'])
            dfindex = indexes[(nodedf[col] == group) &
                              (nodedf['seq_no'] == maxval)]
            mphylist.append(mphysets[dfindex.values[0]])

    return mphylist[0]


def main():
    # Script options definition ----
    parser = OptionParser()
    parser.add_option('-d', '--def', dest='default',
                      help='Default configurations to execute with our data.',
                      action='store_true')
    parser.add_option('-f', '--file', dest='ifile',
                      help='In file',
                      metavar='<path/to/file.txt>')
    parser.add_option('-g', '--groups', dest='groups',
                      help='Groups file (.csv)',
                      metavar='<path/to/file.csv>')
    parser.add_option('-o', '--out', dest='output',
                      help='output directory',
                      metavar='<path/to/folder>')
    (options, args) = parser.parse_args()

    if options.default:
        infile = '../data/0076_108.txt'
        gnmdffile = '../data/0076_norm_groups.csv'
        outdir = '../outputs/'
    else:
        infile = options.ifile
        gnmdffile = options.groups
        outdir = options.output

    ofilenm = infile.rsplit('/', 1)[1].split('.', 1)[0]
    ofile = '%s/%s_dist.csv' % (outdir, ofilenm)

    if not file_exists(ofile):
        create_folder(outdir)

        gnmdf = pd.read_csv(gnmdffile)
        phylome_id = infile.rsplit('/', 1)[1].split('_', 1)[0]

        olist = list()
        for tree in open(infile, 'r'):
            treel = tree.split('\t')
            t = ete3.PhyloTree(treel[3], sp_naming_function=get_species_tag)

            root(t, root_dict[int(phylome_id)])
            t.get_descendant_evol_events()

            annotate_lineages(t, gnmdf, ['Normalising group',
                                         'Vertebrate',
                                         'Metazoan'])

            norm_dict = clade_norm(t, treel[0], 'Normalising group')

            nfactor = norm_dict['norm_factor']

            vert_dict = get_group_mrca(t, treel[0], 'Vertebrate',
                                       'vertebrate', treel[0])
            met_dict = get_group_mrca(t, treel[0], 'Metazoan',
                                      'metazoan', treel[0])

            odict = dict()
            odict['seed'] = treel[0]
            odict['species'] = get_species_tag(treel[0])
            odict['vert_dist'] = vert_dict['node'].get_distance(treel[0])
            odict['met_dist'] = met_dict['node'].get_distance(treel[0])
            odict['seed_dist'] = t.get_distance(treel[0])
            odict['vert_ndist'] = vert_dict['node'].get_distance(treel[0]) / nfactor
            odict['met_ndist'] = met_dict['node'].get_distance(treel[0]) / nfactor
            odict['seed_ndist'] = t.get_distance(treel[0]) / nfactor
            olist.append(odict)

        odf = pd.DataFrame(olist)
        odf.to_csv(ofile, index=False)

    return 0


if __name__ == '__main__':
    main()
