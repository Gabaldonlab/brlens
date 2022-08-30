#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
treefuns.py -- Functions useful to manipulate phylome trees

Requirements:
 - operator
 - scipy
 - numpy
 - ete3

Written by Moisès Bernabeu <moigil.bernabeu@gmail.com>
August 2022
'''

# Import libraries ----
from operator import itemgetter
from scipy import stats
import numpy as np
import ete3
import pandas as pd


# Define functions ----
def get_species(node):
    '''
    Get species name

    Args:
        node (TreeNode): tree node with sequence name

    Returns:
        string: the tree node species label
    '''

    if '_' in node:
        return node.split("_")[1]
    else:
        return node


def read_treeline(line):
    '''
    Function for reading a phylome tree line

    The function gets the line and splits it using the tab delimiter, then
    stores the tree seed, the inference model used, its likelihood and the
    tree as an ete3 PhyloTree object in a dictionary.

    Args:
        line (str): tree line

    Returns:
        dict: containing the seed, the model, the likelihood and the tree

    Raises:
        KeyError: there is not all the information in the line

    '''

    line = line.split('\t')
    tree_dict = dict()

    try:
        tree_dict['seed'] = line[0]
        tree_dict['model'] = line[1]
        tree_dict['likelihood'] = line[2]
        tree_dict['tree'] = ete3.PhyloTree(line[3],
                                           sp_naming_function=get_species)
    except KeyError:
        raise 'The tree line does not contain all the information.'

    return tree_dict


def get_sp2age(tree, seed):
    '''
    Species to age generator from species tree

    The function gets the species tree and the seed species, calculates the
    distances between the seed and the other species, then they are sorted
    using the distance and an age number is assigned, if two species are
    equally distant, the assigned age number is the same.

    Args:
        tree (PhyloTree): species tree
        seed (str):

    Returns:
        dict: where keys are the species and values the
    '''

    # Creating an auxiliary dictionary
    distdict = dict()

    # Iterating through the leaves and storing the distance to the seed
    for tip in tree.get_leaf_names():
        if tip != seed:
            dist = tree.get_distance(tip, seed)
            distdict[tip] = dist

    # Sorting the distance dictionary
    distdict = {k: v for k, v in sorted(distdict.items(),
                                        key=lambda item: item[1])}

    # Creating the species to age output dictionary
    sp2agedic = dict()
    i = 1
    previtem = list(distdict.keys())[0]

    for item in distdict:
        # Checking whether the item has the same distance than previous one
        if distdict[item] == distdict[previtem]:
            i = i
        else:
            i += 1

        sp2agedic[item] = i
        previtem = item

    return sp2agedic


def get_group_species(tab_file):
    '''
    Converts a table with species of the group in lines to a dictionary

    The function gets a table which contains the species names in rows that
    belong to a group (each column). Then a dictionary is computed where the
    key is the group name, and the value is a list of the species.

    Args:
        tab_file (str): directory of the tsv file with the group information

    Returns:
        dict: dictionary with the group as key and the species list as value
    '''

    # Group
    group_species = dict()

    # Reading the table
    df = pd.read_csv(tab_file, sep='\t')

    # Iterating through the dataframe columns
    for item in df:
        # Removing nans if there are
        spl = [x for x in list(df[item]) if str(x) != 'nan']
        group_species[item] = spl

    return group_species


def root(tree, root_dict):
    '''
    Root the tree according to a rooting dictionary

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
        ogsps = [k for k, val in root_dict.items()
                 if val == ogdval and k in tree.get_species()][0]
        ogseq = [seq for seq in tree.get_leaf_names() if ogsps in seq][0]
    else:
        ogseq = tree.get_farthest_leaf()[0].get_leaf_names()[0]

    tree.set_outgroup(ogseq)

    return ogseq


def annotate_tree(tree, key, splist):
    '''
    Tree leaves annotation

    Tree leaves are annotated with the columns

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

        df (DataFrame): pandas dataframe containing the information to
        annotate the tree

        spcol (string): the column containing the species names

        cols (string or list of strings): column or columns containing the
        annotations

    Returns:
        string: the input tree is annotated, the funtion returns 0

    Raises:
        Exception: description
    '''

    # Iterate the leaves
    for leaf in tree.get_leaves():
        # Iterate the dataframe columns
        sp = list(leaf.get_species())[0]
        if sp in splist:
            leaf.add_feature(key, key)
        else:
            leaf.add_feature(key, 'nan')

    return 0


def tree_stats(tree):
    '''
    Get tree branch stats

    From a tree the function retrieves basic numerical information about the
    branches lengths. First, it creates a list where stores all the root to
    tip distances, then calculates the median, the mean, the width and the sum

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

    Returns:
        dictionary: dictionary with the tree statistics

    Raises:
        Exception: description
    '''

    # Getting all the leaves
    ndlf = tree.get_leaves()

    # Retrieving the root to tip distances
    distl = list()
    for leaf in ndlf:
        distl.append(tree.get_distance(leaf))

    # Generating the output dictionary
    nodedict = dict()
    nodedict['leafno'] = len(distl)
    nodedict['median'] = np.median(distl)
    nodedict['mean'] = np.mean(distl)
    nodedict['width'] = tree.get_farthest_leaf()[1]
    nodedict['sum'] = sum(distl)
    nodedict['kurt'] = stats.kurtosis(distl)
    nodedict['skew'] = stats.skew(distl)

    return nodedict


def get_group_mrca(tree, tree_id, feature, value, sp_in=None):
    '''
    Get the greatest monophyletic subtree of a labeled group

    The function retrieves a subtree which all leaves have the value of value
    of the feature indicated. Check all groups with these conditions and
    returns one that maximizes the number of leaves.

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

        tree_id (string): the phylome tree idea

        feature (string): the name of the feature where the MRCA label
        is stored

        value (string): the name of the value of the feature that defines
        the clade

        sp_in (string): the name of the species that has to be inside the
        MRCA group

    Returns:
        dictionary: dictionary with some information about the node and
        the node

    Raises:
        Exception: description
    '''

    # Get the number of leaves in the tree
    tlno = len(tree.get_leaf_names())

    # Get the monophyletic sets containing all species with col attribute
    mphysets = list()
    for st in tree.traverse():
        # Get the tree leaves
        leaves = st.get_leaves()
        stlno = len(leaves)

        # Width of the monophyletic group
        stwdth = st.get_farthest_leaf()[1]

        # Generate list of col attributes in subtree
        feat_list = list()
        for stleag in st.get_leaves():
            feat_list.append(getattr(stleag, feature))

        # If the species to include is None we get the first leaf
        lnames = st.get_leaf_names()
        if sp_in is None:
            sptoincl = lnames[0]
        else:
            sptoincl = sp_in

        # Check the subtree only contains 1 attributes, it has more
        # than 1 leaves, it is not the entire tree and it is not a polytomy
        if (sptoincl in lnames and len(set(feat_list)) == 1 and stlno > 1 and
                stlno != tlno and stwdth != 0):
            # Appending to a list a dictionary with the basic information of
            # the group monophyletic group
            mphy = dict()
            mphy['tree'] = tree_id
            mphy['node'] = st
            mphy['seq_no'] = stlno
            mphy[feature] = feat_list[0]
            mphysets.append((stlno, mphy))

    # Get the maximum leaves monophyletic group for each group
    mphylist = sorted(mphysets, key=itemgetter(0), reverse=True)

    return mphylist[0][1]


def count_dupl_specs(tree):
    cdict = {'D': 0, 'S': 0}
    for node in tree.iter_search_nodes():
        if len(node.get_leaf_names()) > 1:
            cdict[node.evoltype] += 1

    return cdict
