#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
normalisation.py -- Normalisation functions definition script

Here are stored all the functions needed to get the normalised distances

Requirements:

Written by Miosès Bernabeu <moigil.bernabeu@gmail.com>
February 2022
'''

# Import libraries ----
import numpy as np
from scipy import stats
import statistics as st


# Definitions ----
def subtree_tt_ref(tree, tid):
    '''
    Subtree to tip distance statistics

    Subtree is selected according to maximum number of tips and only speciation
    events. Then, the subtree root-to-tip distances are retrieved. The function
    creates a list with all the distances and returns a dictionary with de
    basic descriptive statistics of its distribution.

    Args:
        tree (Phylotree): phylome tree
    # rmean = np.mean(rdist)
    # rmed = np.median(rdist)
    # rskew = stats.skew(rdist)
    # rkurt = stats.kurtosis(rdist)
    # rsd = st.stdev(rdist)
    #
    # rwidth = sum(rdist)
    Returns:
        dict: set of basic descriptive statistics of the lengths in the subtree

    Raises:
        Exception: description
    '''

    mphsptrees = dict()
    mphgtrees = dict()
    for i, subtree in enumerate(tree.traverse()):
        sps = [sp for sp in subtree.get_leaf_names()]
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

    rstats = {'id': tid, 'wdth': wdth, 'mean': mean, 'med': med,
              'skew': skew, 'kurt': kurt, 'sd': sd}

    return rstats


def mrca_tt_ref(tree, tid):
    '''
    MRCA node to tip distance statistics

    The function gets a pair of orthologous sequences and rerieves their MRCA,
    then it calculates the MRCA to each sequence distances, creates a list with
    all the distances and returns a dictionary with de basic descriptive
    statistics of its distribution.

    Args:
        tree (Phylotree): phylome tree

    Returns:
        dict: set of basic descriptive statistics of the lengths in the subtree

    Raises:
        Exception: description
    '''

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

    mean = np.mean(distl)
    med = np.median(distl)
    skew = stats.skew(distl)
    kurt = stats.kurtosis(distl)
    sd = st.stdev(distl)

    rstats = {'id': tid, 'mean': mean, 'med': med, 'skew': skew,
              'kurt': kurt, 'sd': sd}

    return rstats


def root_tt_ref(tree, tid):
    '''
    Root to tip distance statistics

    The function retrieves whole tree root-to-tip distances are retrieved. The
    function creates a list with all the distances and returns a dictionary
    with de basic descriptive statistics of its distribution.

    Args:
        tree (Phylotree): phylome tree

    Returns:
        dict: set of basic descriptive statistics of the lengths in the subtree

    Raises:
        Exception: description
    '''

    distl = list()

    for leaf in tree.get_leaf_names():
        distl.append(tree.get_distance(leaf))

    twdth = tree.get_farthest_leaf()[1]
    mean = np.mean(distl)
    med = np.median(distl)
    skew = stats.skew(distl)
    kurt = stats.kurtosis(distl)
    sd = st.stdev(distl)

    rbs = sum(distl) / twdth

    rstats = {'id': tid, 'twdth': twdth, 'mean': mean, 'med': med,
              'skew': skew,
              'kurt': kurt, 'sd': sd, 'rbs': rbs}

    return rstats


def rbls_ref(tree, tid):
    brls = [x.dist for x in tree.traverse()]
    twdth = tree.get_farthest_leaf()[1]

    rstats = {'id': tid, 'sum_brl': sum(brls), 'twdth': twdth}

    return rstats
