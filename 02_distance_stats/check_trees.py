#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
check_trees.py -- Check trees having bimodal distributions

Requirements:
 - ete3

Written by Moisès Bernabeu <moigil.bernabeu.sci@gmail.com>
December 2022
'''

import ete3
from treefuns import (get_species, annotate_tree, get_group_mrca,
                      get_group_species, get_fist_part_sp, get_species_sptree)
from operator import itemgetter
import pandas as pd


def analyse_tree(nwk, seed, event_species, firstsplit, sp_mono, lg):
    tree = ete3.PhyloTree(nwk, sp_naming_function=get_species)
    annotate_tree(tree, '44', event_species['44'])

    tree.get_descendant_evol_events()

    try:
        primatest = get_group_mrca(tree, seed, '44', seed)

        for leaf in primatest['node'].get_leaf_names():
            subtree = primatest['node'].get_common_ancestor(seed, leaf)
            if get_species(leaf) == sp_mono and subtree.evoltype == 'S':
                # print(seed, subtree.evoltype)
                # print(subtree)
                seqs_mono = [seed, leaf]
                odict = dict()
                odict['seed'] = seed
                odict['subtree'] = subtree
                odict['splist'] = list(subtree.get_species())
                odict['sp_in'] = sp_mono in list(subtree.get_species())
                odict['seqs_mono'] = seqs_mono
                odict['sp_to'] = sp_mono
                odict['lg'] = lg
                odict['mphy'] = len(subtree.get_common_ancestor(seqs_mono).get_species()) <= 2
            else:
                odict = {}
    except IndexError:
        print('Cannot find the Primate subtree')
        odict = {}
        print(tree)

    return odict


def main():
    event_species = get_group_species('../../04_calc_dists/data/node_data.tsv')

    event_species = get_group_species('../../04_calc_dists/data/node_data.tsv')
    sorted_keys = [(k, v, len(v)) for k, v in event_species.items()]
    sorted_keys = sorted(sorted_keys, key=itemgetter(2))
    sorted_keys = [x[0] for x in sorted_keys]

    sptree = ete3.PhyloTree(newick=open('../../04_calc_dists/data/sp_tree.nwk', 'r').read(), sp_naming_function=get_species_sptree)
    firstsplit = get_fist_part_sp(sptree, event_species,
                                  sorted_keys, 'HUMAN')

    lth05m = open('../outputs/bimlth05_macmu_seed.txt', 'r').read().split('\n')
    gth05m = open('../outputs/bimlth05_macmu_seed.txt', 'r').read().split('\n')

    lth05p = open('../outputs/bimlth05_papan_seed.txt', 'r').read().split('\n')
    gth05p = open('../outputs/bimlth05_papan_seed.txt', 'r').read().split('\n')

    odfl = list()

    # i = 0
    for treeline in open('../../04_calc_dists/data/mammal_trees_rooted.nwk', 'r'):
        if treeline != '\n':
            treeline = treeline.split('\t')
            if treeline[0] in lth05m:
                odfd = analyse_tree(treeline[1], treeline[0], event_species,
                                    firstsplit, 'MACMU', 'lower')  # Expected MACMU HUMAN monophyletic
                odfl.append(odfd)

            if treeline[0] in gth05m:
                odfd = analyse_tree(treeline[1], treeline[0], event_species,
                                    firstsplit, 'MACMU', 'greater')  # Expected MACMU HUMAN polyphyletic or absent in the st
                odfl.append(odfd)

            if treeline[0] in lth05p:
                odfd = analyse_tree(treeline[1], treeline[0], event_species,
                                    firstsplit, 'PAPAN', 'lower')  # Expected PAPAN HUMAN monophyletic
                odfl.append(odfd)

            if treeline[0] in gth05p:
                odfd = analyse_tree(treeline[1], treeline[0], event_species,
                                    firstsplit, 'PAPAN', 'greater')  # Expected PAPAN HUMAN polyphyletic or absent in the st
                odfl.append(odfd)

        # i += 1
        # if i > 1000:
        #     break

    pd.DataFrame(odfl).to_csv('../outputs/summary_trees.csv', index=False)


if __name__ == '__main__':
    main()
