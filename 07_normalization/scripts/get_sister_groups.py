#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
get_sister_groups.py -- Retrieve sister groups from a species tree



Requirements: ete3

Written by Miosès Bernabeu <moigil.bernabeu@gmail.com>
February 2022
'''

import ete3
import pandas as pd


def get_sp_sisters(t, ref_sp):
    from_node = t.get_leaves_by_name(ref_sp)[0]
    sist_group = dict()

    n = 0

    fst_sister = from_node.get_sisters()
    if len(fst_sister) > 0:
        sln = fst_sister[0].get_leaf_names()
        for leaf in sln:
            sist_group[leaf] = n
    n += 1

    for i in from_node.get_ancestors():
        sister = i.get_sisters()
        if len(sister) > 0:
            sln = sister[0].get_leaf_names()
            for leaf in sln:
                sist_group[leaf] = n
        n += 1
    return sist_group


def main():
    t = ete3.PhyloTree('../data/0005_sptree.nwk')
    a = get_sp_sisters(t, 'YEAST')
    pd.DataFrame([a]).transpose().to_csv('../../09_stats/data/0005_sister_group.csv')

    t = ete3.PhyloTree('../data/0076_sptree.nwk')
    a = get_sp_sisters(t, 'HUMAN')
    pd.DataFrame([a]).transpose().to_csv('../../09_stats/data/0076_sister_group.csv')


if __name__ == '__main__':
    main()
