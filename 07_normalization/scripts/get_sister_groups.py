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

    sist_list = list()

    n = 0
    m = 1

    fst_sister = from_node.get_sisters()
    if len(fst_sister) > 0:
        sln = fst_sister[0].get_leaf_names()
        for leaf in sln:
            if len(sln) == 1:
                m -= 1
            sist_group = dict()
            sist_group['leaf'] = leaf
            sist_group['group'] = n
            sist_group['mphy_sister'] = m
            sist_list.append(sist_group)

    n += 1
    m += 1

    lineage = from_node.get_ancestors()

    for i in lineage:
        sister = i.get_sisters()
        if len(sister) > 0:
            sln = sister[0].get_leaf_names()
            if len(sln) == 1:
                m -= 1
            for leaf in sln:
                sist_group = dict()
                sist_group['leaf'] = leaf
                sist_group['group'] = n
                if lineage.index(i) != len(lineage) - 2:
                    sist_group['mphy_sister'] = m
                else:
                    sist_group['mphy_sister'] = 'og'
                sist_list.append(sist_group)

        n += 1
        m += 1

    return sist_list


def main():
    t = ete3.PhyloTree('../data/0005_sptree.nwk')
    print(t)
    a = get_sp_sisters(t, 'YEAST')
    pd.DataFrame(a).to_csv('../outputs/0005_sister_group.csv',
                           index=False)

    t = ete3.PhyloTree('../data/0076_sptree.nwk')
    print(t)
    a = get_sp_sisters(t, 'HUMAN')
    pd.DataFrame(a).to_csv('../outputs/0076_sister_group.csv',
                           index=False)


if __name__ == '__main__':
    main()
