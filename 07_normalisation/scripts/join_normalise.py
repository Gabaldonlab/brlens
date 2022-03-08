#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
join_normalise.py -- Phylome partitions join and calculate normalised distances



Requirements: pandas

Written by Moisès Bernabeu <moigil.bernabeu@gmail.com>
March 2022
'''

# Import libraries ----
import pandas as pd
from glob import glob

# Script
def main():
    files = glob('../outputs/*_*_norm.csv')
    print(files)

    ids = list()
    for file in files:
        ids.append(file.rsplit('/', 1)[1].split('_', 1)[0])
    ids = set(ids)

    for phyid in ids:
        print('Parsing: ', phyid)
        normfiles = glob('../outputs/%s*_norm.csv' % phyid)
        distfiles = glob('../outputs/%s*_dist.csv' % phyid)

        for i, file in enumerate(normfiles):
            if i == 0:
                normdf = pd.read_csv(file)
            else:
                normdf = pd.concat([normdf, pd.read_csv(file)])

        for i, file in enumerate(distfiles):
            if i == 0:
                distdf = pd.read_csv(file)
            else:
                distdf = pd.concat([distdf, pd.read_csv(file)])

        tree_meds = normdf[(normdf['subtree'] == 'whole') & normdf['method'] == 'rbl')].groupby('tree').median()
        rbl_nf = tree_meds['med'] / tree_meds.median()['med']

        tree_meds = normdf[(normdf['subtree'] == 'whole') & normdf['method'] == 'mrca')].groupby('tree').median()
        mrca_nf = tree_meds['med'] / tree_meds.median()['med']

        distdf['rbl_ndist'] = distdf.dist / distdf.tree.map(rbl_nf)
        distdf['mrca_ndist'] = distdf.dist / distdf.tree.map(mrca_nf)

        normdf.to_csv('../outputs/%s_norm.csv' % phyid, index=False)
        distdf.to_csv('../outputs/%s_dist.csv' % phyid, index=False)


if __name__ == '__main__':
    main()
