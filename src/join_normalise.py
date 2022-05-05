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
    files = glob('../inphy_outputs/*_*_norm.csv')
    print(files)

    ids = list()
    for file in files:
        ids.append(file.rsplit('/', 1)[1].split('_', 1)[0])
    ids = set(ids)

    for phyid in ids:
        print('Parsing: ', phyid)
        normfiles = glob('../inphy_outputs/%s*_norm.csv' % phyid)
        distfiles = glob('../inphy_outputs/%s*_dist.csv' % phyid)

        for i, file in enumerate(distfiles):
            if i == 0:
                distdf = pd.read_csv(file)
            else:
                distdf = pd.concat([distdf, pd.read_csv(file)])

        distdf.to_csv('../inphy_outputs/%s_dist.csv' % phyid, index=False)

        for i, file in enumerate(normfiles):
            if i == 0:
                normdf = pd.read_csv(file)
            else:
                normdf = pd.concat([normdf, pd.read_csv(file)])

        normdf.to_csv('../inphy_outputs/%s_norm.csv' % phyid, index=False)


if __name__ == '__main__':
    main()
