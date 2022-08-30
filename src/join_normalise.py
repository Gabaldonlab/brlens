#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
join_normalise.py -- Phylome partitions join and calculate normalised distances



Requirements: pandas

Written by Moisès Bernabeu <moigil.bernabeu.sci@gmail.com>
March 2022
'''

# Import libraries ----
import pandas as pd
from glob import glob
from optparse import OptionParser


# Script
def main():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='idir',
                      help=('Directory with multiple get_distances outputs.'),
                      metavar='<file.csv>')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='Output prefix.',
                      metavar='</path/to/dir/prefix> or <prefix>',
                      default='dist_output')
    (options, args) = parser.parse_args()

    files = glob('%s/*.csv' % options.idir)
    print(files)

    i = 0
    for file in files:
        print('Parsing: ', file)

        if i == 0:
            distdf = pd.read_csv(file)
        else:
            distdf = pd.concat([distdf, pd.read_csv(file)])

        i += 1

    distdf.to_csv('%s/%s.csv' % (options.idir, options.prefix), index=False)


if __name__ == '__main__':
    main()
