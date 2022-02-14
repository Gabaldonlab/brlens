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
from glob import glob
import gzip

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *

# Script options definition ----
parser = OptionParser()
parser.add_option('-d', '--def', dest='default',
                  help='Default configurations to execute with our data.',
                  action='store_true')
parser.add_option('-f', '--fold', dest='folder',
                  help='In folder containing proteomes',
                  metavar='<path/to/folder>')
parser.add_option('-o', '--out', dest='ofolder',
                  help='Out folder',
                  metavar='<path/to/folder>')
(options, args) = parser.parse_args()


# Definitions ----
def split_proteome(proteome, ofolder):
    ostr = ''
    protid = 0
    protnm = proteome.rsplit('/', 1)[1].split('.', 1)[0]
    for line in gzip.open(proteome, 'rt'):
        if '>' in line:
            if ostr.count('>') == 500:
                print(protid)
                ofile = gzip.open('%s/%s_%s.txt.gz' %
                                  (ofolder, protnm, str(protid).zfill(3)),
                                  'wt')
                ofile.write(ostr)
                ofile.close()
                ostr = line
                protid += 1
            else:
                ostr += line
        else:
            ostr += line
    print(protid)
    ofile = gzip.open('%s/%s_%s.txt.gz' %
                      (ofolder, protnm, str(protid).zfill(3)), 'wt')
    ofile.write(ostr)
    ofile.close()


def main():
    if options.default:
        files = glob('../data/proteomes/*')
        ofolder = '../data/splitted'
    else:
        files = glob('%s/*' % options.folder)
        ofolder = options.ofolder

    create_folder(ofolder)

    for file in files:
        print(file)
        split_proteome(file, ofolder)


if __name__ == '__main__':
    main()
