#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
uniprot_parse.py -- Title of the script

Written by Moisès Bernabeu <mail@mail.com>
January 2022
'''

# Import libraries ----
import sys
import os

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import *


# Definitions ----
def parse_up(file):
    protinfo = dict()
    for line in open(file):
        if 'DE   ' in line and 'RecName:' in line:
            # Gene Name
            prot = line[14:].split('=', 1)[1].replace(';', '')
            prot = prot.replace('\n', '').rsplit(' ', 1)[0]
            protinfo['prot'] = prot
        elif 'GN   ' in line:
            # Gene abbreviation
            for item in line[5:].split(';'):
                if 'Name=' in item:
                    genenm = item.replace('Name=', '').rsplit(' ', 1)[0]
                    protinfo['gene'] = genenm
        elif 'CC   ' in line:
            if 'FUNCTION: ' in line:
                toapp = 'fun'
                line = line[9:].replace('FUNCTION: ', '')
                protinfo[toapp] = line.replace('\n', '')
            elif 'LOCATION' in line:
                toapp = 'loc'
                line = line[9:].replace('SUBCELLULAR LOCATION: ', '')
                protinfo[toapp] = line.replace('\n', '')
            elif ('-!-' in line or 'Copyrighted' in line
                  or 'Distributed' in line or '--------' in line):
                toapp = 'Others'
                protinfo[toapp] = ''
            else:
                protinfo[toapp] += line[9:].replace('\n', '')

    protinfo = rmkey(protinfo, 'Others')

    locs = protinfo.get('loc')
    locs = locs.split(' Note=', 1)[0].split('. ')

    oline = file.rsplit('/', 1)[1].rsplit('.', 1)[0] + '\t'
    for i, item in enumerate(protinfo):
        if '{' in protinfo[item] and '}' in protinfo[item]:
            oline += rmbwbr(protinfo[item], '{', '}') + '\t'
        else:
            oline += protinfo[item] + '\t'

    for i, loc in enumerate(locs):
        sloc = loc.split(',', 1)[0].split(' ', 1)[0].split('{', 1)[0]
        if i != len(locs):
            oline += sloc + '\t'
        else:
            oline += sloc

    return oline
