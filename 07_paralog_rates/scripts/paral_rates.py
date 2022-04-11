#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
00_script.py -- Title of the script

Brief description

Requirements:

Written by Name <mail@mail.com>
April 2022
'''

# Import libraries ----
import sys
import os
from optparse import OptionParser
from rooted_phylomes import ROOTED_PHYLOMES as root_dict
import ete3
import pandas as pd
from multiprocessing import Process, Manager
from treefuns import get_species, root, annotate_tree, \
    tree_stats, get_group_mrca, count_dupl_specs
from operator import itemgetter

# Path configuration to import utils ----
filedir = os.path.abspath(__file__)
projdir = filedir.rsplit('/', 3)[0]
sys.path.append(projdir)

from utils import file_exists, create_folder


# Definitions ----
def get_paralstats(tree):
    parlist = list()

    # Analyse all the subtrees from the tree
    for st in tree.traverse():
        # Get those subtrees that root in a duplication node
        if (len(parlist) == 0 and len(st.get_leaf_names()) > 1 and
                st.evoltype == 'D'):
            childs = st.get_children()

            # Check both child clades have more than 4 and there are 2 childs
            if len(childs) == 2 and all([len(ch.get_leaf_names()) > 4
                                         for ch in childs]):
                chevtypes = list()
                # Go through childs storing the descendant node evoltytpe
                for child in childs:
                    chevtypes += [nd.evoltype for nd in child.iter_search_nodes()
                                  if 'evoltype' in list(nd.features)]

                # Check all the descendant nodes from the childs are speciation
                # nodes and then store them in a list of tupples
                if set(chevtypes) == {'S'}:
                    for child in st.get_children():
                        chstats = tree_stats(child)
                        chstats = {**chstats, **count_dupl_specs(child)}
                        parlist.append((chstats['median'], chstats))
                        print('added to parlist')

    parlist = sorted(parlist, key=itemgetter(0))

    pardict = {**{'min_' + k: v for k, v in parlist[0][1].items()},
               **{'max_' + k: v for k, v in parlist[1][1].items()}}

    return pardict


def get_ndists(tree, phylome_id, gnmdf):
    treel = tree.split('\t')
    print('Calculating: ', treel[0])
    t = ete3.PhyloTree(treel[3], sp_naming_function=get_species)

    root(t, root_dict[int(phylome_id)])
    t.get_descendant_evol_events()

    annotate_tree(t, gnmdf, 'Proteome', ['Normalising group',
                                         'Vertebrate',
                                         'Metazoan'])

    whole_stats = tree_stats(t)

    norm_group = get_group_mrca(t, treel[0], 'Normalising group', 'A')
    norm_stats = tree_stats(norm_group['node'])

    pardict = get_paralstats(t)

    odict = {**{'whole_' + k: v for k, v in whole_stats.items()},
             **{'whole_' + k: v for k, v in count_dupl_specs(t).items()},
             **{'norm_' + k: v for k, v in norm_stats.items()},
             **{'norm_' + k: v for k, v in count_dupl_specs(norm_group['node']).items()},
             **pardict}

    return odict


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, gnmdf, olist):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.gnmdf = gnmdf
        self.olist = olist

    def run(self):
        try:
            odict = get_ndists(self.tree_row, self.phylome_id, self.gnmdf)
            self.olist.append(odict)
        except IndexError:
            print('Not only speciation paralogs found')


def main():
    # Script options definition ----
    parser = OptionParser()
    parser.add_option('-d', '--def', dest='default',
                      help='Default configurations to execute with our data.',
                      action='store_true')
    parser.add_option('-f', '--file', dest='ifile',
                      help='In file',
                      metavar='<path/to/file.txt>')
    parser.add_option('-g', '--groups', dest='groups',
                      help='Groups file (.csv)',
                      metavar='<path/to/file.csv>')
    parser.add_option('-o', '--out', dest='output',
                      help='output directory',
                      metavar='<path/to/folder>')
    parser.add_option('-c', '--cpu', dest='cpus',
                      help='Number of CPUs', type='int',
                      metavar='<N>')
    (options, args) = parser.parse_args()

    if options.default:
        infile = '../data/0076_108.txt'
        gnmdffile = '../data/0076_norm_groups.csv'
        outdir = '../outputs/'
        cpus = 1
    else:
        infile = options.ifile
        gnmdffile = options.groups
        outdir = options.output
        cpus = options.cpus

    ofilenm = infile.rsplit('/', 1)[1].split('.', 1)[0]
    ofile = '%s/%s_dist.csv' % (outdir, ofilenm)

    if not file_exists(ofile):
        create_folder(outdir)

        # # Lineal code ----
        # gnmdf = pd.read_csv(gnmdffile)
        # phylome_id = infile.rsplit('/', 1)[1].split('_', 1)[0]
        #
        # olist = list()
        # for tree_row in open(infile, 'r'):
        #     if tree_row != '':
        #         try:
        #             odict = get_ndists(tree_row, phylome_id, gnmdf)
        #             olist.append(odict)
        #             print('Written')
        #         except IndexError:
        #             print('Not only speciation paralogs found')
        #
        # # Writing output files
        # odf = pd.DataFrame(list(olist))
        # odf.to_csv(ofile, index=False)

        # Parallel code ----
        gnmdf = pd.read_csv(gnmdffile)
        phylome_id = infile.rsplit('/', 1)[1].split('_', 1)[0]

        with Manager() as manager:
            olist = manager.list()

            processes = list()
            for tree_row in open(infile, 'r'):
                if tree_row != '':
                    if len(processes) >= cpus:
                        done = False
                        while not done:
                            for process in processes:
                                if not process.is_alive():
                                    processes.remove(process)
                                    done = True

                    process = dist_process(tree_row, phylome_id, gnmdf, olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(ofile, index=False)

    return 0


if __name__ == '__main__':
    main()
