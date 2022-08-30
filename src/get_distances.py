#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
treefuns.py -- Functions useful to manipulate phylome trees

Requirements:
 - pandas
 - ete3
 - multiprocessing
 - treefuns

Written by Moisès Bernabeu <moigil.bernabeu@gmail.com>
August 2022
'''

import pandas as pd
import ete3
from multiprocessing import Process, Manager
from treefuns import (read_treeline, root, get_sp2age, annotate_tree,
                      get_group_species, get_group_mrca, tree_stats,
                      count_dupl_specs)
from utils import create_folder, file_exists
from optparse import OptionParser


def event_dist(tree, event_species, normgroup):
    dfd = dict()
    dfd['seed'] = tree['seed']

    tree['tree'].get_descendant_evol_events()

    # Calculating the tree stats and storing in the dictionary
    wholet_stats = tree_stats(tree['tree'])
    tdupl = count_dupl_specs(tree['tree'])
    dfd = {**dfd, **{'tree_' + k: v for k, v in wholet_stats.items()},
           **{'tree_' + k: v for k, v in tdupl.items()}}

    # Getting the normalising factor
    annotate_tree(tree['tree'], 'normalising', event_species[normgroup])
    normt = get_group_mrca(tree['tree'], tree['seed'],
                           'normalising', tree['seed'])
    normtdupl = count_dupl_specs(normt['node'])

    # Calculating the normalising subtree stats and storing in
    # the dictionary
    normt_stats = tree_stats(normt['node'])
    dfd = {**dfd, **{'norm_' + k: v for k, v in normt_stats.items()},
           **{'norm_' + k: v for k, v in normtdupl.items()}}

    dfd['wdth_ratio'] = dfd['norm_width'] / dfd['tree_width']

    # Getting the distances from the seed to the group MRCA
    for group in event_species:
        # Annotating the tree with the group info
        annotate_tree(tree['tree'], group, event_species[group])

        # Obtaining each group MRCA
        groupt = get_group_mrca(tree['tree'], tree['seed'],
                                group, tree['seed'])

        # Calculating group MRCA (event) to seed distance
        dfd[group + '_dist'] = groupt['node'].get_distance(tree['seed'])
        dfd[group + '_ndist'] = (groupt['node'].get_distance(tree['seed'])
                                 / normt_stats['median'])

    return dfd


class dist_process(Process):
    def __init__(self, tree_row, event_species, normgroup, sp2agedic, olist):
        Process.__init__(self)
        self.tree_row = tree_row
        self.event_species = event_species
        self.normgroup = normgroup
        self.sp2agedic = sp2agedic
        self.olist = olist

    def run(self):
        tree = read_treeline(self.tree_row)

        # Rooting the tree
        root(tree['tree'], self.sp2agedic)

        # Calculating distances
        odict = event_dist(tree, self.event_species, self.normgroup)
        self.olist.append(odict)


def main():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='ifile',
                      help=('File with multiple trees in newick, each line has'
                            ' the format: seed\tmodel\tlikelihood\tnewick'),
                      metavar='<file.nwk>')
    parser.add_option('-s', '--sptree', dest='sptree',
                      help='Newick file containing the species tree.',
                      metavar='<file.nwk>')
    parser.add_option('-c', '--clades', dest='cladedf',
                      help=('Species belonging to clades, columns show the '
                            'clades and rows the species belonging to them.'),
                      metavar='<file.tsv>')
    parser.add_option('-l', '--seedsp', dest='seedsp',
                      help='Species\' code for the seed.',
                      metavar='<SPECIES>')
    parser.add_option('-n', '--normgroup', dest='normgroup',
                      help='Normalisation group header in clades dataframe.',
                      metavar='<group_column_name>')
    parser.add_option('-t', '--threads', dest='threads',
                      help='Number of threads.',
                      metavar='<N>', type='int', default=4)
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='Output prefix.',
                      metavar='</path/to/dir/prefix> or <prefix>',
                      default='dist_output')
    parser.add_option('-r', '--redo', dest='redo',
                      help='Ommit done file and redo.',
                      action='store_true', default=False)
    (options, args) = parser.parse_args()

    ifile = options.ifile
    sptree = options.sptree
    cladedf = options.cladedf
    seedsp = options.seedsp
    normgroup = options.normgroup
    cpus = options.threads
    prefix = options.prefix
    redo = options.redo

    if '/' in prefix:
        odir = prefix.rsplit('/', 1)[0]
        create_folder(odir)

    ofile = '%s.csv' % prefix

    if not file_exists(ofile) or redo:
        # Computing the species to age dictionary from the species tree
        sptree = ete3.PhyloTree(sptree)
        sp2agedic = get_sp2age(sptree, seedsp)

        # Getting the speies belonging to interest groups dictionary
        event_species = get_group_species(cladedf)

        # Iterating in parallel thorugh newick file lines and appending to a
        # dictionary list
        with Manager() as manager:
            olist = manager.list()

            processes = list()
            for tree_row in open(ifile, 'r'):
                if tree_row != '':
                    if len(processes) >= cpus:
                        done = False
                        while not done:
                            for process in processes:
                                if not process.is_alive():
                                    processes.remove(process)
                                    done = True

                    process = dist_process(tree_row, event_species, normgroup,
                                           sp2agedic, olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(ofile, index=False)
    else:
        print()


if __name__ == '__main__':
    main()
