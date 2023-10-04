#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
get_t2in_dist.py -- Tip-to-internode distance calculation

This script calculates the distances from the seed species to the
nodes whose MRCA contains the species in the different columns
of the table.

Copyright (C) 2022  Moisès Bernabeu <moigil.bernabeu.sci@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import pandas as pd
import ete3
from multiprocessing import Process, Manager
from treefuns import (read_treeline, root, get_sp2age, annotate_tree,
                      get_group_species, get_group_mrca,
                      tree_stats, count_dupl_specs, get_species_sptree)
from utils import create_folder, file_exists
from optparse import OptionParser
from operator import itemgetter


def get_fist_part_sp(sptree, event_species, sorted_keys, seed):
    # Getting the evolutionary events of the tree (MRCA function requirement)
    sptree.get_descendant_evol_events()

    # Generating the output dictionary
    odict = dict()

    # Iterating thorugh the groups in the event_species dictionary
    for group in sorted_keys:
        # Annotating the tree with the node
        annotate_tree(sptree, group, event_species[group])

        # Getting the MRCA group
        spgroupt = get_group_mrca(sptree, seed, group, seed)

        # Writing in a tupple the species after the first partition
        # in two groups
        i = 0
        ol = list()
        for node in spgroupt['node'].iter_descendants():
            if i < 2:
                ol.append(node.get_species())
                i += 1
            else:
                break

        odict[group] = (ol[0], ol[1])

    return odict


def event_dist(tree, event_species, sorted_keys, normgroup, firstsplit):
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

    try:
        normt = get_group_mrca(tree['tree'], tree['seed'],
                               'normalising', tree['seed'],
                               firstsplit[normgroup])
        normtdupl = count_dupl_specs(normt['node'])
    except IndexError:
        print('Tree %s: cannot compute the normalising group.' % tree['seed'])
        normt = None

    if normt is not None:
        # Calculating the normalising subtree stats and storing in
        # the dictionary
        normt_stats = tree_stats(normt['node'])
        dfd = {**dfd, **{'norm_' + k: v for k, v in normt_stats.items()},
               **{'norm_' + k: v for k, v in normtdupl.items()}}

        dfd['wdth_ratio'] = dfd['norm_width'] / dfd['tree_width']

        # Getting the distances from the seed to the group MRCA
        last_group = [0, 0]
        for group in sorted_keys:
            # Annotating the tree with the group info
            annotate_tree(tree['tree'], group, event_species[group])

            # Obtaining each group MRCA
            try:
                groupt = get_group_mrca(tree['tree'], tree['seed'],
                                        group, tree['seed'],
                                        firstsplit[group])
            except IndexError:
                print('Distance to node %s in tree %s cannot be computed.' %
                      (group, tree['seed']))
                groupt = None

            if groupt is not None:
                # Checking whether the subtrees are the same
                if (last_group[0] != groupt['seq_no'] and
                        last_group[1] != groupt['sp_no']):
                    # Calculating group MRCA (event) to seed distance
                    dist = groupt['node'].get_distance(tree['seed'])
                    ndist = dist / normt_stats['median_r2t']
                    dfd[group + '_dist'] = dist
                    dfd[group + '_ndist'] = ndist

                last_group = [groupt['seq_no'], groupt['sp_no']]

    return dfd


class dist_process(Process):
    def __init__(self, tree_row, event_species, sorted_keys,
                 normgroup, sp2agedic, olist, firstsplit):
        Process.__init__(self)
        self.tree_row = tree_row
        self.event_species = event_species
        self.sorted_keys = sorted_keys
        self.normgroup = normgroup
        self.sp2agedic = sp2agedic
        self.olist = olist
        self.firstsplit = firstsplit

    def run(self):
        tree = read_treeline(self.tree_row)

        # Rooting the tree
        root(tree['tree'], self.sp2agedic)

        # Calculating distances
        odict = event_dist(tree, self.event_species,
                           self.sorted_keys, self.normgroup, self.firstsplit)
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
        sptree = ete3.PhyloTree(sptree, sp_naming_function=get_species_sptree)
        sp2agedic = get_sp2age(sptree, seedsp)

        # Getting the speies belonging to interest groups dictionary
        event_species = get_group_species(cladedf)
        sorted_keys = [(k, v, len(v)) for k, v in event_species.items()]
        sorted_keys = sorted(sorted_keys, key=itemgetter(2))
        sorted_keys = [x[0] for x in sorted_keys]

        # Getting the species tree first split species for each group
        firstsplit = get_fist_part_sp(sptree, event_species,
                                      sorted_keys, seedsp)

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

                    process = dist_process(tree_row, event_species, sorted_keys,
                                           normgroup, sp2agedic, olist,
                                           firstsplit)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(ofile, index=False)


if __name__ == '__main__':
    main()
