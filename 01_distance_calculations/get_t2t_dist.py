#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
get_t2t_dist.py -- Tip-to-tip distances and normalisation

The script gets the sequence to sequence distances and calculates the
normalisation factors. It does not calculate the normalised distance.

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

from multiprocessing import Process, Manager
from optparse import OptionParser
from treefuns import (read_treeline, get_sp2age, annotate_tree,
                      get_group_mrca, tree_stats, get_species,
                      count_dupl_specs, get_group_species)
import pandas as pd
from utils import create_folder, file_exists


# Definitions ----
class dist_process(Process):
    def __init__(self, line, event_species, normgroup, olist):
        Process.__init__(self)
        self.line = line
        self.event_species = event_species
        self.normgroup = normgroup
        self.olist = olist
    

    def get_dists(self, tree, from_seq, to_seq, normfact):
        # Calculating the distances
        dist = tree['tree'].get_distance(from_seq, to_seq)
        ndist = dist / normfact

        # Getting the duplication and speciation events in the path between tips
        st = tree['tree'].get_common_ancestor(from_seq, to_seq)
        stcount = count_dupl_specs(st)

        # Creating and adding data to the output dictionary
        leafdistd = dict()
        leafdistd['tree'] = tree['seed']
        leafdistd['from'] = from_seq
        leafdistd['from_seq'] = get_species(from_seq)
        leafdistd['to'] = to_seq
        leafdistd['to_seq'] = get_species(to_seq)
        leafdistd['MRCA_type'] = st.evoltype
        leafdistd['sp_count'] = stcount['S']
        leafdistd['dup_count'] = stcount['D']
        leafdistd['dist'] = dist
        leafdistd['ndist'] = ndist

        self.olist.append(leafdistd)


    def run(self):
        tree = read_treeline(self.line)
        if (len(tree['tree'].get_species()) > 10 and
            len(tree['tree'].get_leaf_names()) < 3 * len(tree['tree'].get_species())):
            print('Processing:', tree['seed'])

            tree['tree'].set_outgroup(tree['tree'].get_midpoint_outgroup())
            tree['tree'].get_descendant_evol_events()

            annotate_tree(tree['tree'], 'normalising', self.event_species[self.normgroup])

            try:
                normt = get_group_mrca(tree['tree'], tree['seed'],
                                'normalising', tree['seed'])
                normt_stats = tree_stats(normt['node'])
            except IndexError:
                print('Tree %s: cannot compute the normalising group.' % tree['seed'])
                normt_stats = None

            if normt_stats is not None:
                tnames = tree['tree'].get_leaf_names()
                for i, from_seq in enumerate(tnames):
                    for to_seq in tnames[i + 1:]:
                        self.get_dists(tree, from_seq, to_seq, normt_stats['median_r2t'])


def main():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='ifile',
                      help=('File with multiple trees in newick, each line has'
                            ' the format: seed\tmodel\tlikelihood\tnewick'),
                      metavar='<file.nwk>')
    parser.add_option('-c', '--clades', dest='cladedf',
                      help=('Species belonging to clades, columns show the '
                            'clades and rows the species belonging to them.'),
                      metavar='<file.tsv>')
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
    cladedf = options.cladedf
    normgroup = options.normgroup
    cpus = options.threads
    prefix = options.prefix
    redo = options.redo
    
    event_species = get_group_species(cladedf)

    if '/' in prefix:
        odir = prefix.rsplit('/', 1)[0]
        create_folder(odir)

    ofile = '%s.csv' % prefix

    if not file_exists(ofile) or redo:
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

                    process = dist_process(tree_row, event_species, normgroup, olist)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(ofile, index=False)


if __name__ == '__main__':
    main()
