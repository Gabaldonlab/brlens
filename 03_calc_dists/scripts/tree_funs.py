#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
tree_funs.py -- Define the functions to manage trees

Here are defined the functions to root, calculate interest distances and
retrieve data from them.

Requirements: ete3

Written by Name <mail@mail.com>
Month 2022
'''

# Import libraries ----
import ete3
import pandas as pd
from rooted_phylomes import ROOTED_PHYLOMES


# Definitions ----
def get_species_tag(node):
    return node.split("_")[1]


# Read tree
file_path = '03_calc_dists/data/0003_best_trees.txt'
file = open(file_path).read()
phylome_id = file_path.rsplit('/', 1)[1].split('_', 1)[0]
phylome_no = int(file_path.rsplit('/', 1)[1].split('_', 1)[0])

# Read protein file and get its protein
prot_df = pd.read_csv('03_calc_dists/data/0003_all_protein_names.txt', sep='\t')

# For starts for each tree
for tree_df in file.split('\n'):
    tree = tree_df.split('\t')
    t = ete3.PhyloTree(tree[-1], sp_naming_function=get_species_tag)

    # Get the sequences list
    seqsl = t.get_leaf_names()

    if 'Phy' not in tree[0]:
        continue
    elif len(seqsl) != len(set(seqsl)):
        continue
    else:
        pass

    prot = list(prot_df[prot_df['## PhylomeDB Id'] == tree[0]]['Protein name'])

    # arreglar!
    if len(prot) > 1:
        print('Careful! Double seed in database tree')
        prot = prot[0]
    else:
        prot = prot[0]

    # Root
    og = t.get_farthest_oldest_leaf(ROOTED_PHYLOMES[phylome_no])
    ogseq = og.get_leaf_names()[0]
    t.set_outgroup(og)

    # Get duplications
    t.get_descendant_evol_events()

    # Get distances
    for seq in seqsl:
        if seq != tree[0]:
            seed_dist = t.get_distance(seq, tree[0])
        else:
            seed_dist = 'NA'

        if seq != ogseq:
            og_dist = t.get_distance(seq, ogseq)
        else:
            og_dist = 'NA'

        # if seed_dist != 'NA' and og_dist != 'NA':
        #     print('%s\t%s\t%s\t%s\t%s' % (phylome_id, prot, seq,
        #                                   og_dist, seed_dist))
