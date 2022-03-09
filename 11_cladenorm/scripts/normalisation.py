import pandas as pd
import ete3
import numpy as np


def annotate_lineages(tree, gnmdf, col):
    for leaf in tree.get_leaves():
        sp = list(leaf.get_species())[0]
        feat_val = list(gnmdf[col][gnmdf['Proteome'] == sp])[0]
        leaf.add_feature(col, feat_val)

    return 0


def clade_norm(tree, tree_id, cladedf, col):
    annotate_lineages(tree, cladedf, col)
    tlno = len(tree.get_leaf_names())

    mphysets = list()
    for st in tree.traverse():
        leaves = st.get_leaves()

        stwdth = st.get_farthest_leaf()[1]
        print(stwdth)

        feat_list = list()
        for stleag in st.get_leaves():
            feat_list.append(getattr(stleag, col))

        if (len(set(feat_list)) == 1 and len(leaves) > 1 and
                len(leaves) != tlno and st.evoltype == 'S' and stwdth != 0):
            mphy = dict()
            mphy['tree'] = tree_id
            mphy['node'] = st
            mphy['seq_no'] = len(leaves)
            mphy[col] = feat_list[0]
            mphysets.append(mphy)

    nodedf = pd.DataFrame(mphysets)
    indexes = nodedf.index

    mphylist = list()
    for group in set(nodedf[col]):
        if str(group) != 'nan':
            maxval = max(nodedf.loc[nodedf[col] == group]['seq_no'])
            dfindex = indexes[(nodedf[col] == group) &
                              (nodedf['seq_no'] == maxval)]
            mphylist.append(mphysets[dfindex.values[0]])

    nodel = list()
    for nodedict in mphylist:
        nodel.append(norm_factor(nodedict))

    return nodel


def norm_factor(nodedict):
    ndlf = nodedict['node'].get_leaves()
    distl = list()
    for leaf in ndlf:
        distl.append(nodedict['node'].get_distance(leaf))
    nodedict['mrca_to_tip_median'] = np.median(distl)

    return nodedict
