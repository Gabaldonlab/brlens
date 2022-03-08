import pandas as pd


def annotate_lineages(tree, gnmdf, col):
    for leaf in tree.get_leaves():
        sp = list(leaf.get_species())[0]
        feat_val = list(gnmdf[col][gnmdf['PhylomeDB'] == sp])[0]
        leaf.add_feature(col, feat_val)

    return 0


def clade_norm(tree, cladedf, col):
    annotate_lineages(tree, cladedf, col)
    tlno = len(tree.get_leaf_names())

    mphysets = list()
    for st in tree.traverse():
        leaves = st.get_leaves()

        feat_list = list()
        for stleag in st.get_leaves():
            feat_list.append(getattr(stleag, col))

        if len(set(feat_list)) == 1 and len(leaves) > 1 and len(leaves) != tlno:
            mphy = dict()
            mphy['node'] = st
            mphy['seq_no'] = len(leaves)
            mphy[col] = feat_list[0]
            mphysets.append(mphy)

    nodedf = pd.DataFrame(mphysets)
    indexes = nodedf.index

    mphylist = list()
    for group in set(nodedf[col]):
        maxval = max(nodedf.loc[nodedf[col] == group][col])
        dfindex = indexes[(nodedf[col] == group) & (nodedf[col] == maxval)]
        mphylist.append(mphysets[dfindex.values[0]])

    return mphylist
