import ete3
from scipy import stats
import numpy as np
import pandas as pd


def is_rooted(tree):
    # Get the tree's root
    root = tree.get_tree_root()

    # Checks that is actually rooted
    outgroups = root.get_children()
    if len(outgroups) != 2:
        return False
    else:
        return True


def tree_stats(tree):
    '''
    Get tree branch stats

    From a tree the function retrieves basic numerical information about the
    branches lengths. First, it creates a list where stores all the root to
    tip distances, then calculates the median, the mean, the width and the sum

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

    Returns:
        dictionary: dictionary with the tree statistics

    Raises:
        Exception: description
    '''

    if not is_rooted(tree):
        og = tree.get_midpoint_outgroup()
        tree.set_outgroup(og)
    
    tree.get_descendant_evol_events()

    # Retrieving all, internal, tip and root-to-tip branch lengths lists and the supports
    brlenl = list()
    int_brlenl = list()
    supportl = list()
    tip_brlenl = list()
    r2t_distl = list()

    evoltypes = dict()
    evoltypes['S'] = 0
    evoltypes['D'] = 0

    for node in tree.traverse():
        brlenl.append(node.dist)
        if not node.is_leaf():
            supportl.append(node.support)
            int_brlenl.append(node.dist)
        else:
            tip_brlenl.append(node.dist)
            r2t_distl.append(tree.get_distance(node.name))

        if 'evoltype' in node.features:
            evoltypes[node.evoltype] += 1

    treelen = sum(brlenl)

    # Generating the output dictionary
    nodedict = dict()
    nodedict['leafno'] = len(r2t_distl)
    nodedict['spno'] = len(tree.get_species())

    if nodedict['leafno'] > 1:
        nodedict['median_r2t'] = np.median(r2t_distl)
        nodedict['mean_r2t'] = np.mean(r2t_distl)
        nodedict['var_r2t'] = np.var(r2t_distl)
        nodedict['kurt_r2t'] = stats.kurtosis(r2t_distl)
        nodedict['skew_r2t'] = stats.skew(r2t_distl)
        nodedict['median_brlens'] = np.median(brlenl)
        nodedict['mean_brlens'] = np.mean(brlenl)
        nodedict['var_brlens'] = np.var(brlenl)
        nodedict['kurt_brlens'] = stats.kurtosis(brlenl)
        nodedict['skew_brlens'] = stats.skew(brlenl)
        nodedict['median_int_brlens'] = np.median(int_brlenl)
        nodedict['mean_int_brlens'] = np.mean(int_brlenl)
        nodedict['var_int_brlens'] = np.var(int_brlenl)
        nodedict['kurt_int_brlens'] = stats.kurtosis(int_brlenl)
        nodedict['skew_int_brlens'] = stats.skew(int_brlenl)
        nodedict['median_tip_brlens'] = np.median(tip_brlenl)
        nodedict['mean_tip_brlens'] = np.mean(tip_brlenl)
        nodedict['var_tip_brlens'] = np.var(tip_brlenl)
        nodedict['kurt_tip_brlens'] = stats.kurtosis(tip_brlenl)
        nodedict['skew_tip_brlens'] = stats.skew(tip_brlenl)
        nodedict['mean_bs'] = np.mean(supportl)
        nodedict['width'] = tree.get_farthest_leaf()[1]
        nodedict['tree_length'] = treelen
        nodedict['tlen_leafno_ratio'] = nodedict['tree_length'] / nodedict['leafno']
        nodedict['S'] = evoltypes['S']
        nodedict['D'] = evoltypes['D']
        nodedict['duprate'] =  evoltypes['S'] / sum(evoltypes.values())
        if treelen != 0:
            nodedict['treeness'] = sum(int_brlenl) / treelen
        if nodedict['D'] == 0:
            nodedict['single_copy'] = True
        else:
            nodedict['single_copy'] = False

    return nodedict


# Define functions ----
def get_species(node):
    '''
    Get species name

    Args:
        node (TreeNode): tree node with sequence name

    Returns:
        string: the tree node species label
    '''

    if '_' in node:
        return node.split("_")[1]
    else:
        return node


def read_treeline(line):
    '''
    Function for reading a phylome tree line

    The function gets the line and splits it using the tab delimiter, then
    stores the tree seed, the inference model used, its likelihood and the
    tree as an ete3 PhyloTree object in a dictionary.

    Args:
        line (str): tree line

    Returns:
        dict: containing the seed, the model, the likelihood and the tree

    Raises:
        KeyError: there is not all the information in the line

    '''

    line = line.split('\t')
    tree_dict = dict()

    try:
        tree_dict['seed'] = line[0]
        tree_dict['model'] = line[1]
        tree_dict['likelihood'] = line[2]
        tree_dict['tree'] = ete3.PhyloTree(line[3],
                                           sp_naming_function=get_species)
    except KeyError:
        raise 'The tree line does not contain all the information.'

    return tree_dict


def main():
    odfl = []
    for line in open('../../01_distances/data/mammal_trees.nwk', 'r'):
        tree = read_treeline(line)
        treedict = tree_stats(tree['tree'])
        odfl.append({**{'Seed': tree['seed']}, **treedict})
        print(tree['seed'])
    
    odf = pd.DataFrame(odfl)
    odf.to_csv('../outputs/tree_stats.csv', index=False, sep='\t')

    return 0


if __name__ == '__main__':
    main()