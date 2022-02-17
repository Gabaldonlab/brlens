import ete3
from multiprocessing import Process, Manager
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf


def subtree_plot(phylome_id, prot_dict, seed_id, tree):
    mphsptrees = dict()
    mphgtrees = dict()
    for i, subtree in enumerate(tree.traverse()):
        sps = [sp.split('_')[1] for sp in subtree.get_leaf_names()]
        if (len(sps) == len(set(sps)) and len(sps) > 10 and
                subtree.get_farthest_leaf()[1] != 0):
            mphsptrees[len(sps)] = subtree
        elif len(sps) > 10 and subtree.get_farthest_leaf()[1] != 0:
            mphgtrees[len(sps) - len(set(sps))] = subtree

    if len(mphsptrees) > 0:
        rtree = mphsptrees.get(max(mphsptrees.keys()))
    elif len(mphgtrees) > 0:
        rtree = mphgtrees.get(min(mphgtrees.keys()))
    else:
        rtree = tree

    root = rtree.get_common_ancestor(rtree)

    rtldist = list()
    for leaf in rtree.get_leaf_names():
        rtldist.append(rtree.get_distance(root, leaf))

    odict = {'st_%s' % seed_id: rtldist}

    return odict


def mrca_tt_plot(tree, seed_id):
    ln = tree.get_leaf_names()

    distl = list()
    for seq_from in ln:
        for seq_to in ln:
            if seq_from != seq_to:
                subt = tree.get_common_ancestor(seq_from, seq_to)
                fromd = subt.get_distance(seq_from)
                tod = subt.get_distance(seq_to)
                if (fromd != tree.get_distance(seq_from) and
                        tod != tree.get_distance(seq_to) and
                        subt.evoltype == 'S'):
                    distl.extend([fromd, tod])

    odict = {'mrca_%s' % seed_id: distl}

    return odict


def root_tt_plot(tree, seed_id):
    distl = list()
    for leaf in tree.get_leaf_names():
        distl.append(tree.get_distance(leaf))

    odict = {'root_%s' % seed_id: distl}

    return odict


class dist_process(Process):
    def __init__(self, tree_row, phylome_id, prot_dict, olist, plots):
        Process.__init__(self)
        self.tree_row = tree_row
        self.phylome_id = phylome_id
        self.prot_dict = prot_dict
        self.olist = olist
        self.plots = plots

    def run(self):
        tree = self.tree_row.split('\t')
        t = ete3.PhyloTree(tree[3], sp_naming_function=get_species_tag)

        if (len(t.get_species()) > 10 and
                len(t.get_leaf_names()) < 3 * len(t.get_species())):
            print('Calculating: %s, species no.: %s, leaves no.: %s' %
                  (tree[0], len(t.get_species()), len(t.get_leaf_names())))
            root(t, root_dict[int(self.phylome_id)])
            t.get_descendant_evol_events()

            st_ref = subtree_tt_ref(t)
            mrca_ref = mrca_tt_ref(t)
            root_ref = root_tt_ref(t)

            tnames = t.get_leaf_names()

            self.plots.append(mrca_tt_plot(t, tree[0]))
            self.plots.append(root_tt_plot(t, tree[0]))
            self.plots.append(subtree_plot(self.phylome_id, self.prot_dict,
                                           tree[0], t))

            for i, from_sp in enumerate(tnames):
                for to_sp in tnames[i + 1:]:
                    if from_sp != to_sp:
                        leaf_dist = get_dists(t, from_sp, to_sp, tree[0],
                                              self.phylome_id, self.prot_dict,
                                              st_ref, mrca_ref, root_ref)
                        if leaf_dist is not None:
                            self.olist.append(leaf_dist)


def main():
    phylome_id = ifile.rsplit('/', 1)[1].split('_', 1)[0]
    file_id = ifile.rsplit('/', 1)[1].split('.', 1)[0]

    dist_fn = '/'.join([odir, (file_id + '_dist.csv')])

    if not file_exists(dist_fn):
        print('Creating: ', dist_fn)

        create_folder(odir)

        trees = open(ifile, 'r').read().split('\n')
        prot_dict = csv_to_dict(prots, '\t')

        with Manager() as manager:
            olist = manager.list()
            plots = manager.list()

            processes = list()
            for tree_row in trees:
                if tree_row != '':
                    if len(processes) >= cpus:
                        done = False
                        while not done:
                            for process in processes:
                                if not process.is_alive():
                                    processes.remove(process)
                                    done = True

                    process = dist_process(tree_row, phylome_id,
                                           prot_dict, olist, plots)
                    processes.append(process)
                    process.start()

            for process in processes:
                process.join()

            # Writing output files
            odf = pd.DataFrame(list(olist))
            odf.to_csv(dist_fn, index=False)

            # Plotting densities
            mrca_fn = '/'.join([odir, (file_id + '_mrca.pdf')])
            st_fn = '/'.join([odir, (file_id + '_st.pdf')])
            root_fn = '/'.join([odir, (file_id + '_root.pdf')])

            mrca_pdf = backend_pdf.PdfPages(mrca_fn)
            st_pdf = backend_pdf.PdfPages(st_fn)
            root_pdf = backend_pdf.PdfPages(root_fn)

            for plotl in list(plots):
                plotdf = pd.DataFrame(plotl)
                tree_id = list(plotl.keys())[0]
                plotdf.plot.density(color='darkorange',
                                    legend=False)
                plt.title('Density plot for %s' % tree_id)
                plt.axvline(plotdf[tree_id].mean(), color='k',
                            linestyle='dashed', linewidth=1)
                plt.axvline(plotdf[tree_id].median(), color='blue',
                            linestyle='dashed', linewidth=1)

                if 'mrca_' in tree_id:
                    plt.xlabel('MRCA-to-tip distance')
                    mrca_pdf.savefig()
                elif 'st_' in tree_id:
                    plt.xlabel('Subtree MRCA to tip distance')
                    st_pdf.savefig()
                elif 'root_' in tree_id:
                    plt.xlabel('Root to tip distance')
                    root_pdf.savefig()

            mrca_pdf.close()
            st_pdf.close()
            root_pdf.close()
