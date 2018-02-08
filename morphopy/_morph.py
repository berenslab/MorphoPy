import logging

from ._utils.utils import *
from ._utils.check import *
from ._utils.visualize import *
from ._utils.summarize import *
from ._utils.representation import get_persistence_barcode

__all__ = ['Morph']

class Morph(object):

    def __init__(self, data, voxelsize=None, loglevel='info'):

        """
        Initialize Morph object. Load swc as Pandas DataFrame (df_swc). Split all paths on branch point and save as
        df_paths, related information (connection, path length, branch order etc.) are calculated. Other meta data are
        also saved into Morph Object. If voxelsize is provided, a linestack is constructed and dendritic tree density
        is computed.

        Parameters
        ----------
        data: str
            path to the `.swc` file.
        voxelsize: list or array-like
            specify the voxel separation. e.g. [0.665, 0.665, 1].
            If provided, linestack is reconstructed and dendritic tree density map will be computed.
        loglevel: str
            'debug', 'info', 'warning', 'error', 'critical'.
        """

        # logging
        self._logger = get_logger(loglevel)

        # meta data
        self.unit = 'um'
        self.voxelsize = voxelsize

        logging.info('  SWC file: {}\n'.format(data))
        logging.info('  unit: {}'.format(self.unit))
        logging.info('  voxel size: {}\n'.format(self.voxelsize))

        # load data
        G, df_swc = read_swc(data)

        # check data
        logging.info('  ===================  ')
        logging.info('  Checking `.swc`...   \n')

        check_swc(df_swc)

        # split swc into soma, dendrites, axon, etc..
        df_paths = get_df_paths(G)
        df_paths = check_path_connection(df_paths) # find which paths connect to which

        self.df_swc = df_swc
        self.df_paths = df_paths
        self.G = G
        self.df_persistence_barcode = None

    def processing(self):

        """
        Further processing df_paths and get statistics info such as path lengths, branching order into DataFrame.

        Linestack is reconstructed if needed data is given (voxel size of the original image).
        Then dendritic density is calculated based on linestack.
        """

        self.df_paths = get_path_statistics(self.df_paths)
        self.summary_data = get_summary_data(self.df_paths)

        # reconstruct linestack from swc.

        linestack_output = swc_to_linestack(self.df_swc, self.voxelsize)

        if linestack_output is not None:
            self.linestack, self.soma_on_stack, self.coordindate_padding = linestack_output
            self.df_paths = get_path_on_stack(self.df_paths, self.voxelsize, self.coordindate_padding)
            self.density_stack, self.dendritic_center = calculate_density(self.linestack, self.voxelsize)
        else:
            self.linestack = None

    def show_summary(self):

        """
        Print out summary statistics of the cell.

        Parameters
        ----------
        summary: pd.DataFrame
            a pandas DataFrame that contains summary of one type of neurites of the cell.
        """

        logging.info('  Summary of the cell')
        logging.info('  ======================\n')

        summary = self.summary_data.to_dict()

        for n in range(len(summary['type'])):

            neurite_type = summary['type'][n]
            num_path_segments = summary['num_path_segments'][n]
            num_branchpoints = summary['num_branchpoints'][n]
            num_irreducible_nodes = summary['num_irreducible_nodes'][n]
            max_branch_order = summary['max_branch_order'][n]
            average_nodal_angle_deg = summary['average_nodal_angle_deg'][n]
            average_nodal_angle_rad = summary['average_nodal_angle_rad'][n]
            average_local_angle_deg = summary['average_local_angle_deg'][n]
            average_local_angle_rad = summary['average_local_angle_rad'][n]
            average_tortuosity = summary['average_tortuosity'][n]
            real_length_sum = summary['real_length_sum'][n]
            real_length_mean = summary['real_length_mean'][n]
            real_length_median = summary['real_length_median'][n]
            real_length_min = summary['real_length_min'][n]
            real_length_max = summary['real_length_max'][n]
            euclidean_length_sum = summary['euclidean_length_sum'][n]
            euclidean_length_mean = summary['euclidean_length_mean'][n]
            euclidean_length_median = summary['euclidean_length_median'][n]
            euclidean_length_min = summary['euclidean_length_min'][n]
            euclidean_length_max = summary['euclidean_length_max'][n]

            logging.info('  {}\n'.format(neurite_type).upper())
            logging.info('    Number of arbor segments: {}'.format(num_path_segments))
            logging.info('    Number of branch points: {}'.format(num_branchpoints))
            logging.info('    Number of irreducible nodes: {}'.format(num_irreducible_nodes))
            logging.info('    Max branching order: {}\n'.format(max_branch_order))
            # logging.info('    Max Strahler order: {}\n\n'.format(max_strahler_order))

            logging.info('  ## Angle \n')
            logging.info('    Average nodal angle in degree: {:.3f}'.format(average_nodal_angle_deg))
            logging.info('    Average nodal angle in radian: {:.3f}'.format(average_nodal_angle_rad))
            logging.info('    Average local angle in degree: {:.3f}'.format(average_local_angle_deg))
            logging.info('    Average local angle in radian: {:.3f} \n'.format(average_local_angle_rad))

            logging.info('  ## Average tortuosity: {:.3f}\n'.format(average_tortuosity))

            logging.info('  ## Real length (μm)\n')
            # logging.info('  ## Dendritic length\n')
            logging.info('    Sum: {:.3f}'.format(real_length_sum))
            logging.info('    Mean: {:.3f}'.format(real_length_mean))
            logging.info('    Median: {:.3f}'.format(real_length_median))
            logging.info('    Min: {:.3f}'.format(real_length_min))
            logging.info('    Max: {:.3f}\n'.format(real_length_max))

            logging.info('  ## Euclidean length (μm)\n')

            logging.info('    Sum: {:.3f}'.format(euclidean_length_sum))
            logging.info('    Mean: {:.3f}'.format(euclidean_length_mean))
            logging.info('    Median: {:.3f}'.format(euclidean_length_median))
            logging.info('    Min: {:.3f}'.format(euclidean_length_min))
            logging.info('    Max: {:.3f}\n'.format(euclidean_length_max))

            logging.info('  ======================\n')

    def show_morph(self, view='xy', plot_axon=True, plot_basal_dendrites=True, plot_apical_dendrites=True):

        df_paths = self.df_paths.copy()
        fig, ax = plt.subplots(1, 1, figsize=(12,12))
        ax = plot_morph(ax, df_paths, view, plot_axon, plot_basal_dendrites, plot_apical_dendrites)

        return fig, ax

    def show_threeviews(self, plot_axon=True, plot_basal_dendrites=True, plot_apical_dendrites=True):

        df_paths = self.df_paths.copy()

        fig, ax = plt.subplots(1, 3, figsize=(18,6))

        ax0 = plot_morph(ax[0], df_paths, 'xy', plot_axon, plot_basal_dendrites, plot_apical_dendrites)
        ax1 = plot_morph(ax[1], df_paths, 'xz', plot_axon, plot_basal_dendrites, plot_apical_dendrites)
        ax2 = plot_morph(ax[2], df_paths, 'yz', plot_axon, plot_basal_dendrites, plot_apical_dendrites)

        return fig, [ax0, ax1, ax2]

    def show_animation(self):

        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.animation as animation
        from IPython.display import HTML

        df_paths = self.df_paths.copy()

        soma = df_paths[df_paths.type == 1].path[0][0]
        axon = df_paths[df_paths.type == 2]
        basal_dendrites = df_paths[df_paths.type == 3]
        apical_dendrites = df_paths[df_paths.type == 4]

        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection='3d')

        def init():

            # soma
            ax.scatter(0,0,0, s=280, color='grey')

            # basal dendrites

            if len(basal_dendrites)>0:
                bdcolors_idx = np.linspace(0, 200, max(basal_dendrites.branch_order)+1).astype(int)
                bdcolors = np.vstack(plt.cm.Reds_r(bdcolors_idx))[:, :3]

                for row in basal_dendrites.iterrows():

                    path_id = row[0]
                    path = row[1]['path'] - soma
                    order = row[1]['branch_order']
                    bpt = path[0]

                    ax.plot(path[:, 0], path[:, 1], path[:, 2], color=bdcolors[int(order)])
                    ax.scatter(bpt[0], bpt[1], bpt[2], color=bdcolors[int(order)], zorder=1)

            # apical dendrites
            if len(apical_dendrites)>0:
                adcolors_idx = np.linspace(0, 200, max(apical_dendrites.branch_order)+1).astype(int)
                adcolors = np.vstack(plt.cm.Purples_r(adcolors_idx))[:, :3]

                for row in apical_dendrites.iterrows():

                    path_id = row[0]
                    path = row[1]['path'] - soma
                    order = row[1]['branch_order']
                    bpt = path[0]

                    ax.plot(path[:, 0], path[:, 1], path[:, 2], color=adcolors[int(order)])
                    ax.scatter(bpt[0], bpt[1], bpt[2], color=adcolors[int(order)], zorder=1)

            # axon
            if len(axon)>0:
                acolors_idx = np.linspace(0, 200, max(axon.branch_order)+1).astype(int)
                acolors = np.vstack(plt.cm.Blues_r(acolors_idx))[:, :3]

                for row in axon.iterrows():

                    path_id = row[0]
                    path = row[1]['path'] - soma
                    order = row[1]['branch_order']
                    bpt = path[0]

                    ax.plot(path[:, 0], path[:, 1], path[:, 2], color=acolors[int(order)])
                    ax.scatter(bpt[0], bpt[1], bpt[2], color=acolors[int(order)], zorder=1)

            lim_max = int(np.ceil((np.vstack(df_paths.path.as_matrix()) - soma).max() / 20) * 20)
            lim_min = int(np.floor((np.vstack(df_paths.path.as_matrix()) - soma).min() / 20) * 20)

            lim = max(abs(lim_max), abs(lim_min))

            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_zlim(-lim, lim)

            return fig,

        def animate(i):
            # azimuth angle : 0 deg to 360 deg
            ax.view_init(elev=10, azim=i*4)
            return fig,

        logging.info('  Generating animation. It might take some time...')

        ani = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=90, interval=150, blit=True)

        plt.close()

        return HTML(ani.to_html5_video())

    def show_persistence_diagram(self, axon=True, basal_dendrites=True, apical_dendrites=True):
        """
        Plots the persistence diagram of the neuron. Persistence is a concept from topology that defines invariant
        structures. Its clearer definition can be found in
         - S. Chepushtanova, T. Emerson, E.M. Hanson, M. Kirby, F.C. Motta, R. Neville, C. Peterson, P.D. Shipman, and
L. Ziegelmeier. Persistence images: An alternative persistent homology representation. CoRR, abs/1507.06217, 2015.
         - 	arXiv:1603.08432
         - Li, Yanjie, et al.
        "Metrics for comparing neuronal tree shapes based on persistent homology." PloS one 12.8 (2017): e0182184.


        :param axon: boolean (default True). When set to False the axonal branches are excluded.
        :param basal_dendrites: boolean (default True). When set to False the basal dendritic branches are excluded.
        :param apical_dendrites: boolean (default True). When set to False the apical dendritic branches are excluded.
        :return: fig, ax
        """

        if self.df_persistence_barcode is None:
            self.df_persistence_barcode = get_persistence_barcode(self.G)

        index = (self.df_persistence_barcode.type == 1)

        if axon:
            index |= self.df_persistence_barcode.type == 2
        if basal_dendrites:
            index |= self.df_persistence_barcode.type == 3
        if apical_dendrites:
            index |= self.df_persistence_barcode.type == 4

        plotting_data = self.df_persistence_barcode[index]

        fig, ax = plt.subplots(1, 3, figsize=(12, 4))

        plot_persistence_diagram(plotting_data, ax[0])
        plot_persistence_image_2d(plotting_data, ax[1])
        plot_persistence_image_1d(plotting_data, ax[2])

        return fig, ax


