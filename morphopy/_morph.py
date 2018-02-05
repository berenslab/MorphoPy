import logging

from ._utils.utils import *
from ._utils.check import *
from ._utils.visualize import *
from ._utils.summarize import *


__all__ = ['Morph']

class Morph(object):

    def __init__(self, data, voxelsize=None, loglevel='info'):

        """
        Initialize Morph object. Load swc as Pandas DataFrame (df_swc). Split all paths on branch point and save as
        df_paths, related information (connection, path length, branch order etc.) are calculated. Other meta data are
        also saved into Morph Object. If voxelszie is provided, a linestack is constructed and dendritic tree density
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
            num_dendritic_segments = summary['num_path_segments'][n]
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
            logging.info('    Number of dendritic arbor segment: {}'.format(num_dendritic_segments))
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

    def show_morph(self, view='xy', plot_axon=True, plot_dendrites=True):
        
        df_paths = self.df_paths.copy()
        fig, ax = plt.subplots(1, 1, figsize=(12,12))
        ax = plot_morph(ax, df_paths, view, plot_axon, plot_dendrites)

        return fig, ax
    
    def show_threeviews(self, plot_axon=True, plot_dendrites=True):
        
        df_paths = self.df_paths
        
        fig, ax = plt.subplots(1, 3, figsize=(18,6))
        
        ax0 = plot_morph(ax[0], df_paths, 'xy', plot_axon, plot_dendrites)
        ax1 = plot_morph(ax[1], df_paths, 'xz', plot_axon, plot_dendrites)
        ax2 = plot_morph(ax[2], df_paths, 'yz', plot_axon, plot_dendrites)
    
    def show_animation(self):
        
        import matplotlib.animation as animation
        from IPython.display import HTML
        
        df_paths = self.df_paths
        
        soma = df_paths[df_paths.type == 1].path[0][0]
        axon = df_paths[df_paths.type == 2]
        dendrites = df_paths[df_paths.type == 3]
            
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection='3d')

        def init():
            
            ax.scatter(0,0,0, s=280, color='grey')

            dcolors_idx = np.linspace(0, 255, max(df_paths.branch_order)+1).astype(int)
            dcolors = np.vstack(plt.cm.Reds_r(dcolors_idx))[:, :3]

            for row in dendrites.iterrows():

                path_id = row[0]
                path = row[1]['path'] - soma
                order = row[1]['branch_order']
                bpt = path[0]     

                ax.plot(path[:, 0], path[:, 1], path[:, 2], color=dcolors[int(order)])
                ax.scatter(bpt[0], bpt[1], bpt[2], color=dcolors[int(order)], zorder=1)

            acolors_idx = np.linspace(0, 255, max(df_paths.branch_order)+1).astype(int)
            acolors = np.vstack(plt.cm.Blues_r(dcolors_idx))[:, :3]

            for row in axon.iterrows():

                path_id = row[0]
                path = row[1]['path'] - soma
                order = row[1]['branch_order']
                bpt = path[0]     

                ax.plot(path[:, 0], path[:, 1], path[:, 2], color=acolors[int(order)])
                ax.scatter(bpt[0], bpt[1], bpt[2], color=acolors[int(order)], zorder=1)

            lim_max = int(np.ceil((np.vstack(df_paths.path.as_matrix()) - soma).max() / 20) * 20)
            # lim_min = int(np.floor(np.vstack(self.df_paths.path.as_matrix()).min() / 20) * 20)

            ax.set_xlim(-lim_max, lim_max)
            ax.set_ylim(-lim_max, lim_max) 
            ax.set_zlim(-lim_max, lim_max) 
            
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

    # def show_threeviews(self, save_to=None):

    #     """
    #     Plot cell morphology in three views (Top and two sides).

    #     Parameters
    #     ----------
    #     save_to: str
    #         Path the figure saved to. e.g. "./figure/threeviews.png"

    #     """

    #     import matplotlib.pyplot as plt
    #     from matplotlib_scalebar.scalebar import ScaleBar

    #     plt.figure(figsize=(16,16))
    #     ax1 = plt.subplot2grid((4,4), (0,1), rowspan=3, colspan=3)
    #     ax2 = plt.subplot2grid((4,4), (0,0), rowspan=3, colspan=1)
    #     ax3 = plt.subplot2grid((4,4), (3,1), rowspan=1, colspan=3)
    #     ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1, colspan=1)  

    #     df_paths = self.df_paths
    #     dendrites = df_paths[df_paths.type != 1]
    #     soma = df_paths[df_paths.type == 1].path[0][0]

    #     lims = find_lims(dendrites)

    #     plot_skeleton(ax2, dendrites, soma, 2, 0, lims)
    #     plot_skeleton(ax3, dendrites, soma, 1, 2, lims)
    #     plot_skeleton(ax1, dendrites, soma, 1, 0, lims)
    #     scalebar = ScaleBar(1, units=self.unit, location='lower left', box_alpha=0)
    #     ax1.add_artist(scalebar)    
    #     ax4.axis('off')

    #     if save_to is not None:
    #         plt.savefig(save_to)
    
    # def show_density(self): 

    #     """
    #     Plot cell morphology on dendritic density map.
    #     """

    #     try:
    #         density_stack = self.density_stack
    #         voxelsize = self.voxelsize
    #     except:
    #         logging.info('No density stack. Please provide the voxel sizes of the `.swc` file.')
    #         return None
        
    #     import matplotlib.pyplot as plt
    #     from matplotlib_scalebar.scalebar import ScaleBar

    #     linestack = self.linestack
    #     dendritic_center = self.dendritic_center
    #     soma_on_stack = self.soma_on_stack
        
    #     plt.figure(figsize=(16, 16))
    #     plt.imshow(density_stack.sum(2), cmap=plt.cm.gnuplot2_r, origin='lower')
    #     plt.scatter(dendritic_center[1], dendritic_center[0], color='g', marker='*', s=180, label='Dendritic Center')
    #     plt.scatter(soma_on_stack[1], soma_on_stack[0], color='r',  marker='*', s=180, label='Soma')
        
    #     linestack_xy = linestack.sum(2)
    #     linestack_xy[linestack_xy !=0] = 1
    #     linestack_xy = np.ma.masked_array(linestack_xy, ~linestack.any(2))
    #     plt.imshow(linestack_xy, origin='lower', cmap=plt.cm.binary)
        
    #     plt.legend(frameon=False)

    #     scalebar = ScaleBar(voxelsize[0], units=self.unit, location='lower left', box_alpha=0)
    #     plt.gca().add_artist(scalebar)   
        
    #     plt.axis('off')