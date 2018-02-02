import logging

from ._utils.utils import *
from ._utils.check import *
from ._utils.visualize import *


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

        df_paths = get_path_statistics(self.df_paths)

        # reconstruct linestack from swc.

        linestack_output = swc_to_linestack(self.df_swc, self.voxelsize)

        if linestack_output is not None:
            self.linestack, self.soma_on_stack, self.coordindate_padding = linestack_output
            self.df_paths = get_path_on_stack(self.df_paths, self.voxelsize, self.coordindate_padding)
            self.density_stack, self.dendritic_center = calculate_density(self.linestack, self.voxelsize)
        else:
            self.linestack = None

    def summary(self):

        """
        The summary of the cell morphology. 
        """
        
        df_paths = self.df_paths.copy()
        
        soma = df_paths[df_paths.type == 1]
        axon = df_paths[df_paths.type == 2]
        dend_basal = df_paths[df_paths.type == 3]
        dend_apical = df_paths[df_paths.type == 4]
        
        axon_summary = get_summary_of_paths(axon)
        dend_basal_summary = get_summary_of_paths(dend_basal)
        dend_apical_summary = get_summary_of_paths(dend_apical)

        labels = [
                'type',
                'num_path_segments',
                'num_branchpoints',
                'num_irreducible_nodes',
                'max_branch_order',
                'average_nodal_angle_deg',
                'average_nodal_angle_rad',
                'average_local_angle_deg',
                'average_local_angle_rad',
                'average_tortuosity',
                'real_length_sum',
                'real_length_mean',
                'real_length_median',
                'real_length_min',
                'real_length_max',
                'euclidean_length_sum',
                'euclidean_length_mean',
                'euclidean_length_median',
                'euclidean_length_min',
                'euclidean_length_max',
                ]
        
        neurites = [axon_summary,dend_basal_summary,dend_apical_summary]
        df_summary = pd.DataFrame.from_records([n for n in neurites if n is not None], columns=labels)
            
        return df_summary

    # def summary(self, save_to=None,  print_results=True):

    #     """
    #     Print out summary of the cell morphology. 

    #     Parameters
    #     ----------
    #     save_to: str
    #         path and filename of the output json file. 

    #     print_results: bool
    #         if True, print out summary using logging.
    #     """

    #     # branch order / number of branch points

    #     branchpoints = np.vstack(self.df_paths.connect_to_at)
    #     branchpoints = unique_row(branchpoints)
    #     num_branchpoints = len(branchpoints)

    #     max_branch_order = max(self.df_paths.branch_order)

    #     terminalpaths = self.df_paths.path[self.df_paths.connected_by.apply(len) == 0].as_matrix()
    #     terminalpoints = np.vstack([p[-1] for p in terminalpaths])
    #     num_terminalpoints = len(terminalpoints)

    #     outerterminals = get_outer_terminals(terminalpoints)

    #     num_irreducible_nodes = num_branchpoints + num_terminalpoints

    #     num_dendritic_segments = len(self.df_paths)

    #     # path length

    #     dendritic = self.df_paths['real_length']
    #     dendritic_sum = dendritic.sum()
    #     dendritic_mean = dendritic.mean()
    #     dendritic_median = dendritic.median()
    #     dendritic_min = dendritic.min()
    #     dendritic_max = dendritic.max()

    #     euclidean = self.df_paths['euclidean_length']
    #     euclidean_sum = euclidean.sum()
    #     euclidean_mean = euclidean.mean()
    #     euclidean_median = euclidean.median()
    #     euclidean_min = euclidean.min()
    #     euclidean_max = euclidean.max()

    #     tortuosity = dendritic / euclidean
    #     average_tortuosity = np.mean(tortuosity)      

    #     # node angles
    #     average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad = get_average_angles(self.df_paths)

    #     summary = {
    #         "general": {
    #             "number_of_dendritic_segments": int(num_dendritic_segments),
    #             "number_of_branch_points": int(num_branchpoints),
    #             "number_of_irreducible_nodes": int(num_irreducible_nodes),
    #             "max_branch_order": int(max_branch_order),
    #             # "max_strahler_order": int(max_strahler_order),
    #         },
    #         "angle": {
    #             "average_nodal_angle_in_degree": average_nodal_angle_deg,
    #             "average_nodal_angle_in_radian": average_nodal_angle_rad,
    #             "average_local_angle_in_degree": average_local_angle_deg,
    #             "average_local_angle_in_radian": average_local_angle_rad,
    #         },
    #         "length": {
    #             "tortuosity": average_tortuosity,
    #             "dendritic": {
    #                 "sum": dendritic_sum,
    #                 "mean": dendritic_mean,
    #                 "median": dendritic_median,
    #                 "min": dendritic_min,
    #                 "max": dendritic_max,
    #             },
    #             "euclidean":{
    #                 "sum": euclidean_sum,
    #                 "mean": euclidean_mean,
    #                 "median": euclidean_median,
    #                 "min": euclidean_min,
    #                 "max": euclidean_max,
    #             }
    #         }

    #     }

    #     # dendritic tree density
    #     if self.linestack is not None:

    #         if self.voxelsize is None:
    #             voxelsize = [1,1,1]
    #         else:
    #             voxelsize = self.voxelsize

    #         dendritic_center = self.dendritic_center

    #         soma_coordinates = self.soma_on_stack

    #         asymmetry = np.sqrt(((soma_coordinates[:2] - dendritic_center)**2).sum())

    #         all_termianls_to_dendritic_center = np.sqrt(np.sum((terminalpoints[:,:2] - dendritic_center) ** 2,  1)) * voxelsize[0]
    #         out_terminals_to_dendritic_center = np.sqrt(np.sum((outerterminals[:, :2]- dendritic_center) **2, 1)) * voxelsize[0]

    #         typical_radius = np.mean(all_termianls_to_dendritic_center)
    #         outer_radius = np.mean(out_terminals_to_dendritic_center)

    #         import cv2
    #         dendritic_area = cv2.contourArea(outerterminals[:, :2].astype(np.float32)) * voxelsize[0]**2/1000

    #         summary.update({"density": {
    #             "asymmetry": asymmetry,
    #             "outer_radius": outer_radius,
    #             "typical_radius": typical_radius,
    #             "dendritic_area": dendritic_area
    #         }})

    #     # print and save
    #     if print_results:
    #         print_summary(summary)

    #     if save_to is not None:

    #         logging.info('  Writing to {}'.format(save_to))

    #         import json
    #         with open(save_to, 'w') as f:
    #             json.dump(summary, f)

    #     # return summary


    def show_threeviews(self, save_to=None):

        """
        Plot cell morphology in three views (Top and two sides).

        Parameters
        ----------
        save_to: str
            Path the figure saved to. e.g. "./figure/threeviews.png"

        """


        import matplotlib.pyplot as plt
        from matplotlib_scalebar.scalebar import ScaleBar

        plt.figure(figsize=(16,16))
        ax1 = plt.subplot2grid((4,4), (0,1), rowspan=3, colspan=3)
        ax2 = plt.subplot2grid((4,4), (0,0), rowspan=3, colspan=1)
        ax3 = plt.subplot2grid((4,4), (3,1), rowspan=1, colspan=3)
        ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1, colspan=1)  

        df_paths = self.df_paths
        dendrites = df_paths[df_paths.type != 1]
        soma = df_paths[df_paths.type == 1].path[0][0]

        lims = find_lims(dendrites)

        plot_skeleton(ax2, dendrites, soma, 2, 0, lims)
        plot_skeleton(ax3, dendrites, soma, 1, 2, lims)
        plot_skeleton(ax1, dendrites, soma, 1, 0, lims)
        scalebar = ScaleBar(1, units=self.unit, location='lower left', box_alpha=0)
        ax1.add_artist(scalebar)    
        ax4.axis('off')

        if save_to is not None:
            plt.savefig(save_to)
    
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