import os
import pickle
import logging

from ._utils import *

__all__ = ['Morph']

class Morph(object):

    def __init__(self, data, unit='um', voxelsize=None, threshold=30, loglevel='INFO'):

        # logging

        self._logger = logging.getLogger()
        
        if loglevel == 'INFO':
            self._logger.setLevel(logging.INFO)
            self._logger.info('  Logging: ON')
        elif loglevel == 'DEBUG':
            self._logger.setLevel(logging.DEBUG)
            self._logger.info('  Logging: ON')
        elif loglevel == 'WARNING':
            self._logger.setLevel(logging.WARNING)
        elif loglevel == 'ERROR':
            self._logger.setLevel(logging.ERROR)
        elif loglevel == 'CRITICAL':
            self._logger.setLevel(logging.CRITICAL)
        else:
            self._logger.setLevel(logging.INFO)
            logging.info('  Please enter a valid logging mode (DEBUG, INFO, WARNING, ERROR, CRITICAL).')
            self._logger.setLevel(logging.ERROR)

        # meta data

        self.voxelsize = voxelsize
        self.unit_swc = unit
        filetype = data.split('/')[-1].split('.')[-1]
        
        # load data
        if filetype == 'swc':
            df_paths, soma, unit_df = read_swc(data, unit, voxelsize)
        else:
            logging.info('`.{}` is not supported yet.'.format(filetype))
            return None

        df_paths = update_df_paths(df_paths, soma) # find which paths connnect to which
        paths_to_fix = detect_messiness(df_paths, threshold) # fix some rare irregularities.
        df_paths = clean_messiness(df_paths, paths_to_fix)
        df_paths = get_better_paths(df_paths, soma) 
        df_paths = cleanup_better_paths(df_paths)
        df_paths = get_path_statistics(df_paths)
        df_paths = get_sorder(df_paths) # get schaler order

        self.df_paths = df_paths
        self.unit_df = unit_df
        self.soma = soma
        

        # linestack | pixel

        linestack_output = swc_to_linestack(data, self.unit_swc, voxelsize)

        if linestack_output is not None:
            self.linestack, self.soma_on_stack, self.coordindate_padding = linestack_output
            self.df_paths = get_path_on_stack(self.df_paths, self.voxelsize, self.coordindate_padding)
            self.density_stack, self.dendritic_center = calculate_density(self.linestack, voxelsize)
        else:
            self.linestack = None


    def to_swc(self, filename='morph.swc', save_to='./'):
    
        logging.info('  Writing to {}'.format(save_to + '/' + filename))

        soma_xyz = self.df_paths.loc[0].path[0]
        soma = np.array([1, 1, soma_xyz[0], soma_xyz[1], soma_xyz[2], 0, -1])

        swc_arr = soma.reshape(1,len(soma))

        for row in self.df_paths.iterrows():    
            path_id = row[0]
            path = row[1]['path'][1:]

            connect_to = row[1]['connect_to']
            connect_to_at = row[1]['connect_to_at']

            swc_path = np.column_stack([np.ones(len(path)) * 3, path])
            swc_path = np.column_stack([np.arange(len(swc_arr)+1, len(path)+len(swc_arr)+1), swc_path])
            swc_path = np.column_stack([swc_path, np.zeros(len(path))])
            swc_path = np.column_stack([swc_path, swc_path[:, 0]-1])

            swc_path[0][-1] = np.where((swc_arr[:, 2:5] == np.array(connect_to_at)).all(1))[0]+1

            swc_arr = np.vstack([swc_arr, swc_path])

        df_swc = pd.DataFrame(swc_arr)
        df_swc.index = np.arange(1, len(df_swc)+1)
        df_swc.columns = [['ID', 'Type', 'x', 'y', 'z', 'Raidus', 'PID']]
        df_swc[['ID', 'Type', 'PID']] = df_swc[['ID', 'Type', 'PID']].astype(int)

        df_swc.to_csv(save_to + '/' + filename, sep=' ', index=None, header=None)


    def summary(self, save_to=None,  print_results=True):

        # branch order / number of branch points

        branchpoints = np.vstack(self.df_paths.connect_to_at)
        branchpoints = unique_row(branchpoints)
        num_branchpoints = len(branchpoints)

        max_branch_order = max(self.df_paths.corder)
        max_strahler_order = max(self.df_paths.sorder)

        terminalpaths = self.df_paths.path[self.df_paths.connected_by.apply(len) == 0].as_matrix()
        terminalpoints = np.vstack([p[-1] for p in terminalpaths])
        num_terminalpoints = len(terminalpoints)

        outerterminals = get_outer_terminals(terminalpoints)

        num_irreducible_nodes = num_branchpoints + num_terminalpoints

        num_dendritic_segments = len(self.df_paths)

        # path length

        dendritic = self.df_paths['dendritic_length']
        dendritic_sum = dendritic.sum()
        dendritic_mean = dendritic.mean()
        dendritic_median = dendritic.median()
        dendritic_min = dendritic.min()
        dendritic_max = dendritic.max()

        euclidean = self.df_paths['euclidean_length']
        euclidean_sum = euclidean.sum()
        euclidean_mean = euclidean.mean()
        euclidean_median = euclidean.median()
        euclidean_min = euclidean.min()
        euclidean_max = euclidean.max()

        tortuosity = dendritic / euclidean
        average_tortuosity = np.mean(tortuosity)      

        # node angles
        average_nodal_angle_deg, average_nodal_angle_rad, average_local_angle_deg, average_local_angle_rad = get_average_angles(self.df_paths)

        summary = {
            "general": {
                "number_of_dendritic_segments": int(num_dendritic_segments),
                "number_of_branch_points": int(num_branchpoints),
                "number_of_irreducible_nodes": int(num_irreducible_nodes),
                "max_branch_order": int(max_branch_order),
                "max_strahler_order": int(max_strahler_order),
            },
            "angle": {
                "average_nodal_angle_in_degree": average_nodal_angle_deg,
                "average_nodal_angle_in_radian": average_nodal_angle_rad,
                "average_local_angle_in_degree": average_local_angle_deg,
                "average_local_angle_in_radian": average_local_angle_rad,
            },
            "length": {
                "tortuosity": average_tortuosity,
                "dendritic": {
                    "sum": dendritic_sum,
                    "mean": dendritic_mean,
                    "median": dendritic_median,
                    "min": dendritic_min,
                    "max": dendritic_max,
                },
                "euclidean":{
                    "sum": euclidean_sum,
                    "mean": euclidean_mean,
                    "median": euclidean_median,
                    "min": euclidean_min,
                    "max": euclidean_max,
                }
            }

        }

        if self.linestack is not None:

            if self.voxelsize is None:
                voxelsize = [1,1,1]
            else:
                voxelsize = self.voxelsize

            dendritic_center = self.dendritic_center

            soma_coordinates = self.soma_on_stack

            asymmetry = np.sqrt(((soma_coordinates[:2] - dendritic_center)**2).sum())

            all_termianls_to_dendritic_center = np.sqrt(np.sum((terminalpoints[:,:2] - dendritic_center) ** 2,  1)) * voxelsize[0]
            out_terminals_to_dendritic_center = np.sqrt(np.sum((outerterminals[:, :2]- dendritic_center) **2, 1)) * voxelsize[0]

            typical_radius = np.mean(all_termianls_to_dendritic_center)
            outer_radius = np.mean(out_terminals_to_dendritic_center)

            import cv2
            dendritic_area = cv2.contourArea(outerterminals[:, :2].astype(np.float32)) * voxelsize[0]**2/1000

            summary.update({"density": {
                "asymmetry": asymmetry,
                "outer_radius": outer_radius,
                "typical_radius": typical_radius,
                "dendritic_area": dendritic_area
            }})

        if print_results:
            print_summary(summary, self.unit_df)

        if save_to is not None:
            import json
            with open(save_to, 'w') as f:
                json.dump(summary, f)

        # return summary


    def show_threeviews(self, order='c'):

        import matplotlib.pyplot as plt
        from matplotlib_scalebar.scalebar import ScaleBar

        plt.figure(figsize=(16,16))
        ax1 = plt.subplot2grid((4,4), (0,1), rowspan=3, colspan=3)
        ax2 = plt.subplot2grid((4,4), (0,0), rowspan=3, colspan=1)
        ax3 = plt.subplot2grid((4,4), (3,1), rowspan=1, colspan=3)
        ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1, colspan=1)  

        df_paths = self.df_paths
        soma = self.soma

        # xylims, zlims = find_lims(df_paths)
        # lims = (xylims, zlims)
        # maxlims = (np.max(np.vstack(df_paths.path), 0)[1:]).astype(int) 
        # maxlims = np.hstack([maxlims[0], maxlims]) + 30
        # minlims = (np.min(np.vstack(df_paths.path), 0)[1:]).astype(int)

        lims = find_lims(df_paths)

        plot_skeleten(ax2, df_paths, soma, 2, 0, order, lims)
        plot_skeleten(ax3, df_paths, soma, 1, 2, order, lims)
        plot_skeleten(ax1, df_paths, soma, 1, 0, order, lims)
        scalebar = ScaleBar(1, units=self.unit_df, location='lower left', box_alpha=0)
        ax1.add_artist(scalebar)    
        ax4.axis('off')

        plt.show()
    
    def show_density(self):

        try:
            density_stack = self.density_stack
            voxelsize = self.voxelsize
        except:
            logging.info('No density stack. Please provide the voxel sizes of the `.swc` file.')
            return None
        
        import matplotlib.pyplot as plt
        from matplotlib_scalebar.scalebar import ScaleBar

        linestack = self.linestack
        dendritic_center = self.dendritic_center
        soma_on_stack = self.soma_on_stack
        
        plt.figure(figsize=(16, 16))
        plt.imshow(density_stack.sum(2), cmap=plt.cm.gnuplot2_r, origin='lower')
        plt.scatter(dendritic_center[1], dendritic_center[0], color='g', marker='*', s=180, label='Dendritic Center')
        plt.scatter(soma_on_stack[1], soma_on_stack[0], color='r',  marker='*', s=180, label='Soma')
        
        linestack_xy = linestack.sum(2)
        linestack_xy[linestack_xy !=0] = 1
        linestack_xy = np.ma.masked_array(linestack_xy, ~linestack.any(2))
        plt.imshow(linestack_xy, origin='lower', cmap=plt.cm.binary)
        
        plt.legend(frameon=False)

        scalebar = ScaleBar(voxelsize[0], units=self.unit_df, location='lower left', box_alpha=0)
        plt.gca().add_artist(scalebar)   
        
        plt.axis('off')