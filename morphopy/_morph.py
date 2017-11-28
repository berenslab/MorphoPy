import os
import pickle
import logging

from ._utils import *

__all__ = ['Morph']

class Morph(object):

    def __init__(self, data, unit='um', voxelsize=None, imagesize=None, threshold=30, loglevel='INFO'):

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


        self.voxelsize = voxelsize
        self.imagesize = imagesize
        self.unit_swc = unit
        # df_paths | real length

        df_swc, unit_df = read_swc(data, unit, voxelsize)
        soma, neurites = get_soma(df_swc)
        df_paths = get_df_paths(df_swc) 
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

        self.linestack = swc2linestack(data, self.unit_swc, imagesize, voxelsize)

        if self.linestack is not None:
            self.densitystack, self.dendritic_center = calculate_density(self.linestack, voxelsize)


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


    def cell_statistics(self):

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

        # Summary
        logging.info('  Statistics of the cell')
        logging.info('  ======================\n')

        logging.info('  # Meta Infomation\n')
        logging.info('    Number of dendritic arbor segment: {}'.format(num_dendritic_segments))
        logging.info('    Number of branch points: {}'.format(num_branchpoints))
        logging.info('    Number of irreducible nodes: {}\n'.format(num_irreducible_nodes))

        logging.info('    Max branching order: {}'.format(max_branch_order))
        logging.info('    Max Strahler order: {}\n\n'.format(max_strahler_order))
        
        logging.info('  # Angle \n')
        logging.info('    Average nodal angle in degree: {:.3f}'.format(average_nodal_angle_deg))
        logging.info('    Average nodal angle in radian: {:.3f} \n'.format(average_nodal_angle_rad))
        logging.info('    Average local angle in degree: {:.3f}'.format(average_local_angle_deg))
        logging.info('    Average local angle in radian: {:.3f} \n'.format(average_local_angle_rad))
        
        logging.info('  ## Average tortuosity: {:.3f}\n'.format(average_tortuosity))
        
        logging.info('  ## Dendritic length ({})\n'.format(self.unit_df))
        # logging.info('  ## Dendritic length\n')
        logging.info('     Sum: {:.3f}'.format(dendritic_sum))
        logging.info('     Mean: {:.3f}'.format(dendritic_mean))
        logging.info('     Median: {:.3f}'.format(dendritic_median))
        logging.info('     Min: {:.3f}'.format(dendritic_min))
        logging.info('     Max: {:.3f}\n'.format(dendritic_max))
        
        logging.info('  ## Euclidean length ({})\n'.format(self.unit_df))
        # logging.info('  ## Euclidean length\n')

        logging.info('     Sum: {:.3f}'.format(euclidean_sum))
        logging.info('     Mean: {:.3f}'.format(euclidean_mean))
        logging.info('     Median: {:.3f}'.format(euclidean_median))
        logging.info('     Min: {:.3f}'.format(euclidean_min))
        logging.info('     Max: {:.3f}\n'.format(euclidean_max))

        # density related 

        if self.linestack is not None:

            if self.voxelsize is None:
                voxelsize = [1,1,1]
            else:
                voxelsize = self.voxelsize

            dendritic_center = self.dendritic_center

            if self.unit_swc == 'um':
                soma_coordinates = self.soma / self.voxelsize
            else:
                soma_coordinates = self.soma

            asymmetry = np.sqrt(((soma_coordinates[0][:2] - dendritic_center)**2).sum())

            all_termianls_to_dendritic_center = np.sqrt(np.sum((terminalpoints[:,:2] - dendritic_center) ** 2,  1)) * voxelsize[0]
            out_terminals_to_dendritic_center = np.sqrt(np.sum((outerterminals[:, :2]- dendritic_center) **2, 1)) * voxelsize[0]

            typical_radius = np.mean(all_termianls_to_dendritic_center)
            outer_radius = np.mean(out_terminals_to_dendritic_center)

            import cv2
            dendritic_area = cv2.contourArea(outerterminals[:, :2].astype(np.float32)) * voxelsize[0]**2/1000

            logging.info('  # Density related ({})\n'.format(self.unit_df))
            logging.info('    Asymmetry: {:.3f}'.format(asymmetry))
            logging.info('    Outer Radius: {:.3f}'.format(outer_radius))
            logging.info('    Typical Radius : {:.3f}\n'.format(typical_radius))
            logging.info('    Dendritic Area: {:.3f} ×10\u00b3 um\u00b2\n\n'.format(dendritic_area))

    def threeviews(self, order_type='c'):

        import matplotlib.pyplot as plt
        from matplotlib_scalebar.scalebar import ScaleBar
        
        plt.figure(figsize=(16,16))
        ax1 = plt.subplot2grid((4,4), (0,1), rowspan=3, colspan=3)
        ax2 = plt.subplot2grid((4,4), (0,0), rowspan=3, colspan=1)
        ax3 = plt.subplot2grid((4,4), (3,1), rowspan=1, colspan=3)
        ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1, colspan=1)  

        df_paths = self.df_paths
        soma = self.soma

        xylims, zlims = find_lims(df_paths)

        # maxlims = (np.max(np.vstack(df_paths.path), 0)[1:]).astype(int) 
        # maxlims = np.hstack([maxlims[0], maxlims]) + 30
        # minlims = (np.min(np.vstack(df_paths.path), 0)[1:]).astype(int)

        plot_skeleten(ax2, df_paths, soma, 2, 0, order_type, lims)
        plot_skeleten(ax3, df_paths, soma, 1, 2, order_type, lims)
        plot_skeleten(ax1, df_paths, soma, 1, 0, order_type, lims)
        scalebar = ScaleBar(1, units=self.unit_df, location='lower left', box_alpha=0)
        ax1.add_artist(scalebar)    
        ax4.axis('off')
    