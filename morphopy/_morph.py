import os
import pickle
import logging

from ._utils import *

__all__ = ['Morph']

class Morph(object):

    def __init__(self, data, voxelsize=None, imagesize=None, threshold=20, loglevel='INFO'):

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

        df_swc = read_swc(data)
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
        self.soma = soma


    def topviews(self, highlight=[], order_type='c'):
        
        import matplotlib.pyplot as plt
        # %matplotlib inline    
        
        df_paths = self.df_paths
        soma = self.soma

        if order_type == 'c':
            colors = plt.cm.viridis.colors
            colors_idx = np.linspace(0, 255, max(df_paths.corder)+1).astype(int)
        elif order_type == 's':
            colors = plt.cm.viridis_r.colors
            colors_idx = np.linspace(0, 255, max(df_paths.sorder)+1).astype(int)
        
        plt.figure(figsize=(12,12))
        plt.scatter(soma[0][1], soma[0][0], s=280, color='grey')
        for row in df_paths.iterrows():
            
            path_id = row[0]
            
            path = row[1]['path']
            bpt = path[0]
            if order_type == 'c':
                order = row[1]['corder']
                plt.plot(path[:, 1], path[:, 0], color=colors[colors_idx[order]])
                plt.scatter(bpt[1], bpt[0], color=colors[colors_idx[order]], zorder=1)

            elif order_type == 's':
                order = row[1]['sorder']
                plt.plot(path[:, 1], path[:, 0], color=colors[colors_idx[order-1]])
                plt.scatter(bpt[1], bpt[0], color=colors[colors_idx[order-1]], zorder=1)

            if path_id in highlight:
                plt.plot(path[:, 1], path[:, 0], color='red')


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