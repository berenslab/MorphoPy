import numpy as np
import matplotlib.pyplot as plt

def plot_morph(ax, df_paths, view, plot_axon, plot_dendrites):
    
    if view == 'xy':
        axis0 = 0
        axis1 = 1
    elif view == 'xz':
        axis0 = 0
        axis1 = 2
    elif view == 'yz':
        axis0 = 1
        axis1 = 2
    
    
    soma = df_paths[df_paths.type == 1].path[0][0]
    axon = df_paths[df_paths.type == 2]
    dendrites = df_paths[df_paths.type == 3]
        

    ax.scatter(0, 0, s=280, color='grey')
    
    if plot_dendrites:
        
        dcolors_idx = np.linspace(0, 255, max(dendrites.branch_order)+1).astype(int)
        dcolors = np.vstack(plt.cm.Reds_r(dcolors_idx))[:, :3]

        for row in dendrites.iterrows():

            path_id = row[0]
            path = row[1]['path'] - soma
            order = row[1]['branch_order']
            bpt = path[0]     

            dend_plot = ax.plot(path[:, axis0], path[:, axis1], color=dcolors[int(order)])
            ax.scatter(bpt[axis0], bpt[axis1], color=dcolors[int(order)], zorder=1)


    
    if plot_axon:
        
        acolors_idx = np.linspace(0, 255, max(axon.branch_order)+1).astype(int)
        acolors = np.vstack(plt.cm.Blues_r(acolors_idx))[:, :3]
        
        for row in axon.iterrows():

            path_id = row[0]
            path = row[1]['path'] - soma
            order = row[1]['branch_order']
            bpt = path[0]     

            axon_plot = ax.plot(path[:, axis0], path[:, axis1], color=acolors[int(order)])
            ax.scatter(bpt[axis0], bpt[axis1], color=acolors[int(order)], zorder=1)

    lim_max = int(np.ceil((np.vstack(df_paths.path.as_matrix()) - soma).max() / 20) * 20)
#     lim_min = int(np.floor((np.vstack(self.df_paths.path.as_matrix()) - soma).min() / 20) * 20)

    ax.set_xlim(-lim_max, lim_max)
    ax.set_ylim(-lim_max, lim_max)
    
    ax.set_title('{}'.format(view))
    
    ax.axis('off')

    return ax

# def find_lims(df_paths):
#     """

#     :param df_paths:
#     :return:
#     """
#     points = np.vstack(df_paths.path)
#     maxlims = np.max(points, 0)
#     minlims = np.min(points, 0)
#     xylims = np.hstack([maxlims[:2], minlims[:2]])
#     zlims = np.hstack([maxlims[2], minlims[2]])
#     if (xylims >= 0).all():
#         xylims = np.array([0, max(xylims).astype(int) + 30])
#     else:
#         xylims = np.array([-max(abs(xylims)).astype(int) - 30, max(abs(xylims)).astype(int) + 30])

#     if (zlims >= 0).all():
#         zlims = np.array([-10, max(zlims).astype(int) + 30])
#     else:
#         zlims = np.array([-max(abs(zlims)).astype(int) - 30, max(abs(zlims)).astype(int) + 30])
        
#     return xylims, zlims

# def plot_skeleton(ax, df_paths, soma, axis0, axis1, lims):
#     """

#     :param ax:
#     :param df_paths:
#     :param soma:
#     :param axis0:
#     :param axis1:
#     :param lims:
#     :return:
#     """

#     colors = plt.cm.viridis.colors
#     colors_idx = np.linspace(0, 255, max(df_paths.branch_order)+1).astype(int)
        
#     ax.scatter(soma[axis0], soma[axis1], s=280, color='grey')
#     for row in df_paths.iterrows():

#         path_id = row[0]

#         path = row[1]['path']
#         bpt = path[0]

#         order = row[1]['branch_order']
#         ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[int(order)]])
#         ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[int(order)]], zorder=1)
    
#     xylims, zlim = lims

#     if axis0 == 2 and axis1 == 0: # ax2
#         ax.set_xlim(zlim[0], zlim[1])
#         ax.set_ylim(xylims[0], xylims[1])   

#     elif axis0 == 1 and axis1 == 2: # ax3
#         ax.set_xlim(xylims[0], xylims[1])
#         ax.set_ylim(zlim[0], zlim[1])
     

#     elif axis0 == 1 and axis1 == 0: # ax1
#         ax.set_xlim(xylims[0], xylims[1])
#         ax.set_ylim(xylims[0], xylims[1])

#     ax.axis('off')