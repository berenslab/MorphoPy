import numpy as np
import matplotlib.pyplot as plt

def find_lims(df_paths):
    """

    :param df_paths:
    :return:
    """
    points = np.vstack(df_paths.path)
    maxlims = np.max(points, 0)
    minlims = np.min(points, 0)
    xylims = np.hstack([maxlims[:2], minlims[:2]])
    zlims = np.hstack([maxlims[2], minlims[2]])
    if (xylims >= 0).all():
        xylims = np.array([0, max(xylims).astype(int) + 30])
    else:
        xylims = np.array([-max(abs(xylims)).astype(int) - 30, max(abs(xylims)).astype(int) + 30])

    if (zlims >= 0).all():
        zlims = np.array([-10, max(zlims).astype(int) + 30])
    else:
        zlims = np.array([-max(abs(zlims)).astype(int) - 30, max(abs(zlims)).astype(int) + 30])
        
    return xylims, zlims

def plot_skeleton(ax, df_paths, soma, axis0, axis1, order_type, lims):
    """

    :param ax:
    :param df_paths:
    :param soma:
    :param axis0:
    :param axis1:
    :param order_type:
    :param lims:
    :return:
    """

    if order_type == 'c':
        colors = plt.cm.viridis.colors
        colors_idx = np.linspace(0, 255, max(df_paths.corder)+1).astype(int)
    elif order_type == 's':
        colors = plt.cm.viridis_r.colors
        colors_idx = np.linspace(0, 255, max(df_paths.sorder)+1).astype(int)
        
    ax.scatter(soma[axis0], soma[axis1], s=280, color='grey')
    for row in df_paths.iterrows():

        path_id = row[0]

        path = row[1]['path']
        bpt = path[0]
        if order_type == 'c':
            order = row[1]['corder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[int(order)]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[int(order)]], zorder=1)

        elif order_type == 's':
            order = row[1]['sorder']
            ax.plot(path[:, axis0], path[:, axis1], color=colors[colors_idx[int(order)-1]])
            ax.scatter(bpt[axis0], bpt[axis1], color=colors[colors_idx[int(order)-1]], zorder=1)
    
    xylims, zlim = lims

    if axis0 == 2 and axis1 == 0: # ax2
        ax.set_xlim(zlim[0], zlim[1])
        ax.set_ylim(xylims[0], xylims[1])   

    elif axis0 == 1 and axis1 == 2: # ax3
        ax.set_xlim(xylims[0], xylims[1])
        ax.set_ylim(zlim[0], zlim[1])
     

    elif axis0 == 1 and axis1 == 0: # ax1
        ax.set_xlim(xylims[0], xylims[1])
        ax.set_ylim(xylims[0], xylims[1])

    ax.axis('off')