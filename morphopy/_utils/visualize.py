import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_morph(ax, df_paths, view, plot_axon, plot_basal_dendrites, plot_apical_dendrites):

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
    basal_dendrites = df_paths[df_paths.type == 3]
    apical_dendrites = df_paths[df_paths.type == 4]

    ax.scatter(0, 0, s=280, color='grey')

    if plot_basal_dendrites and len(basal_dendrites)>0:

        bdcolors_idx = np.linspace(0, 200, max(basal_dendrites.branch_order)+1).astype(int)
        bdcolors = np.vstack(plt.cm.Reds_r(bdcolors_idx))[:, :3]

        for row in basal_dendrites.iterrows():

            path_id = row[0]
            path = row[1]['path'] - soma
            order = row[1]['branch_order']
            bpt = path[0]

            dend_plot = ax.plot(path[:, axis0], path[:, axis1], color=bdcolors[int(order)])
            ax.scatter(bpt[axis0], bpt[axis1], color=bdcolors[int(order)], zorder=1)

    if plot_apical_dendrites and len(apical_dendrites)>0:

        adcolors_idx = np.linspace(0, 200, max(apical_dendrites.branch_order)+1).astype(int)
        adcolors = np.vstack(plt.cm.Purples_r(adcolors_idx))[:, :3]

        for row in apical_dendrites.iterrows():

            path_id = row[0]
            path = row[1]['path'] - soma
            order = row[1]['branch_order']
            bpt = path[0]

            dend_plot = ax.plot(path[:, axis0], path[:, axis1], color=adcolors[int(order)])
            ax.scatter(bpt[axis0], bpt[axis1], color=adcolors[int(order)], zorder=1)

    if plot_axon and len(axon)>0:

        acolors_idx = np.linspace(0, 200, max(axon.branch_order)+1).astype(int)
        acolors = np.vstack(plt.cm.Blues_r(acolors_idx))[:, :3]

        for row in axon.iterrows():

            path_id = row[0]
            path = row[1]['path'] - soma
            order = row[1]['branch_order']
            bpt = path[0]

            axon_plot = ax.plot(path[:, axis0], path[:, axis1], color=acolors[int(order)])
            ax.scatter(bpt[axis0], bpt[axis1], color=acolors[int(order)], zorder=1)

    lim_max = int(np.ceil((np.vstack(df_paths.path.as_matrix()) - soma).max() / 20) * 20)
    lim_min = int(np.floor((np.vstack(df_paths.path.as_matrix()) - soma).min() / 20) * 20)

    lim = max(abs(lim_max), abs(lim_min))

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax.set_title('{}'.format(view))

    ax.axis('off')

    return ax


def plot_persistence_diagram(data):

    fig = plt.figure()
    ax = fig.gca()
    ax.scatter(data['birth'], data['death'], s=4, c='k', alpha=.4)
    ax.plot([0, np.max(data['birth'])], [0, np.max(data['death'])])
    sns.despine()
    ax.xlabel('birth [dist from soma in um]')
    ax.ylabel('death [dist from soma in um]')

    return fig, ax