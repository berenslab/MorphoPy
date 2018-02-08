import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

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


def plot_persistence_diagram(data, ax):

    ax.scatter(data['birth'], data['death'], s=4, c='k', alpha=.4)
    ax.plot([0, np.max(data['birth'])], [0, np.max(data['death'])])
    sns.despine()
    ax.set_xlabel('birth [dist from soma in um]')
    ax.set_ylabel('death [dist from soma in um]')
    ax.set_title('persistence diagram')

    return ax


def plot_persistence_image_1d(data, ax):
    """
    Plots the persistence as a 1 dimensional persistence image as defined in _Metrics for comparing neuronal tree
    shapes based on persistent homology_ Y. Li, D. Wang, G. Ascoli et al. , 2017.

    The 1D persistence image is defined as a sum of Gaussian kernels located at the time of birth of each branch and
    weighted by its lifetime (|birth - death|).
    Formally:
        $p_D(x) = \sum_{i=1}^{k} |y_i - x_i|\cdot K_t(x,x_i)$
        where $K_t(x,x_i)$ denotes a Gaussian kernel centered at $x_i$ with width $t$. Here $t$ is chosen to be $50$ as
        in the original paper.
    :param data: pandas.DataFrame holding the persistence data.
    :return: figure and axis of the plot
    """
    steps = 100
    y = np.zeros((steps,))
    x = np.linspace(0, np.max(data['birth']), steps)
    t = 50

    for k, p in data.iterrows():
        m = np.abs(p['birth'] - p['death'])
        y += m * norm.pdf(x, loc=p['birth'], scale=t)

    ax.plot(x, y)
    ax.set_xlabel('birth [dist from soma in um]')
    ax.set_ylabel('persistence p_D(x)')
    sns.despine()
    ax.set_title('persistence image 1D \n (Wang et al., 2017)')

    return ax


def plot_persistence_image_2d(data, ax):
    sns.kdeplot(data['birth'], data['death'], ax=ax)
    ax.scatter(data['birth'], data['death'], s=4, c='k', alpha=.4)
    ax.set_xlabel('birth [dist from soma in um]')
    ax.set_ylabel('death [dist from soma in um]')
    sns.despine()
    ax.set_title('persistence image 2D \n(Kanari et al., 2016) \n https://arxiv.org/abs/1603.08432')
    return ax
