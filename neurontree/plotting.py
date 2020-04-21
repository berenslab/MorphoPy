from neuronTree import NeuronTree
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar


def show_threeview(nt, fig=None, **kwargs):
    """
    Plots all three two-dimensional projections of the neuron in nt.
    :param nt: NeuronTree, neuron to be plotted
    :param fig: (optional) pass the figure handle from outside if you want more control.
    :param kwargs: arguments that can be passed to the NeuronTree.draw_2D() function.
    """

    if not fig:
        fig = plt.figure(figsize=(16, 16))

    ax1 = plt.subplot2grid((4, 4), (0, 1), rowspan=3, colspan=3)
    ax2 = plt.subplot2grid((4, 4), (0, 0), rowspan=3, colspan=1)
    ax3 = plt.subplot2grid((4, 4), (3, 1), rowspan=1, colspan=3)
    ax4 = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=1)

    nt.draw_2D(fig, ax=ax2, projection='zy', **kwargs)
    nt.draw_2D(fig, ax=ax3, projection='xz', **kwargs)
    nt.draw_2D(fig, ax=ax1, projection='xy', **kwargs)

    scalebar = ScaleBar(1, units='um', location='lower left', box_alpha=0)
    ax1.add_artist(scalebar)
    ax4.axis('off')