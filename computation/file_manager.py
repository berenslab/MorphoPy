import pandas as pd
import neurontree.NeuronTree as nt
import neurontree.utils as utils


def load_swc_file(filename=None, nxversion=1):
    """
    This function imports a swc file to a pandas dataframe and then uses utils to standardize it.
    This dataframe is then used to initialize a new NeuronTree object and returns it to the caller.
    :param filename: path to swc file as string
    :param nxversion: version of installed networkX passed to NeuronTree
    :return: NeuronTree object
    """
    f = lambda x: float(x.replace(",", "."))
    swc = pd.read_csv(filename, delim_whitespace=True, comment='#', converters={'x': f, 'y': f, 'z': f, 'radius': f},
                      names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)

    # check returned object if import was successful
    if swc is None or type(swc) != pd.DataFrame or swc.size == 0:
        raise ValueError('No points found in swc file, please check format of data!')

    swc = utils.get_standardized_swc(swc)
    my_tree = nt.NeuronTree(swc=swc, nxversion=nxversion)
    check_neurontree(my_tree)
    return my_tree


def check_neurontree(neurontree=None):
    """
    This function checks a neurontree
    :param filename: path to swc file as string
    :param nxversion: version of installed networkX passed to NeuronTree
    :return: NeuronTree object
    """
    if not neurontree.has_root():
        raise ValueError('Graph has no root defined!')
    if not neurontree.is_connected():
        raise ValueError('Graph has disconnected nodes!')
