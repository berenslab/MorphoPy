### TEST READING METHODS ####

import networkx as nx
from morphopy.computation.file_manager import load_swc_file


def test_read_swc_returned_fileformat():

    filepath = '../data/Image001-005-01.CNG.swc'
    G = load_swc_file(filepath)

    assert(isinstance(G.get_graph(), nx.DiGraph)), "read_swc() should return a graph as networkx.DiGraph"

