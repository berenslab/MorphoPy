import pandas as pd
import configparser as cp
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

def read_config(configfile=None):
    """
    This function reads a config file and returns two dictionaries with all found values
    :param configfile: path to the configfile
    :return: (config_params) dictionary with global params and normalization bounds with key, value pairs
                            None will be returned if no file or wrong file was found
    """
    try:
        # read from configfile if available
        if configfile is not None:
            cfg = cp.ConfigParser()
            cfg.read(configfile)

            config_params = {}
            if cfg.has_option("global", "proj_axes"):
                config_params['proj_axes'] = cfg.get("global", "proj_axes")
            if cfg.has_option("global", "n_bins"):
                config_params['n_bins'] = cfg.getint("global", "n_bins")
            if cfg.has_option("global", "normed"):
                config_params['normed'] = cfg.getboolean("global", "normed")
            if cfg.has_option("global", "smooth"):
                config_params['smooth'] = cfg.getboolean("global", "smooth")
            if cfg.has_option("global", "sigma"):
                config_params['sigma'] = cfg.getint("global", "sigma")

            if cfg.has_section("norm_bound"):
                config_params["r_max_x"] = cfg.getfloat("norm_bound", "r_max_x", fallback=0)
                config_params["r_max_y"] = cfg.getfloat("norm_bound", "r_max_y", fallback=0)
                config_params["r_max_z"] = cfg.getfloat("norm_bound", "r_max_z", fallback=0)
                config_params["r_min_x"] = cfg.getfloat("norm_bound", "r_min_x", fallback=0)
                config_params["r_min_y"] = cfg.getfloat("norm_bound", "r_min_y", fallback=0)
                config_params["r_min_z"] = cfg.getfloat("norm_bound", "r_min_z", fallback=0)

            return config_params
        else:
            return None
    except FileNotFoundError:
        print(f'Failure in computing density map: Config file not found or not readable: {configfile}')
        return None


