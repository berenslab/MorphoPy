import logging

def check_soma(types):

    """
    Check if soma exist.
    """

    logging.info('  Has soma (Type 1)?')
    if 1 in types:
        logging.info('\tYes\n')
        return True
    else:
        logging.info('\tNo. The first point in .swc file is used as soma.')
        logging.info('\tThis is not always accurate, please check the file.\n')
        return False

def check_axon(types):

    """
    Check if axon exist.
    """

    logging.info('  Has axon (Type 2)?')
    if 2 in types:
        logging.info('\tYes\n')
        return True
    else:
        logging.info('\tNo\n')
        return False

def check_basal_dendrites(types):

    """
    Check if (basal) dendrites exist.
    """

    logging.info('  Has basal dendrites (Type 3)?')
    if 3 in types:
        logging.info('\tYes\n')
        return True
    else:
        logging.info('\tNo\n')
        return False

def check_apical_dendrites(types):

    """
    Check if apical dendrites exist.
    """

    logging.info('  Has apical dendrites (Type 4)?')
    if 4 in types:
        logging.info('\tYes\n')
        return True
    else:
        logging.info('\tNo\n')
        return False

def check_others(types):
    """
    Check if undefined or custom types exist.
    """
    logging.info('  Has undefined or custom types (Type 0/5)?')
    if 0 in types or 5 in types:
        logging.info('\tYes\n')
        return True
    else:
        logging.info('\tNo\n')
        return False


def check_swc(df_swc):

    """
    Check if if the input swc valid or not.
    """

    types = df_swc.type.unique()

    hassoma = check_soma(types)
    hasaxon = check_axon(types)
    hasbasaldendrites = check_basal_dendrites(types)
    hasapicaldenderites = check_apical_dendrites(types)
    hasothers = check_others(types)
