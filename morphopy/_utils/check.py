import logging

def check_soma(types):

    """
    Check if soma exist.
    """

    logging.info('  Has soma points?')
    if 1 in types:
        logging.info('\t Yes')
        return True
    else:
        logging.info('\t No')
        return False
    
def check_axon(types):
    
    """
    Check if axon exist.
    """    
    
    logging.info('  Has axon?')
    if 2 in types:
        logging.info('\t Yes')
        return True
    else:
        logging.info('\t No')
        return False

def check_basal_dendrites(types):
    
    """
    Check if (basal) dendrites exist.
    """    
    
    logging.info('  Has basal dendrites?')
    if 3 in types:
        logging.info('\t Yes')
        return True
    else:
        logging.info('\t No')
        return False
    
def check_apical_dendrites(types):
    
    """
    Check if apical denderites exist.
    """    
    
    logging.info('  Has apical dendrites?')
    if 4 in types:
        logging.info('\t Yes')
        return True
    else:
        logging.info('\t No')
        return False

def check_swc(df_swc):

    """
    Check if if the input swc valid or not.
    """

    types = df_swc.type.unique()

    check_soma(types)
    check_axon(types)
    check_basal_dendrites(types)
    check_apical_dendrites(types)    