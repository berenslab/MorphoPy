# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from unittest import mock

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
MOCK_MODULES = ['numpy', 'scipy', 'scipy.io', 'scipy.ndimage', 'scipy.ndimage.filters', 'scipy.sparse', 'scipy.spatial', 'scipy.interpolate', 'scipy.signal', 'matplotlib', 'matplotlib.colors', 'matplotlib.pyplot', 'matplotlib.artist', 'matplotlib.font_manager', 'matplotlib.rcsetup', 'matplotlib.patches', 'matplotlib_scalebar.scalebar', 'matplotlib.offsetbox', 'pandas', 'seaborn', 'sklearn', 'sklearn.decomposition']
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()

sys.modules['scipy.stats'] = mock.Mock(rv_continuous=object)

sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'Morphopy'
copyright = '2020, Sophie Laturnus, Adam von Daranyi, Ziwei Huang and Philipp Berens'
author = 'Sophie Laturnus, Adam von Daranyi, Ziwei Huang and Philipp Berens'

# The full version, including alpha/beta/rc tags
release = '0.7'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
    html_theme_options = {
        "collapse_navigation": False
    }

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

htmlhelp_basename = 'MorphoPydoc'

intersphinx_mapping = {'python': ('https://docs.python.org/3', None)}