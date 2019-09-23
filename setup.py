from setuptools import setup, find_packages

setup(name='morphopy',
      version='0.2.0',
      description='Analyze morphological data from ImageJ/Simple Neurite Tracer',
      author='Ziwei Huang',
      author_email='huang-ziwei@outlook.com',
      url='https://github.com/berenslab/MorphoPy',
      packages=find_packages(),
      install_requires=[
        "numpy>=1.13.1",
        "pandas==0.20.3",
        "scipy",
        "matplotlib",
        "matplotlib-scalebar",
        "networkx<2.0"
      ],
     )
