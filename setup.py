from setuptools import setup, find_packages

setup(name='morphopy',
      version='0.1.0',
      description='Analyze morphological data from ImageJ/Simple Neurite Tracer',
      author='Ziwei Huang',
      author_email='huang-ziwei@outlook.com',
      url='https://github.com/berenslab/MorphoPy',
      packages=find_packages(),
      install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "matplotlib-scalebar",
      ],
     )
