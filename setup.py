import os
from setuptools import setup, find_packages
packagePathFile = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rgt/packagePathFile.txt"),"w")
packagePathFile.write(os.path.dirname(os.path.abspath(__file__))+"/")
packagePathFile.close()
version = "0.0.1"
README = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.txt")
long_description = open(README).read() + "nn"
setup(name="RGT",
      version=version,
      description=("Toolkit to perform common analysis of regulatory genomics data"),
      long_description=long_description,
      classifiers=[
          'Topic :: Scientific/Engineering :: Bio-Informatic',
          'Topic :: Scientific/Engineering :: Artificial Intelligence'
          ],
      keywords='motif',
      author='Ivan G. Costa, Eduardo G. Gusmao, Manuel Allhoff',
      author_email='software@costalab.org',
      license='GPL',
      packages=find_packages(),
      package_data={'rgt': ['packagePathFile.txt']},
      entry_points = {
          'console_scripts': [
              'rgt-motifanalysis = rgt.motifanalysis.main:main'
          ]
      },
      install_requires=['numpy>=1.4.0','scipy>=0.7.0','Biopython>=1.60','pandas==0.7.1','fisher>=0.1.0','statsmodels>=0.4.0','HTML>=0.04','matplotlib>=1.1.0'] # TODO Put the correct versions
      )
