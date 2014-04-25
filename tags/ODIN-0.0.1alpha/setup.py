import os
import sys
from setuptools import setup, find_packages

#################################################
# Data Config
#################################################

# Paths
data_dir = os.path.join(os.getenv('HOME'), "rgtdata")
for i in range(0,len(sys.argv)):
    if("--rgt-data-path" in sys.argv[i]):
        data_dir = sys.argv[i].split("=")[1]
        if("~/" in data_dir): data_dir = os.path.join(os.getenv('HOME'), data_dir[2:])
        elif("~" in data_dir): data_dir = os.getenv('HOME')
        sys.argv.remove(sys.argv[i])
script_dir = os.path.dirname(os.path.abspath(__file__))
os.system("mkdir -p "+data_dir)

# Creating data.config
data_config_file_name = os.path.join(data_dir, "data.config")
data_config_file = open(data_config_file_name,"w")
data_config_file.write("[GenomeData]\n")
data_config_file.write("organism: hg19\n")
data_config_file.write("genome: genome.fa\n")
data_config_file.write("chromosome_sizes: chrom.sizes\n")
data_config_file.write("association_file: association_file.bed\n\n")
data_config_file.write("[MotifData]\n")
data_config_file.write("pwm_dataset: PWM\n")
data_config_file.write("logo_dataset: logo\n\n")
data_config_file.close()

# Creating data.config.path
data_config_path_file_name = os.path.join(script_dir,"rgt","data.config.path")
data_config_path_file = open(data_config_path_file_name,"w")
data_config_path_file.write(data_config_file_name)
data_config_path_file.close()

# Copying data from package folder to installation folder
folderNameList = ["hg19","mm9"]
fileNameList = ["association_file.bed","chrom.sizes"]
for e1 in folderNameList:
    os.system("mkdir -p "+os.path.join(data_dir, e1))
    for e2 in fileNameList:
        os.system("cp "+os.path.join(script_dir,"data",e1,e2)+" "+os.path.join(data_dir, e1))

#################################################
# EG Data config - deprecated
#################################################
packagePathFile = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rgt/packagePathFile.txt"),"w")
packagePathFile.write(os.path.dirname(os.path.abspath(__file__))+"/")
packagePathFile.close()

#################################################
# Setup
#################################################

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
      package_data={'rgt': ['packagePathFile.txt','data.config.path']},
      entry_points = {
          'console_scripts': [
              'rgt-motifanalysis = rgt.motifanalysis.main:main', 'rgt-ODIN = rgt.ODIN.ODIN:main'
          ]
      }
      #install_requires=['numpy>=1.4.0','scipy>=0.7.0','Biopython>=1.60','pandas==0.7.1','fisher>=0.1.0','statsmodels>=0.4.0','HTML>=0.04','matplotlib>=1.1.0'] # TODO Put the correct versions
      )

# Removing data.config.path
os.system("rm "+data_config_path_file_name)

# Modifying permissions
currUser = os.getenv('HOME').split("/")[-1]
os.system("chown -R "+currUser+" "+os.path.join(data_dir,"."))
os.system("chgrp -R "+currUser+" "+os.path.join(data_dir,"."))
os.system("chmod -R 770 "+data_dir)


