# Code availability and Installation via Github

RGT code is deposited is github (https://github.com/ManuelAllh/reg-gen/). You can download and install RGT core classes as following.

```
git checkout https://github.com/ManuelAllh/reg-gen.git
cd reg-gen
sudo python setup.py install
```

Detailled installation instructions are found in:
http://www.regulatory-genomics.org/rgt/download-installation

# Download

Although you can download the developer version via github, we also provide RGT stable releases. Please find the latest stable distributions of RGT in the following link:

https://drive.google.com/file/d/0B77RY6Xty6pSa05CR0pMWVlhdGc/view

# Basic Installation via Command Line

The RGT will automatically install all python package requirements. However, we highly recommend that you pre-install Cython, numpy and scipy. You can do this with the following commands:

```
pip install cython
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
```

To install the rgt package change to the RGT code folder and simply type:

```
sudo python setup.py install
```

# Installation Requirements

If you experience any problems with the installation regarding the package requirements, you can try to install these packages individually (please visit their website listed bellow for more information). The package requirements for the most current version is:

```
Python >= 2.7
Numpy >= 1.4.0
Scipy >= 0.7.0
Pysam >= 0.7.5
```

# Customized Installation

## Prefix scheme (installation without sudo)

In the introduction it was explained how to install the package in a fast way if you have superuser privileges. However, if you don’t have such privileges (for instance, you are running the tool in a server) you can install the package manually using the prefix scheme. This can be easily performed (for python version X.Y) by typing the following command:

```
export PATH=$PATH:~/app/bin
export PYTHONPATH=$PYTHONPATH:~/app/lib/pythonX.Y/site-packages
mkdir -p ~/app/bin ~/app/lib/pythonX.Y/site-packages
python setup.py install --prefix=~/app
```

The above command will install the package to your ~/app folder. You can change the location of this folder to whichever you want. Note that in the above commands we made some modifications to your PYTHONPATH and PATH variables. The PYTHONPATH variable should contain all the paths to python libraries and the PATH variable should contain the paths to all binaries you want to execute.

## Further installation options

Furthermore, we provide some extra options for RGT packages installation:

1. Data Path: RGT uses data such as the genome, genome size, default hidden Markov models, etc. All these data will be installed by default to ~/rgtdata. However you may want to modify such path by using “--rgt-data-path” option:

```
python setup.py install --prefix=~/app --rgt-data-path=~/app
```

Some RGT functionalities require data sets that are too big to be included in the RGT distribution package. However they can be easily obtained by running a few scripts available in the RGT distribution package. Please refer to

http://www.regulatory-genomics.org/rgt/rgt-data-folder/

for more information on how to download all data sets used by some RGT functionalities as well as to learn how to use your own data sets within the RGT data structure.

2. RGT Tool: . In order to install only a specific tool (e.g. HINT) please use the option “--rgt-tool” as in the example below:

```
python setup.py install --prefix=~/app --rgt-tool=hint
```

You can install more than one tool at a time by simply using a comma-separated list:

```
python setup.py install --prefix=~/app --rgt-tool=hint,motifanalysis,ODIN
```

If this option is not used, all RGT tools will be installed.

### Binaries

Some RGT functionalities use scripts from UCSC utilities website (linux x86). The required scripts are going to be automatically copied to the default installation location binaries (“bin”) folder. Therefore, you must change your $PATH variable (UNIX systems) to include the default installation binary folder. If you used the “prefix scheme” installation your “RGT INSTALLATION BIN FOLDER” will be “<PREFIX_INSTALLATION_FOLDER>/rgt-x.y.z/bin”, where x.y.z is the RGT distribution you have downloaded.  If you used the simple sudo installation scheme then please type the following command:

```
python -c "import site; print([e.replace(\"dist\",\"site\") for e in site.getsitepackages()]+site.getsitepackages());"
```

Search in all locations retrieved with the command above for the following folder: “rgt-x.y.z/bin”, where x.y.z is the RGT distribution you have downloaded. This folder is your “RGT INSTALLATION BIN FOLDER”. To add it to your $PATH variable please execute:

```
export $PATH=$PATH:<RGT_INSTALLATION_BIN_FOLDER>
```

