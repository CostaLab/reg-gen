# Installation

These instructions are only to be attempted if you know what you are doing and the basic instructions failed. Please contact us on our Support group if you experience any problems. Please notice that RGT only supports Python 3 now.

## Linux

The recommended way of installing RGT is via [pip](https://pip.pypa.io/en/stable/), the main python package manager. After installing pip, type in this order:

```shell
pip install cython numpy scipy
pip install RGT
```

This will install the full RGT suite with all python dependencies. A few more libraries will be needed, but should be installed via your own distribution’s package manager.

### Fail in installation?

Please try the following code if you didn’t install these libraries yet:

```shell
sudo apt-get install libcurl4-gnutls-dev zlib1g-dev
```
Then repeat the installation of RGT again.

Some errors are also due to older pip versions, so make sure to keep that updated:
```shell
pip install --user pip -U
```

## MacOS X
Before installing via pip, we need to setup the environment. First of all, you need the basic command line utilities for compilation:
```shell
xcode-select --install
```
See [this guide](https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/) for more visual details. After this, you must install the Homebrew version of Python:
```shell
brew install python boost llvm
```

Finally, we can use pip to install RGT. We suggest you **do not use the –user
option**, otherwise you’ll have to add another directory to your PATH.
```shell
pip install cython numpy scipy
pip install RGT --no-binary RGT
```

<!-- ### Installation with Apple silicon 
When you are using Mac with Apple Silicon chip (M1/2), the installation can be complicated, since they are based on arm64.

We have 

there are two scenarios,  -->

## Windows
Sorry. We don’t support windows environment this moment.

## Remove multiple installation
If you have multiple installation of RGT accidentally, please repeat the following command until no more RGT installation is detected.
```shell
pip uninstall RGT
```

Then please follow the instruction above for installation again.
