## Download & Installation

These instructions are only to be attempted if you know what you are doing and the basic instructions failed. Please contact us on our Support group if you experience any problems. Please notice that RGT only supports Python 3 now.

### Linux Installation

The recommended way of installing RGT is via [pip](https://pip.pypa.io/en/latest/installation/), the main python package manager. After installing pip, type in this order:

```shell
pip install --user cython numpy scipy
pip install --user RGT --no-binary RGT
```

This will install the full RGT suite with all python dependencies. A few more libraries will be needed, but should be installed via your own distribution’s package manager.