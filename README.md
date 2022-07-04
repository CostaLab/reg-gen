[![Stars](https://img.shields.io/github/stars/CostaLab/reg-gen?logo=GitHub&color=yellow)](https://github.com/CostaLab/reg-gen/stargazers)
[![PyPI](https://img.shields.io/pypi/v/rgt?logo=PyPI)](https://pypi.org/project/RGT/)
[![PyPIDownloads](https://static.pepy.tech/badge/rgt)](https://static.pepy.tech/badge/rgt)
[![Docs](https://readthedocs.org/projects/reg-gen/badge/?version=latest)](https://reg-gen.readthedocs.io)

# RGT - Regulatory Genomics Toolbox

RGT is an open source Python 3.6+ library for analysis of regulatory genomics. RGT is programmed in an oriented object fashion and its core classes provide functionality for handling regulatory genomics data.

The toolbox is made of a core library and several tools:

* [HINT](https://reg-gen.readthedocs.io/en/latest/hint/introduction.html): ATAC-seq/DNase-seq footprinting method
* [THOR](https://reg-gen.readthedocs.io/en/latest/thor/introduction.html):
ChIP-Seq differential peak caller
* [Motif Analysis](https://reg-gen.readthedocs.io/en/latest/motif_analysis/introduction.html): TBFS match and enrichment
* [RGT-Viz](https://reg-gen.readthedocs.io/en/latest/rgt-viz/introduction.html): Visualization tool
* [TDF](https://reg-gen.readthedocs.io/en/latest/tdf/introduction.html): DNA/RNA triplex domain finder

See https://reg-gen.readthedocs.io for documentation and tutorials.

# Installation

The quickest and easiest way to get RGT is to to use pip. First some dependencies:

```shell
pip install --user cython numpy scipy
```

Then install the full RGT suite with all other dependencies:
```shell
pip install --user RGT
```

Detailed installation instructions and basic problem solving can be found on our [website](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html).

