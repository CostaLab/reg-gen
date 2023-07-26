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

# Installation with conda

We recommend using conda to manange the python environmont to avoid issues.

You can install conda from [here](https://docs.conda.io/en/latest/miniconda.html)

Once you successfully installed conda, first create a specific environment:

```shell
conda create -n rgt python=3.9
```

Then activate your environment and install the full RGT suite with all other dependencies:
```shell
conda activate rgt
pip install RGT
```

Detailed installation instructions and basic problem solving can be found on our [website](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html).

Please also consider citing our main paper if you used any sub-tools from RGT:
```
@article{li2023rgt,
  title={RGT: a toolbox for the integrative analysis of high throughput regulatory genomics data},
  author={Li, Zhijian and Kuo, Chao-Chung and Ticconi, Fabio and Shaigan, Mina and Gehrmann, Julia and Gusmao, Eduardo Gade and Allhoff, Manuel and Manolov, Martin and Zenke, Martin and Costa, Ivan G},
  journal={BMC bioinformatics},
  volume={24},
  number={1},
  pages={1--12},
  year={2023},
  publisher={BioMed Central}
}
```


