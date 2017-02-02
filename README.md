# RGT - Regulatory Genomics Toolbox

RGT is an open source python library for analysis of regulatory genomics. RGT is programmed in an oriented object fashion and its core classes provide functionality for handling regulatory genomics data.

This library has been used for implementation of several tools as ChIP-Seq differential peak callers ([ODIN](http://www.regulatory-genomics.org/odin-2/) and [THOR](http://www.regulatory-genomics.org/thor-2/)), DNase-Seq footprinting method ([HINT](http://www.regulatory-genomics.org/hint/)) and the visualization tool [RGT-Viz](http://www.regulatory-genomics.org/rgt-viz/).

# Quick installation

The quickest and easiest way to get RGT is to to use pip:

```
pip install RGT
```

This will install the full RGT suite with all dependencies.

Alternatively, you can clone this repository:
```
git clone https://github.com/CostaLab/reg-gen.git
```

or download a specific [release](https://github.com/CostaLab/reg-gen/releases), then proceed to manual installation:

```
cd reg-gen
sudo python setup.py install
```

Detailed installation instructions and basic problem solving can be found at:

http://www.regulatory-genomics.org/rgt/download-installation
