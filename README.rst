RGT - Regulatory Genomics Toolbox
=================================

RGT is an open source python library for analysis of regulatory
genomics. RGT is programmed in an oriented object fashion and its core
classes provide functionality for handling regulatory genomics data.

This library has been used for implementation of several tools such as
ChIP-Seq differential peak caller
`THOR <http://www.regulatory-genomics.org/thor-2/>`__), DNase-Seq
footprinting method
(`HINT <http://www.regulatory-genomics.org/hint/>`__) and the
visualization tool
`RGT-Viz <http://www.regulatory-genomics.org/rgt-viz/>`__.

Quick installation
==================

The quickest and easiest way to get RGT is to to use pip:

::

    pip install --user RGT

This will install the full RGT suite with all dependencies.

Alternatively, you can clone this repository:

::

    git clone https://github.com/CostaLab/reg-gen.git

or download a specific
`release <https://github.com/CostaLab/reg-gen/releases>`__, then proceed
to manual installation:

::

    cd reg-gen
    sudo python setup.py install --user

Detailed installation instructions and basic problem solving can be
found at:

http://www.regulatory-genomics.org/rgt/download-installation
