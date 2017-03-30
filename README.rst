RGT - Regulatory Genomics Toolbox
=================================

RGT is an open source python library for analysis of regulatory
genomics. RGT is programmed in an oriented object fashion and its core
classes provide functionality for handling regulatory genomics data.

The toolbox is made of a `core library <http://www.regulatory-genomics.org/rgt/>`__ and several tools:

* `THOR <http://www.regulatory-genomics.org/thor-2/>`__: ChIP-Seq differential peak caller, replaces
  `ODIN <http://www.regulatory-genomics.org/odin-2/>`__

* `Motif Analysis <http://www.regulatory-genomics.org/motif-analysis/>`__: TBFS match and enrichment

* `HINT <http://www.regulatory-genomics.org/hint/>`__: DNase-Seq footprinting method

* `RGT-Viz <http://www.regulatory-genomics.org/rgt-viz/>`__: Visualization tool

* `TDF <http://www.regulatory-genomics.org/tdf/>`__: DNA/RNA triplex domain finder

Installation
============

The quickest and easiest way to get RGT is to to use pip:

::

    pip install --user RGT

This will install the full RGT suite with all dependencies.
If you have errors, this is usually due to three dependencies
that you can install beforehand:

::

    pip install --user cython numpy scipy

Alternatively, you can clone this repository:

::

    git clone https://github.com/CostaLab/reg-gen.git

or download a specific
`release <https://github.com/CostaLab/reg-gen/releases>`__, then proceed
to manual installation:

::

    cd reg-gen
    python setup.py install --user

Detailed installation instructions and basic problem solving can be
found at:

http://www.regulatory-genomics.org/rgt/download-installation
