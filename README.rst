RGT - Regulatory Genomics Toolbox
=================================

.. class:: no-web no-pdf

|pypi| |dev_build| |coverage|

RGT is an open source Python 3.6+ library for analysis of regulatory
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

**Python 2 is no longer supported.**

The quickest and easiest way to get RGT is to to use pip. First some dependencies:

::

    pip install --user cython numpy scipy

Then install the full RGT suite with all other dependencies:

::

    pip install --user RGT


Alternatively (but not recommended), you can clone this repository:

::

    git clone https://github.com/CostaLab/reg-gen.git

or download a specific
`release <https://github.com/CostaLab/reg-gen/releases>`__, then proceed
to manual installation:

::

    cd reg-gen
    python setup.py install --user

Detailed installation instructions and basic problem solving can be
found `on our website <http://www.regulatory-genomics.org/rgt/download-installation>`__.

For any issues, please write to our `support mailing list <https://groups.google.com/forum/#!forum/rgtusers>`__.

.. |pypi| image:: https://img.shields.io/pypi/v/rgt.svg?label=latest%20release
    :target: https://pypi.python.org/pypi/rgt
    :alt: Latest version released on PyPi

.. |mast_build| image:: https://img.shields.io/travis/CostaLab/reg-gen.svg?branch=master&label=master
    :target: https://travis-ci.org/CostaLab/reg-gen
    :alt: Build status of the master branch

.. |dev_build| image:: https://img.shields.io/travis/CostaLab/reg-gen.svg?branch=develop&label=develop
    :target: https://travis-ci.org/CostaLab/reg-gen
    :alt: Build status of the develop branch

.. |coverage| image:: https://img.shields.io/coveralls/CostaLab/reg-gen/develop.svg?label=coverage
    :target: https://coveralls.io/r/CostaLab/reg-gen?branch=develop
    :alt: Test coverage
