.. RGT documentation master file, created by
   sphinx-quickstart on Mon Jul  4 13:14:44 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

RGT - Regulatory Genomics Toolbox
===================================

RGT is an open source Python 3.6+ library for analysis of regulatory genomics. 
It is programmed in an oriented object fashion and its core classes provide functionality for handling regulatory genomics data.

This project is developed and maintained by `Institute for Computational Genomics <http://www.costalab.org/>`__ , RWTH University Hostpital. 
If you have any questions/comments/problems, please open an issue on `Github <https://github.com/CostaLab/reg-gen>`__.

The toolbox is made of a core library and several tools:

- `HINT <http://www.regulatory-genomics.org/hint/>`__: ATAC-seq/DNase-seq footprinting method
- `THOR <http://www.regulatory-genomics.org/thor-2/>`__: ChIP-Seq differential peak caller
- `Motif Analysis <http://www.regulatory-genomics.org/motif-analysis/>`__: TBFS match and enrichment
- `RGT-Viz <http://www.regulatory-genomics.org/rgt-viz/>`__: Visualization tool
- `TDF <http://www.regulatory-genomics.org/tdf/>`__: DNA/RNA triplex domain finder


.. toctree::
   :caption: RGT
   :maxdepth: 2
   :hidden:

   RGT/installation.md
   RGT/setup_data.md

.. toctree::
   :caption: HINT
   :maxdepth: 2
   :hidden:
   
   HINT/introduction.md
   HINT/tutorial-dc.md

.. toctree::
   :caption: THOR
   :maxdepth: 2
   :hidden:
   
   THOR/introduction.md

.. toctree::
   :caption: Motif analysis
   :maxdepth: 2
   :hidden:
   
   MotifAnalysis/introduction.md

.. toctree::
   :caption: TDF
   :maxdepth: 2
   :hidden:
   
   TDF/introduction.md

.. toctree::
   :caption: RGT-Viz
   :maxdepth: 2
   :hidden:
   
   RGT-Viz/introduction.md
