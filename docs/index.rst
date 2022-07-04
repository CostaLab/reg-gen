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

- `HINT <https://reg-gen.readthedocs.io/en/latest/hint/introduction.html>`__: ATAC-seq/DNase-seq footprinting method
- `THOR <https://reg-gen.readthedocs.io/en/latest/thor/introduction.html>`__: ChIP-Seq differential peak caller
- `Motif Analysis <https://reg-gen.readthedocs.io/en/latest/motif_analysis/introduction.html>`__: TBFS match and enrichment
- `RGT-Viz <https://reg-gen.readthedocs.io/en/latest/rgt-viz/introduction.html>`__: Visualization tool
- `TDF <https://reg-gen.readthedocs.io/en/latest/tdf/introduction.html>`__: DNA/RNA triplex domain finder


.. toctree::
   :caption: RGT
   :maxdepth: 2
   :hidden:

   rgt/installation.md
   rgt/setup_data.md

.. toctree::
   :caption: HINT
   :maxdepth: 2
   :hidden:
   
   hint/introduction.md
   hint/tutorial-dendritic-cell.md
   hint/tutorial-single-cell.md

.. toctree::
   :caption: Motif analysis
   :maxdepth: 2
   :hidden:
   
   motif_analysis/introduction.md
   motif_analysis/tutorial.md
   motif_analysis/tool_usage.md
   motif_analysis/additional_motif_data.md

.. toctree::
   :caption: THOR
   :maxdepth: 2
   :hidden:
   
   thor/introduction.md

.. toctree::
   :caption: TDF
   :maxdepth: 2
   :hidden:
   
   tdf/introduction.md

.. toctree::
   :caption: RGT-VIZ
   :maxdepth: 2
   :hidden:
   
   rgt-viz/introduction.md
