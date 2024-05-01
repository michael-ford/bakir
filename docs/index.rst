.. kir-annotator documentation master file, created by
   sphinx-quickstart on Thu Feb 29 21:04:04 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BAKIR's documentation!
=========================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   alignment_variants
   annotation_utils
   common
   gene_identification
   io_utils
   main
   mapper_wrappers
   mapping_utils

.. automodule:: kir_annotator
   :members:
   :undoc-members:
   :show-inheritance:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Command Line Interface
----------------------

.. argparse::
   :module: kir_annotator.io_utils
   :func: parse_arguments
   :prog: kir-annotator