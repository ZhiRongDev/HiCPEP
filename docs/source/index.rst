.. 
    HiCPAP documentation master file, created by
    sphinx-quickstart on Tue Mar 26 23:09:02 2024.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

Welcome to HiCPAP's documentation!
==================================

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    Installation <install>
    Quick start <quick_start>
    Python API <api>
    Code for paper <code_for_paper>
    data_store directory structure <data_store>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

How to build this document
==========================

The document is created with `Sphinx documentation <https://www.sphinx-doc.org/en/master/>`_. 
All the HTML pages will be placed in the ``build`` directory.

.. code:: bash

    cd HiCPAP/docs
    make html