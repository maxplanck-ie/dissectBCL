.. dissectBCL documentation master file, created by
   sphinx-quickstart on Mon Sep  5 14:49:49 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dissectBCL's documentation!
======================================

**dissectBCL** is a demultiplexing pipeline for Illumina sequencing data, loosely based on `TWTWTH <https://github.com/maxplanck-ie/TheWhoTheWhatTheHuh>`_.

The workflow is relatively simple:

 1. a new flowcell is found
 2. `parkour <https://github.com/maxplanck-ie/parkour2>`_ is queried and a sampleSheet is read
 3. `BCLConvert <https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html>`_ is used to demultiplex into samples.
 4. Re-name/organize project & samples, run a number of QC metrics.
 5. data is transferred to PI specific directories or uploaded with `fexsend <https://fex.rus.uni-stuttgart.de/usecases/fexsend.html>`_

Overview
========

.. toctree::
   :maxdepth: 3

   usage
   overview
