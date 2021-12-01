![flake8](https://github.com/WardDeb/dissectBCL/actions/workflows/flake/badge.svg)

# dissectBCL
demultiplexing pipeline.

bin/dissect calls dissectBCL/dissect.py, which is the main workflow.

Structure:  

dissectBCL module:  
  - classes.py: flowcell / sampleSheet class definitions
  - misc.py: various helper functions
  - fakeNews.py: home of the log system, as well as every 'communication' functions (email, API, ... ?)
  - drHouse: Future home for debug/diagnose barcode issues.
  - demux.py: Determine mismatches, create bclConvert demuxSheet and bclConvert runner
  - postmux.py: future home of the postprocessing.

tests:
  - future home for function tests


Install:

>  pip install ./
