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

Requirements

some software (defined in dissectBCL.ini) is not included in the installation.
You should have:
 - fastqc
 - multiqc
 - splitFastq
 - fastq_screen

 defined in the ini file.