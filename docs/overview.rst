workflow
========

Find unprocessed flowcells.
---------------------------

Upon execution, the pipeline will look for unprocessed flowcells. 
An unprocessed flowcell has two characteristics:

 - a directory in *config[Dirs][baseDir]* that contains an *RTAComplete.txt* file.
 - no matching directories under *config[Dirs][outputDir]* that contain either a *fastq.made* file, or a *communication.done* file

Note that we split up a flowcell in lanes whenever we can (you can usually set a higher MisMatchIndex that way, retrieving more reads/sample).
This means that in *config[Dirs][baseDir]* we can have flowcell directory:

.. code-block:: console
    
    220101_A00000_0000_AXXXXXXXXX

and in *config[Dirs][outputDir]*

.. code-block:: console
    
    220101_A00000_0000_AXXXXXXXXX_lanes_1
    220101_A00000_0000_AXXXXXXXXX_lanes_2

if *fastq.made* or *communication.done* exists in **either** the first or the second folder, 220101_A00000_0000_AXXXXXXXXX will be set as **processed**.
This is important in case you need to re-run flowcells: All flags need to be removed for **both** folders.

demultiplex unprocessed flowcells.
----------------------------------

Once an unprocessed flowcell is found, a couple of steps are done before we start demultiplexing:

 1. Initiate a logfile *config[Dirs][logDir]*
 2. create a *flowcell class*
 3. create a *sampleSheet class*
 4. run bcl-convert
 5. run post processing steps
     1. fastqc
     2. fastqc_screen
     3. clumpify
     4. multiqc
 6. copy over the data to the *periphery* or upload to fexsend
 7. create a *drhouse class* and communicate results.

classes
-------

flowcell class
^^^^^^^^^^^^^^

This will be initiated almost immediately after a new flowcell is found.
It will contain the paths taken from the config file and set the paths to samplesheet.csv & runInfo.xml.
A check is ran that all these paths exist, and runinfo is parsed as well.

samplesheet class
^^^^^^^^^^^^^^^^^

After the flowcell class is made, the samplesheet class is created.
Here:

- the samplesheet is read
- parkour is queried
- the samplesheet df & parkour df are merged
- the decision is made to split up a flowcell in lanes or not


drhouse class
^^^^^^^^^^^^^

For every lane outputted, a drHouse class is made that contains information about the run:

- undetermined number of reads
- undetermined barcodes
- diskspace free
- etc...

Here the final email to notify about the run is also created.
