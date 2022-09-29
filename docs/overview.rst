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
     2. kraken
     3. clumpify
     4. multiqc
 6. copy over the data to the *periphery* or upload to fexsend
 7. create a *drhouse class* and communicate results.


kraken
------

A little note on kraken. The custom database is prepared using the *contam* executable.
Inside `src/tools/prep_contaminome.py` some data structures are defined:

 1. ignore_chrs (dictionary)
 2. rrna_mask (list of tuples)
 3. taxmap (dictionary)

Ignore chrs have 'vulgar names' (defined in contaminome.yml) as keys, and a list of contig IDs as values.
These headers will *not* be included in the final database.

rrna mask is a list of tuples of which each tuple contains two vulgar names.
The idea is that vulgar name 1 sequences will be masked in vulgar name 0.
For example, in *(vulgar name 0, vulgar name 1)* a blast database will be created of *vulgar name 0*. Sequences in *vulgar name 1* will be blasted against this database, and resulting hits will be masked inside *vulgar name 0*.
The idea behind this is to ensure that hits against *vulgar name 1* won't give hits to *vulgar name 0* as well.
This is used to mask rRNA sequences, but can actually be used for anything.

Finally the taxmap is a dictionary with taxonomy names as keys, and a list with `[taxid, parent taxid, taxonomic rank]` as values.
A custom one is created as to not rely on a full NCBI taxonomy database dump.

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
