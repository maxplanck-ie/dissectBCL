workflow
========

Find unprocessed flowcells.
---------------------------

Upon execution, the pipeline will look for unprocessed flowcells. 
An unprocessed flowcell has two characteristics:

 - a directory in *config[Dirs][baseDir]* that contains an *RTAComplete.txt* file and an *CopyComplete.txt* file.
 - no matching directories under *config[Dirs][outputDir]* that contain either a *fastq.made* file, or a *communication.done* file

Note that we split up a flowcell in lanes whenever we can (you can usually set a higher MisMatchIndex that way, retrieving more reads/sample).
This means that in *config[Dirs][baseDir]* we can have flowcell directory:

.. code-block:: console
    
    220101_A00000_0000_AXXXXXXXXX

and in *config[Dirs][outputDir]*

.. code-block:: console
    
    220101_A00000_0000_AXXXXXXXXX_lanes_1
    220101_A00000_0000_AXXXXXXXXX_lanes_2

if *fastq.made* or *communication.done* exists in **either** the first or the second folder, 220101_A00000_0000_AXXXXXXXXX will be considered as **processed**.
This is important in case you need to re-run flowcells: All flags need to be removed for **both** folders.

demultiplex unprocessed flowcells.
----------------------------------

Once an unprocessed flowcell is found, a couple of steps are performed.

 1. Initiate a logfile *config[Dirs][logDir]*
 2. create the *flowcell class*
 3. prepConvert() - determine mismatches and masking.
 4. demux() - run demultiplexing with bclconvert.
 5. postmux() - run renaming of projects, clumping, fastqc, kraken, multiqc and md5sum calculation.
 6. fakenews() - upload project via fexsend (if applicable), collate quality metrics, create and send email.
 7. organiseLogs() - dump out configs and settings to the outLanes.

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