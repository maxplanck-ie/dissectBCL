Usage
=====

Installation
------------

There are a couple of pre-requisites

 1. `conda <https://docs.conda.io/en/latest/miniconda.html>`_
 2. `BCL-convert <https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html>`_ and/or `Bases2fastq <https://docs.elembio.io/docs/bases2fastq/introduction/>`_ in case you have aviti data.
 3. A running `parkour <https://github.com/maxplanck-ie/parkour2>`_ instance with API access.

To install dissectBCL, first clone the repository:

.. code-block:: console

    git clone https://github.com/maxplanck-ie/dissectBCL
    cd dissectBCL

next create the conda environment. By default this will be named dissectBCL

.. code-block:: console

    conda env create -f env.yml --name dissectBCL

activate the environment and pip install dissectBCL

.. code-block:: console

    conda activate dissectBCL
    pip install .

The next task is to prepare the contamination database. These files will get downloaded so make sure you have internet access.
Note that the contaminome.yml file is included in the repository. Note that it takes a while to create this database, and that you need quite a lot of memory (~20GB).
The final footprint of this database is around 30GB. The output directory must already exist. If the contaminome database is already present, rerun with ``--force`` only when you want to replace it.

.. code-block:: console

    contam --threads 10 -c contaminome.yml -o /path/to/folder
    contam --threads 10 --force -c contaminome.yml -o /path/to/folder

Next task is to set up the ini file. A template is provided in the repository (*dissectBCL.ini*). It is *necessary* that all variables are filled in appropriately.
By default the pipeline expects this ini file to be available under:

.. code-block:: console

    ~/configs/dissectBCL_prod.ini

Although you can always override this location using the executable.

**NB: if you are running Parkour, the config file will contain a 
plaintext password, so set permissions appropriately!**

Quickstart
----------

Starting the pipeline is as simple as:

.. code-block:: console

    dissect

or with a custom configfile location:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini

Forcing to run a specific flowcell can also be done via the command line. Use
``--sequencer`` to select the platform when processing a specific directory:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini -f /full/path/to/flowcell/directory -s illumina
    dissect -c /path/to/dissectBCL.ini -f /full/path/to/flowcell/directory -s aviti

Lane splitting
~~~~~~~~~~~~~~

dissectBCL decides per flow cell whether to demultiplex each lane separately
(``decideSplit``). Splitting is preferred: fewer barcodes per demultiplex run
allows more index mismatches. The decision is written to the log
(``decide_lanesplit returns True/False``). Splitting is switched off when:

* a sample is loaded on more than one lane
  (``No lane splitting: >= 1 sample in multiple lanes.``),
* a project is loaded on more than one lane
  (``No lane splitting: >= 1 project in multiple lanes.``), or
* the sample sheet lists fewer lanes than the run has
  (``No lane splitting: 1 lane listed, N found.``).

If splitting is off but two samples share the same index (index 1, or index 1
plus index 2), it is switched back on to avoid a barcode clash
(``Found a sample clash even though laneSplit == False.``).

Use ``--forcelanesplit`` (``-F``) to skip this decision and always split per
lane. The flag is off unless passed.

When to force the split
^^^^^^^^^^^^^^^^^^^^^^^

The automatic decision cannot tell whether a project spread over several lanes
holds the *same* library in each lane or *different* libraries. Force the split
when the lanes carry different libraries with different index layouts, for
example a Multiome project with one pool per lane.

Example: a flow cell has one project in lanes 1 and 2. Lane 1 is a standard
dual-index pool, lane 2 is a Multiome ATAC pool. Without ``-F`` the log shows
``No lane splitting: >= 1 project in multiple lanes.`` and
``decide_lanesplit returns False``. dissectBCL then demultiplexes both lanes
in one run (``<flowcell>_lanes_1_2``). Because lane 2 contains scATAC, the
whole run gets the scATAC mask ``Y101;I8N2;U24;Y101``. That mask does not fit
the 10 bp indices of lane 1, and BCLConvert exits with code 1 (``Demux -
BCLConvert exit 1``).

Re-running with ``-F`` demultiplexes each lane on its own
(``<flowcell>_lanes_1``, ``<flowcell>_lanes_2``). The mask is determined per
lane, from that lane's samples only. Masks below are for a run with 101 cycle
reads, 10 i7 cycles and 24 i5 cycles; other recipes give other numbers:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini -f /full/path/to/flowcell/directory -s illumina -F

====== =========================== ============================
Lane   Samples contain             Mask
====== =========================== ============================
1      nothing special             ``Y101;I10;I10N14;Y101``
2      an scATAC protocol          ``Y101;I8N2;U24;Y101``
====== =========================== ============================

Where the mask comes from
^^^^^^^^^^^^^^^^^^^^^^^^^

The sample sheet is generated by Parkour; besides the indices it holds the
library protocol in ``Description``. dissectBCL also pulls library details from
Parkour by flow cell ID (``pullParkour``) and merges them onto the sample sheet
(Illumina: by ``Sample_ID``, ``Sample_Name`` and ``Sample_Project``; Aviti: by
``SampleName`` and ``Project``). Two columns change the mask (``detMask``):

* ``indexType`` (the index type chosen in Parkour for the library): if any
  sample contains ``NuGEN Ovation SoLo RNA-Seq System``, the index read is
  masked as index plus UMI (``I<n>U<m>``).
* ``Description`` (the library protocol): if any sample
  contains one of ``scATAC-Seq 10xGenomics``, ``NextGEM_Multiome_ATAC``,
  ``Next GEM Single Cell ATAC`` or ``MUXscATAC-seq v3``, the scATAC mask is
  used.

Both checks use ``any`` over the samples that are demultiplexed together. One
special library therefore changes the mask for every sample in the run. With
lane splitting each lane is checked on its own, so a special library affects
only its lane.

Where ``Description`` is read from depends on the platform. Illumina keeps the
sample sheet's column, i.e. the protocol at the time the sample sheet was
generated. Aviti drops the manifest's column and uses the one pulled live from
Parkour. If the protocol or index type is wrong in Parkour, the mask is wrong
too. Fix it in Parkour, not in dissectBCL; on Illumina also regenerate the
sample sheet before re-running, or the old protocol is used.

The complete option list is available from ``dissect --help``.


API
---

As of 0.3.0 a flow cell can be processed purely over a python shell:

.. code-block:: python

    from dissectBCL.dissect import createFlowcell
    f = createFlowcell("/path/to/config.ini", "/path/to/flowcell/", "illumina")
    f.prepConvert()
    f.demux()
    f.postmux()
    f.fakenews()
    f.organiseLogs()

The third argument in createFlowcell is the sequencer type, and can be either 'illumina' or 'aviti'.
By default the logs are printed to stdout, but you can move them to a file as well.

.. code-block:: python

    from dissectBCL.dissect import createFlowcell
    f = createFlowcell("/path/to/config.ini", "/path/to/flowcell/", 'aviti', logFile = "/path/to/logfile")
    f.prepConvert()
    f.demux()
    f.postmux()
    f.fakenews()
    f.organiseLogs()


Hands-on
--------

Say a flow cell has been processed. A first point of entry would be to look at the email received:

- All samples have good 'actual' vs. 'requested' ratios (~=1)?
- what's the percentage of undetermined reads?
- what are the top unknown barcodes?
- how are we doing on space?
- are the kraken2 organism and parkour organism the same?


Next, have a look at the multiqc files (1 per project). These get copied over into *config[Dirs][bioinfoCoreDir]*.
Important here are:

- phred scores
- read composition
- detailed kraken2 report

If everything looks fine, touch *fastq.made* into the lane folders and let `BigRedButton <https://github.com/maxplanck-ie/BigRedButton>`_ do its job.

We assume that end users can access the files in the *periphery* by group rights, not with user rights.
'releasing' data in this case just chmod to 750.
The folders in the *periphery* can be released by running:

.. code-block:: console

    wd40 rel /path/to/outLane/folder

The release changes permissions to 750, and pushes back to parkour that the flow cell has been released. Some PIs in parkour are external collaborators (with a contract) of inside PIs. Their data normally skips the periphery (``/data/<pi>``) and goes to the FEX server. For such a project, ship it to the periphery of the inside PI first, with:

.. code-block:: console

    wd40 rel /path/to/outLane/folder --force=PROJECT,PI

The other ``wd40`` subcommands are documented in :ref:`wd40` and can be
listed with ``wd40 --help`` or ``wd40 help``. In particular,
``wd40 reset`` prepares an outLane for a redemux and ``wd40 fex`` uploads a
project as an RO-Crate archive.

Finally, you can notify the end user with the email functionality. Run the
command from the outLane directory and use ``email --help`` for the complete
option list:

.. code-block:: console

    email -h

Barcode issues
^^^^^^^^^^^^^^
Often, the biggest issues encountered will be wrong barcodes. An indication of this can be:

- low actual vs requested ratios
- high undetermined indices

Entry points here would be the email received, cross-referenced with
``outlanefolder/Reports/Top_Unknown_Barcodes.csv`` and the platform's
sample-sheet/run-manifest file (``demuxSheet.csv`` for Illumina or
``manifest/RunManifest.csv`` for Aviti).

Identify what (and if) changes can be made, back up the generated sheet, and
make the changes accordingly. For a supported reset, use:

.. code-block:: console

    wd40 reset /path/to/outLane/folder

After confirmation, ``wd40 reset`` removes the demultiplexing output and
completion flags while keeping the sample sheet or run manifest. It also
handles both Illumina and Aviti output layouts. Rerun dissectBCL after the
reset; an existing sheet is used as provided.

Issues with Parkour verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this case (which is rare as it's caused by changing the certificate provider and it is not commonly listed), the certificate issuer is not recognized as dissect throws this error:

.. code-block:: console

    requests.exceptions.SSLError: HTTPSConnectionPool(host='parkourURL', port=443): Max retries exceeded with url: 
    /api/analysis_list/analysis_list/?flowcell_id=XXXXXXXXX (Caused by SSLError(SSLCertVerificationError(1, '[SSL: CERTIFICATE_VERIFY_FAILED] certificate 
    verify failed: unable to get local issuer certificate (_ssl.c:1007)')))

then, the new certificate needs to be added in the system (i.e. for CentOS 7, copying it to /etc/pki/ca-trust/source/anchors/ and run "update-ca-trust").

Finally, the cert field under the parkour header in the configuration file needs to point to the file copied in the anchors directory.
    

Other issues
^^^^^^^^^^^^
It can happen that the pipeline just crashes. A point of entry there would be to have a look at the log files. These are written per flowcell.
The folder in which these are written is specified in the ini file, per platform, as *config[Dirs][flowLogDir_illumina]* / *config[Dirs][flowLogDir_aviti]*.
Warnings in the log file usually correspond to what module is invoked, and Info tags show what is actually being done. 
Cross-referencing this information with the code can give you information on where to start debugging.
