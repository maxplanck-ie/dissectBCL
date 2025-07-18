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
Note that the contaminome.yml file is included in the repository. Note that it takes a while to create this database, and that you need quite a lot of memory (~20GB)
The final footprint of this database is around 30GB. 

.. code-block:: console

    contam --threads 10 -c contaminome.yml -o /path/to/folder

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

Forcing to run a specific flowcell can also be done via the command line:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini -f /full/path/to/flowcell/directory


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

The release changes permissions to 750, and pushes back to parkour that the flow cell has been released.
Finally, you can notify the end user with the email functionality.

.. code-block:: console

    email -h

Barcode issues
^^^^^^^^^^^^^^
Often, the biggest issues encountered will be wrong barcodes. An indication of this can be:

- low actual vs requested ratios
- high undetermined indices

Entry points here would be the email received, cross-referenced with outlanefolder/Reports/Top_Unknown_Barcodes.csv and outlanefolder/demuxSheet.csv

Identify what (and if) changes can be made, backup the generated demuxSheet, and make changes accordingly.
After the changes have been made in the demuxSheet:

- remove the project/FASTQC folders in the periphery
- remove the project/FASTQC folders in the outlane folder(s)

remove all the flags:

- analysis.done
- bclconvert.done
- communication.done
- fastq.made
- postmux.done
- renamed.done

and rerun dissectBCL. Note that an existing demuxSheet in the folder won't be overwritten but used as provided.

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
The folder in which these are written is specified in the ini file *config[Dirs][flowLogDir]*. 
Warnings in the log file usually correspond to what module is invoked, and Info tags show what is actually being done. 
Cross-referencing this information with the code can give you information on where to start debugging.
