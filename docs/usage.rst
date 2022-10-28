Usage
=====

Installation
------------

There are a couple of pre-requisites

 1. `conda <https://docs.conda.io/en/latest/miniconda.html>`_
 2. `BCL-convert <https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html>`_
 3. A running `parkour <https://github.com/maxplanck-ie/parkour2>`_ instance with API access.


To install dissectBCL, first clone the repository:

.. code-block:: console

    git clone https://github.com/maxplanck-ie/dissectBCL
    cd dissectBCL

next create the conda environment. By default this will be named dissectBCL

.. code-block:: console

    conda env create -f env.yml --name dissectBCL

or if you have mamba installed:

.. code-block:: console

    mamba env create -f env.yml --name dissectBCL

activate the environment and pip install dissectBCL

.. code-block:: console

    conda activate dissectBCL
    pip install .

The next task is to prepare the contaminome database. These files will get downloaded so make sure to run the command on a node that has internet access.
Note that the contaminome.yml file is included in the repository. Note that it takes a while to create this database, and that you need quite a lot of memory.
The final footprint of this database is around 20GB.

.. code-block:: console

    contam -c contaminome.yml -o /path/to/folder

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

Note that to run persistently on a server, dissect should be run from tmux or using `nohup dissect &`. 
All other customisations have to happen in the ini file. Once running, the pipeline will check every hour. 
Forcing a check can be done by killing the HUP:

.. code-block:: console

    killall -HUP dissect


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

    wd40 rel

in the outlane folders.
once this is done, the end user can be notified using

.. code-block:: console

    email

Barcode issues
^^^^^^^^^^^^^^
Often, the biggest issues encountered will be wrong barcodes. An indication of this can be:

- low actual vs requested ratios
- high undetermined indices

Entry points here would be the email received, cross-referenced with outlanefolder/Reports/Top_Unknown_Barcodes.csv and outlanefolder/demuxSheet.csv
You could get additional information by running

.. code-block:: console 

    wd40 diag

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

and rerun dissectBCL. Note that an existing demuxSheet in the folder won't be overwritten, allowing you to jump in.

Other issues
^^^^^^^^^^^^
It can happen that the pipeline just crashes. A point of entry there would be to have a look at the log files. These are written per flowcell.
The folder in which these are written is specified in the ini file *config[Dirs][flowLogDir]*. 
Warnings in the log file usually correspond to what module is invoked, and Info tags show what is actually being done. 
Cross-referencing this information with the code can give you information on where to start debugging.
