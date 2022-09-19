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

    conda env create -f env.yaml

or if you have mamba installed:

.. code-block:: console

    mamba env create -f env.yaml

activate the environment and pip install dissectBCL

.. code-block:: console

    conda activate dissectBCL
    pip install .

Next task is to set up the ini file. A template is provided in the repository (*dissectBCL.ini*). It is *necessary* that all variables are filled in appropriately.
By default the pipeline expects this ini file to be available under:

.. code-block:: console

    ~/configs/dissectBCL_prod.ini

Although you can always override this location using the executable.

quickstart
----------

Starting the pipeline is as simple as:

.. code-block:: console

    dissect

or with a custom configfile location:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini

That's it, all other customisations have to happen in the ini file. 
Once running, the pipeline will check every hour. Forcing a check can be done by killing the HUP:

.. code-block:: console

    killall -HUP dissect


hands-on
--------

Say a flow cell has been processed. A first point of entry would be to look at the email received:

- All samples have good 'gotten' vs. 'requested' ratio's (~=1) ?
- what's the percentage of undetermined reads ?
- what are the top unknown barcodes ?
- how are we doing on space ?
- are the fqScreen organism and parkour organism the same ?

Next, have a look at the multiqc files (1 per project). These get copied over into *config[Dirs][bioinfoCoreDir]*.
Important here are:

- phred scores
- read composition
- detailed fastq_screen report

If everything looks fine, touch *fastq.made* into the lane folders and let `BigRedButton <https://github.com/maxplanck-ie/BigRedButton>`_ do it's job.

Afterwards, the folders in the *periphery* can be released by running:

.. code-block:: console

    wd40 rel

in the outlane folders.
once this is done, the end user can be notified using

.. code-block:: console

    email

Barcode issues
^^^^^^^^^^^^^^
Often, the biggest issues encountered will be wrong barcodes. An indication of this can be:

- low got vs requested ratios
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

and rerun dissectBCL.