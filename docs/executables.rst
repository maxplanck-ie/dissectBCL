executables
===========

Installing dissectBCL will result in a couple of executables being added into your path. These are defined as entry points in the *setup.cfg* file.


#. :ref:`dissect <dissect>`
#. :ref:`wd40 <wd40>`
#. :ref:`email <email>`
#. :ref:`contam <contam>`

Help can be called for every executable using:

.. code-block:: console

    executable --help

.. _dissect:

dissect
^^^^^^^

dissect is the main pipeline function, and only has 2 optional arguments, which specifies the path to the :ref:`config.ini <config.ini>` file and (optional) a path to a flowcell to process.

.. code-block:: console

    dissect
    dissect -c config.ini
    dissect -f /path/to/flowcell

If no argument is specified, dissect looks for the config file in this path:

.. code-block:: console

    ~/configs/dissectBCL_prod.ini

.. _wd40:

wd40
^^^^

wd40 is a set of 'helper' tools, intended to make life slightly easier. At this point there are 3 helper functions. These are mainly work in progress, so don't expect to much from these (yet).

#. rel

rel
---

*rel* can be used to open up group permission (750) for :ref:`internal projects <Internals>`.
It can either be ran without arguments, which assumes the current working directory is a processed flow cell folder (written in the :ref:`outputDir <Dirs>`), or you can specify the path as a positional argument:

.. code-block:: console

    wd40 rel path/to/folder

Upon execution the fraction of files that have their rights changed will be printed per folder within a project.

.. _email:

email
^^^^^

*email* notifies the :ref:`internal <Internals>` end user of a released project. It takes a project folder as a positional argument, and has a couple of other options:

#. --configfile: path to a configfile (default = ~/configs/dissectBCL_prod.ini)
#. --notGood: flag that omits 'quality was good' string in the email.
#. --analysis: flag that specifies that `BRB <https://github.com/maxplanck-ie/BigRedButton>` did an analysis for this project
#. --cc: argument to include another email address in cc.
#. --comment: include a string in the email
#. --fromPerson: Name of the person taking care of the data. (e.g. Max)
#. --fromEmail: Email address of the person taking care of the data.
#. --fromSignature: path to a txt file with an email signature
#. --toEmail: email of the receiver.
#. --toName: name of the receiver.

The end user will be inferred by either setting it explicitely (--toEmail), or if not specified by querying parkour.
Since this command is used quite often, it can be beneficial to alias this command to something relevant for you:

.. code-block:: console

    email is aliased to `email --fromPerson Max --fromEmail mustermann@uni.de --fromSignature /path/to/max/signature.txt `

In which case an email could be sent with:

.. code-block:: console

    email Project_200_doe_john
    email --comment "This data is contaminated" Project_200_doe_john
    email --analysis --comment "This data is contaminated, but also analysed!" Project_200_doe_john

.. _contam:

contam
^^^^^^

*contam* is an executable that builds a kraken2 database from a yaml file. It has two required arguments:

#. -c / --contaminome: path to a contaminome yaml file.
#. -o / --outputdir: output directory to write the database into.

and one optional argument:

#. -t / --threads: number of threads (default = 15)


Note that we use a 'custom' taxonomical hierarchy, to simplify the output and to make sure we don't have to download the full taxdump database from NCBI.
It's organised as followed:

.. code-block:: console

    root (1) (no rank)
    |---|alive (2) (no rank)
    |---|---|eukaryote (4) (domain)
    |---|---|---|humangrp (9) (family)
    |---|---|---|mousegrp (10) (family)
    |---|---|---|flygrp (11) (family)
    |---|---|---|eugrp (13) (family)
    |---|---|prokaryote (5) (domain)
    |---|---|---|pseudomonasgrp (12) (family)
    |---|---|---|progrp (14) (family)
    |---|non-alive (3) (no rank)
    |---|---|phage (6) (domain)
    |---|---|virus (7) (domain)
    |---|---|vector (8) (domain)

And all the specific organisms are either part of a domain or a family.

For eukaryotes, the mitochondrial genome is excluded from the genome, and rRNA sequences are masked.
This is hardcoded in the prep_contaminome.py file, under the ignore_chrs dictionary and rrna_mask list, respectively.
In case you deviate from the provided contaminome.yml file, make sure to update these two variables if necessary.
If you update the contaminome.yml file, you *have* to update the taxmap dictionary, which has following structure:

`vulgarname: [taxid, parent_taxid, taxonomic level]`