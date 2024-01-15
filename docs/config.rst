.. _config.ini:

configfile
==========

Headers
^^^^^^^

The configfile is an `ini <https://en.wikipedia.org/wiki/INI_file>`_ file which is separated in different blocks:

#. :ref:`Dirs <Dirs>`
#. :ref:`Internals <Internals>`
#. :ref:`parkour <parkour>`
#. :ref:`software <software>`
#. :ref:`misc <misc>`
#. :ref:`communication <communication>`

.. _Dirs:

Dirs
----

The *Dirs block* defines path information to important directories.

#. baseDir: the base directory where the sequencer writes the output into.
#. outputDir: the directory where demultiplexing will be performed.
#. flowLogDir: the directory where dissectBCL will write its log files into.
#. seqFacDir: the directory where the sequencing facility has access to. Lightweight QC files will be written here.
#. piDir: The base directory that holds each principal investigator's (PI) folder (See :ref:`PIs <PIs>`).

.. _Internals:

Internals
---------

The *Internals block* defines which PI is internal. Upon completion, projects are either copied into the 'periphery' or uploaded via fexsend so external users can download the project.
Inside this block there are two elements:

.. _PIs:

#. PIs: a list of principal investigators.
#. seqDir: the directory inside a PI's directory where the sequencing data can be deposited.
#. fex: Boolean that indicates if an external project (PI not in PIs list) should be packed as a tar and uploaded using fexsend.

If a project is from an internal PI, it will be copied over into:

piDir/PI/seqDir

Note that multiple seqDirs per PI are allowed. For example if seqDir = sequencing_data:

#. sequencing_data
#. sequencing_data1
#. sequencing_data2

can exist, and the latest (e.g. the one with the highest number) will be used to copy over the data.


.. _parkour:

parkour
-------

The *parkour block* contains all necessary information to communicate with `parkour <https://github.com/maxplanck-ie/parkour2>`.
Note that this block contains sensitive information.

#. pullURL: the URL to pull flowcell information from. Is parkoururl/api/analysis_list/analysis_list
#. pushURL: the URL to push flowcell statistics to. Is parkoururl/api/run_statistics/upload
#. user: the username for API requests
#. pw: the password for API requests
#. cert: the pem certificate for API requests
#. URL: the URL to Parkour2, `https://` is implicit!

.. _software:

software
--------

The *software block* contains paths to all the necessary software and files that are *NOT* included in the conda installation.

#. bclconvert: path to the bcl-convert executable
#. fastqc_adapters: a (custom) list of adapters used by fastqc.
#. kraken2db: path to your kraken database (created with `contam`, or sourced from `elsewhere <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>`)

.. _misc:

misc
----

the *misc block* only contains one item (for now), which is the png file used in the custom multiqc header:

#. mpiImg: path to jpg file.

.. _communication:

communication
-------------

The *communication block* has four elements, all of which are related to email communication by the pipeline.

#. fromAddress: the e-mail address where the emails come from.
#. host: the email `host <https://docs.python.org/3/library/smtplib.html>`
#. finishedTo: email address(es) to send a notification upon completion of a flowcell. If multiple emails, these are comma separated.
#. bioinfoCore: email address of the core unit, where error messages go to.


example
^^^^^^^

.. code-block:: console

    [Dirs]
    baseDir=/path/to/bcl/folder
    outputDir=/path/to/fastq/output/folder
    flowLogDir=/path/to/log/folder
    seqFacDir=/path/to/share/qc/with/facility
    piDir=/base/with/enduser/folders
    bioinfoCoreDir=/path/to/share/qc/with/core

    [Internals]
    PIs=[pi1,pi2,pi3,pi4,pi5]
    seqDir=seqfolderstr

    [parkour]
    pullURL=parkour.pull.url/api/analysis_list/analysis_list
    pushURL=parkour.push.url/api/run_statistics/upload
    user=parkourUser
    password=parkourPw
    cert=/path/to/cert.pem
    URL=parkour.domain.tld

    [software]
    bclconvert=/path/to/bclconvert
    fastqc_adapters=/path/to/fastqc_adapters.txt
    kraken2db=/path/to/kraken2_contaminome/contaminomedb

    [misc]
    mpiImg=/path/to/multiqc_headerimg.jpg

    [communication]
    deepSeq=email@seqfacility.de
    bioinfoCore=email@bioinfocore.de
    fromAddress=sender@dissectbcl.de
    host=hostmail.address.de
