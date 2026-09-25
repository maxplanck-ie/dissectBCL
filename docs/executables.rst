executables
===========

Installing dissectBCL adds four executables to ``PATH``. Their entry points
are defined in the ``[project.scripts]`` table of ``pyproject.toml``:

* ``dissect``
* ``wd40``
* ``email``
* ``contam``

Every executable provides a help option. Use ``-h`` or ``--help`` with
``dissect``, ``email``, and ``contam``; ``wd40`` provides the same options on
its top-level command and on every subcommand.

.. _dissect:

dissect
^^^^^^^

``dissect`` is the main polling demultiplexing pipeline.

.. code-block:: console

    dissect [OPTIONS]

Options
-------

* ``-c PATH, --configfile PATH`` — existing configuration file. The default
  is ``~/configs/dissectBCL_prod.ini``.
* ``-f PATH, --flowcellpath PATH`` — process this flowcell directory instead
  of waiting for a new flowcell. The path is not required to exist when the
  command starts, but it must identify a flowcell directory when processing.
* ``-s {illumina,aviti}, --sequencer {illumina,aviti}`` — restrict processing
  to the selected platform and read only that platform's configuration keys.
* ``-F, --forcelanesplit`` — force lane splitting even when the sample sheet
  does not request it.
* ``-h, --help`` — show the command help and exit.

``dissect`` checks for new flowcells hourly and continues running until it is
stopped. Use ``--sequencer`` together with ``--flowcellpath`` when processing
a specific Illumina or Aviti directory, for example:

.. code-block:: console

    dissect -c /path/to/dissectBCL.ini -f /path/to/flowcell -s aviti
    dissect -c /path/to/dissectBCL.ini -F

.. _wd40:

wd40
^^^^

``wd40`` provides operational helper commands for a processed flowcell.

.. code-block:: console

    wd40 [OPTIONS] COMMAND [ARGS]...

Top-level options
-----------------

* ``--configpath PATH`` — existing configuration file. The default is
  ``~/configs/dissectBCL_prod.ini``.
* ``--debug / --no-debug`` (also ``-d / -n``) — enable or disable debug
  logging. The default is ``--no-debug``.
* ``--version`` — show the installed dissectBCL version and exit.
* ``-h, --help`` — show the ``wd40`` help and exit.

The top-level ``--help`` and ``--version`` options, as well as ``wd40 help``,
work without a config file. The other ``wd40`` subcommands load the
configuration before running and therefore need an existing, valid file.

``wd40 rel``
------------

Release a finished flowcell to the periphery and report its paths to Parkour2.
It sets the expected group ownership and mode 750 on the flowcell, project,
FASTQC, and Analysis folders. The command warns when BigRedButton has not
created ``analysis.done``.

.. code-block:: console

    wd40 rel [FLOWCELL] [--force PROJECT,PI]

* ``FLOWCELL`` — processed flowcell/outLane directory. The default is ``./``.
* ``-h, --help`` — show the ``rel`` help and exit.
* ``--force PROJECT,PI`` — force-ship only
  ``Project_<PROJECT>_<user>_<PI>`` to the latest sequencing-data volume for
  the explicit ``PI`` before releasing that project. The project ID must be
  numeric and the PI is a simple directory name.

Examples:

.. code-block:: console

    wd40 rel /path/to/outLane
    wd40 rel /path/to/outLane --force=4070,iovino

``wd40 reset``
--------------

Reset a processed outLane so its sample sheet or run manifest can be edited
and the lane can be demultiplexed again. The command lists the files it will
delete and asks for confirmation (the confirmation defaults to no). It keeps
``demuxSheet.csv`` for Illumina or ``manifest/RunManifest.csv`` for Aviti, and
removes demultiplexing output and completion flags. If neither expected
manifest is present, it aborts; if there is nothing to remove, it reports that
there is nothing to reset.

.. code-block:: console

    wd40 reset [OUTLANE]

* ``OUTLANE`` — outLane directory. The default is ``./``.
* ``-h, --help`` — show the ``reset`` help and exit.

``wd40 fex``
------------

Upload a dissectBCL project to FEX as an RO-Crate archive. The project name
must have the form ``Project_<request_id>_<user>_<PI>``. The archive is sent
to ``fexsend`` without first being written to disk.

.. code-block:: console

    wd40 fex PROJECT [OPTIONS]

* ``PROJECT`` — existing ``Project_<request_id>_<user>_<PI>`` directory.
* ``--parkour-url URL`` — override the Parkour base URL. When omitted, the URL
  from the configuration's ``[parkour] URL`` setting is used.
* ``-h, --help`` — show the ``fex`` help and exit.

``fexsend`` must be available in ``PATH`` (or at the standard user-local
fallback). For example:

.. code-block:: console

    wd40 fex Project_3358_Hohl_Manke
    wd40 fex --parkour-url https://parkour-test.ie-freiburg.mpg.de Project_3358_Hohl_Manke

``wd40 help``
-------------

List the available ``wd40`` subcommands and a short description of when to
use each one.

.. code-block:: console

    wd40 help

* ``-h, --help`` — show the ``help`` command help and exit.

.. _email:

email
^^^^^

Send a notification about one or more finished projects. Run ``email`` from
the flowcell outLane directory: the directory name is used to locate the
sequencing-data folder in the message.

.. code-block:: console

    email [OPTIONS] PROJECT [PROJECT ...]

Options
-------

* ``--configfile CONFIGFILE`` — configuration file. The default is
  ``~/configs/dissectBCL_prod.ini``; the file must be readable and valid.
* ``--notGood`` — omit the statement that the sequencing quality was good.
* ``--analysis`` — mention that BigRedButton performed an analysis.
* ``--cc ADDRESS [ADDRESS ...]`` — add one or more CC recipients.
* ``--comment COMMENT`` — add a comment string, or the contents of a file.
* ``--fromPerson NAME`` — name of the sender.
* ``--fromEmail ADDRESS`` — sender email address. The sender receives a BCC.
* ``--fromSignature PATH`` — optional signature file.
* ``--toEmail ADDRESS`` — recipient email address. If ``--toEmail`` or
  ``--toName`` is omitted, the recipient is looked up in Parkour.
* ``--force-to=EMAIL,PI`` — use the specified recipient and sequencing-data
  PI instead of Parkour contact lookup, resolve the latest matching
  ``sequencing_data*`` directory for that PI, and append the FEX download URL
  retrieved with ``fexsend -l`` to the comment block. Use the exact form
  ``EMAIL,PI``; for example,
  ``--force-to=mendelevich@ie-freiburg.mpg.de,iovino``.
* ``--toName NAME`` — name of the recipient.
* ``-h, --help`` — show the command help and exit.

One or more project directories may be supplied. ``--fromPerson`` and
``--fromEmail`` are required at runtime. The configured
``[communication] bioinfoCore`` address is always BCC'd, and the configured
``[communication] host`` is used to send the message. When multiple projects
are listed, the contact for the first project receives the email and all
project IDs are included in the body. Projects delivered externally via FEX
are not supported by the normal path-resolution mode.

Examples:

.. code-block:: console

    email --fromPerson Core --fromEmail core@example.org Project_200_doe_john
    email --comment "This data is contaminated" Project_200_doe_john
    email --force-to=mendelevich@ie-freiburg.mpg.de,iovino Project_200_doe_iovino

Because the sender options are commonly reused, they can be supplied through
an alias. For example:

.. code-block:: console

    alias email='email --fromPerson Core --fromEmail core@example.org --fromSignature /path/to/signature.txt'

.. _contam:

contam
^^^^^^

Build a Kraken2 contaminome database from a YAML file.

.. code-block:: console

    contam [OPTIONS] -c CONTAMINOME -o OUTPUTDIR

Options
-------

* ``-c PATH, --contaminome PATH`` — required YAML contaminome specification.
* ``-o PATH, --outputdir PATH`` — required, existing output directory.
* ``-t THREADS, --threads THREADS`` — number of worker threads. The default
  is ``15``.
* ``-f, --force`` — remove and recreate an existing
  ``<outputdir>/contaminomedb``. Without this option, an existing database
  causes the command to stop.
* ``-h, --help`` — show the command help and exit.

Example:

.. code-block:: console

    contam --threads 10 -c contaminome.yml -o /path/to/existing/folder
    contam --threads 10 --force -c contaminome.yml -o /path/to/existing/folder

The build requires network access to the genome URLs in the YAML file and
the external tools ``makeblastdb``, ``blastn``, ``bedtools maskfasta``, and
``kraken2-build``. Existing downloaded genome files are reused on later runs;
``--force`` controls replacement of the contaminome database itself.

Screening notes
---------------

The custom taxonomical hierarchy and chromosome-filtering behavior used by
``contam`` are described in the :ref:`workflow overview <kraken>`. The
screening configuration and PlusPF escalation behavior are described in the
:ref:`screening configuration <screening>`.
