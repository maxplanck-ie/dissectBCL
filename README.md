[![Documentation Status](https://readthedocs.org/projects/dissectbcl/badge/?version=latest)](https://dissectbcl.readthedocs.io/en/latest/?badge=latest)
[![Lint](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml/badge.svg)](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml)
![Pytest](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/pytest.yml/badge.svg)

# dissectBCL

Demultiplexing pipeline for illumina data (novaseq/miseq/nextseq). Continuation of Devon Ryan's [TWTWTWTW](https://github.com/maxplanck-ie/TheWhoTheWhatTheHuh).

## Installation.

Clone this repository and run the install script. It creates a version-named
conda env (e.g. `dissect_v1.0.3`) from the latest git tag and installs
dissectBCL into it, so the reported version always matches what's running
and older releases stay available side by side.

 > git clone git@github.com:maxplanck-ie/dissectBCL.git  
 > cd dissectBCL  
 > ./install_dissect.sh  
 > conda activate dissect_v1.0.3  

`git` is required in the environment (see `env.yml`) since `pip install` uses setuptools-scm to derive the package version from git metadata.

## Running.

Fill in the dissectBCL.ini file appropriately. By default the config file is expected to be in ~/configs/dissectBCL_prod.ini.

 > dissect

or 

 > dissect -c /path/to/config.ini

or

 > dissect -f /path/to/flowcell.ini

## Docs.

Documentation is available [here](https://dissectbcl.readthedocs.io/en/latest/).
