[![Documentation Status](https://readthedocs.org/projects/dissectbcl/badge/?version=latest)](https://dissectbcl.readthedocs.io/en/latest/?badge=latest)
[![Lint](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml/badge.svg)](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml)
![Pytest](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/pytest.yml/badge.svg)

# dissectBCL

Demultiplexing pipeline for illumina data (novaseq/miseq/nextseq). Continuation of Devon Ryan's [TWTWTWTW](https://github.com/maxplanck-ie/TheWhoTheWhatTheHuh).

## Installation.

Clone this repository, create the environment and pip install

 > git clone git@github.com:maxplanck-ie/dissectBCL.git  
 > cd dissectBCL  
 > conda env create -f env.yml --name dissectBCL  
 > conda activate dissectBCL  
 > pip install ./  

> [!NOTE]
> If you get "LookupError: setuptools-scm was unable to detect version". Check that git is available on `$PATH` and, on the base conda environment you've just activated.

## Running.

Fill in the dissectBCL.ini file appropriately. By default the config file is expected to be in ~/configs/dissectBCL_prod.ini.

 > dissect

or 

 > dissect -c /path/to/config.ini

or

 > dissect -f /path/to/flowcell.ini

## Docs.

Documentation is available [here](https://dissectbcl.readthedocs.io/en/latest/).
