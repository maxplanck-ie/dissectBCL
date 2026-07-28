[![Documentation Status](https://readthedocs.org/projects/dissectbcl/badge/?version=latest)](https://dissectbcl.readthedocs.io/en/latest/?badge=latest)
[![Lint](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml/badge.svg)](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/lint.yml)
![Pytest](https://github.com/maxplanck-ie/dissectBCL/actions/workflows/pytest.yml/badge.svg)

# dissectBCL

Demultiplexing pipeline for short read genome sequencing data (illumina or aviti.)

## Installation

Clone this repository and run the install script. It creates a version-named
conda env (e.g. `dissect_v1.0.3`) from the latest git tag and installs
dissectBCL into it, so the reported version always matches what's running
and older releases stay available side by side.

 > git clone git@github.com:maxplanck-ie/dissectBCL.git  
 > cd dissectBCL  
 > ./install_dissect.sh  
 > VERSION=$(git tag -l | tail -1)
 > conda activate dissect_${VERSION}

## Running

Fill in the dissectBCL.ini file appropriately. By default the config file is expected to be in ~/configs/dissectBCL_prod.ini.

 > dissect

or 

 > dissect -c /path/to/config.ini

or

 > dissect -f /path/to/flowcell.ini

## Docs

Documentation is available [here](https://dissectbcl.readthedocs.io/en/latest/).
