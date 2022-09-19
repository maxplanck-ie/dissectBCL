# dissectBCL

Demultiplexing pipeline for illumina data (novaseq/miseq/nextseq). Continuation of Devon Ryan's [TWTWTWTW](https://github.com/maxplanck-ie/TheWhoTheWhatTheHuh).

## Installation.

Clone this repository, create the environment and pip install

 > git clone git@github.com:maxplanck-ie/dissectBCL.git  
 > cd dissectBCL  
 > conda create -f env.yml  
 > conda activate dissectBCL  
 > pip install ./  

Fill in the dissectBCL.ini file appropriately. By default the config file is expected to be in ~/configs/dissectBCL_prod.ini.

## Running.

 > dissect

or 

 > dissect -c /path/to/config.ini
