# Contributing

## Issues

Please double check if any similar issue already exists before opening one.

## Pull requests

There are 2 eternal branches:

 - main
 - prod

No development happens on the prod branch and the prod branch should never be merged back into main.
Releases should be tagged from the prod branch, after newly integrated feature(s) are considered stable.

Feature branches are merged into main, and require 3 passing status checks:

 - docs build passing  
 - flake8 passing
 - pytests passing

Note that flake8 tests are done via a pre-commit test as well.
Note that pytests are still under development and prone to change.

