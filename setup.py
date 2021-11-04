import setuptools

setuptools.setup(
    name= 'dissectBCL',
    version = "0.0.1",
    author = 'WardDeb',
    author_email = 'deboutte@ie-freiburg.mpg.de',
    description = "MPI-IE's demux pipeline. Continuation of 'TheWhoTheWhatTheHuh' from DpRyan. ",
    package_dir = {"": "dissectBCL"},
    python_requires = "==3.9.2",
    install_requires = [
        'rich==10.12.0',
        'matplotlib==3.4.3',
        'flake8==3.9.2'
    ]
)
