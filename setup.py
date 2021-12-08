import setuptools

setuptools.setup(
    name='dissectBCL',
    version="0.0.1",
    author='WardDeb',
    author_email='deboutte@ie-freiburg.mpg.de',
    description="MPI-IE's demux pipeline."
                "Continuation of 'TheWhoTheWhatTheHuh' from DpRyan. ",
    scripts=['bin/dissect'],
    packages=["dissectBCL"],
    python_requires=">3.7",
    install_requires=[
        'flake8==3.9.2',
        'matplotlib==3.4.3',
        'pandas==1.3.4',
        'rich==10.12.0',
        'requests==2.26.0',
        'coverage==6.2',
        'pytest==6.2.5',
        'tabulate==0.8.9',
        'pandoc',
        'pypandoc'
    ]
)
