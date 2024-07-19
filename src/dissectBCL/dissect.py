from dissectBCL.misc import getNewFlowCell
from dissectBCL.misc import getConf
from dissectBCL.flowcell import flowCellClass
from importlib.metadata import version
import logging
import os
from pathlib import Path
from rich import print, inspect
import rich_click as click
from time import sleep
import sys

@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"]
    )
)
@click.option(
   "-c",
   "--configfile",
   type=click.Path(exists=True),
   required=False,
   default=os.path.expanduser('~/configs/dissectBCL_prod.ini'),
   help='specify a custom ini file.',
   show_default=True
)
@click.option(
    "-f",
    "--flowcellpath",
    required=False,
    default=None,
    help='specify a full path to a flow cell to process. Should be pointing to a directory written by an Illumina sequencer'
)
def dissect(configfile, flowcellpath):
    '''
    define config file and start main dissect function.
    '''
    print("This is dissectBCL version {}".format(version("dissectBCL")))
    print("Loading conf from {}".format(configfile))
    config = getConf(configfile)
    main(config, flowcellpath)


def main(config, flowcellpath):
    '''
    every hour checks for a new flow cell.
    if new flowcell:
        - initiate log
        - create flowcellClass
        - create sampleSheetClass
        - prepconvert, demux, postmux
        - QC & communication.
    '''

    # Set pipeline.
    while True:
        # Reload setlog
        flowcellName, flowcellDir = getNewFlowCell(config, flowcellpath)
        if flowcellName:

            # Define a logfile.
            logFile = Path(config['Dirs']['flowLogDir'], flowcellName + '.log')

            # initiate log
            logging.basicConfig(
                filename=logFile,
                level="DEBUG",
                format="%(levelname)s    %(asctime)s    %(message)s",
                filemode='a',
                force=True
            )
            # Set flowcellname in log.
            logging.info(f"Log Initiated - flowcell:{flowcellName}, filename:{logFile}")
            print(f"Logfile set as {logFile}")
            # Include dissectBCL version in log
            logging.info(f"dissectBCL - version {version("dissectBCL")}")
            # Include software versions in log
            for lib in config['softwareVers']:
                logging.debug(f"{lib} = {config['softwareVers'][lib]}")

            # Create class.
            flowcell = flowCellClass(name=flowcellName, bclPath=flowcellDir, logFile=logFile, config=config)
            # Run workflow
            flowcell.prepConvert()
            flowcell.demux()
            flowcell.postmux()
            flowcell.fakenews()
            flowcell.organiseLogs()
            inspect(flowcell)
        else:
            print("No flowcells found. Go back to sleep.")
            sleep(60*60)

def createFlowcell(config, fpath, logFile = None):
    config = getConf(config)
    flowcellName, flowcellDir = getNewFlowCell(config, fpath)
    if not logFile:
        logging.basicConfig(
            stream=sys.stdout,
            level="DEBUG",
            format="%(levelname)s    %(asctime)s    %(message)s",
            filemode='a',
            force=True
        )
        logFile = 'STDOUT'
    else:
        logging.basicConfig(
            filename=logFile,
            level="DEBUG",
            format="%(levelname)s    %(asctime)s    %(message)s",
            filemode='a',
            force=True
        )
    return flowCellClass(name=flowcellName, bclPath=flowcellDir, logFile=logFile, config=config)