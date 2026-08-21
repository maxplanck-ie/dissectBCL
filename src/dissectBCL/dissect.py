import logging
import sys
from importlib.metadata import version
from pathlib import Path
from time import sleep

import rich_click as click
from rich import print

from dissectBCL.flowcell import flowCellClass
from dissectBCL.misc import getConf, getNewFlowCell


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True),
    required=False,
    default=Path("~/configs/dissectBCL_prod.ini").expanduser(),
    help="specify a custom ini file.",
    show_default=True,
)
@click.option(
    "-f",
    "--flowcellpath",
    required=False,
    default=None,
    help="specify a full path to a flow cell to process. Should be pointing to a directory written by an Illumina sequencer",
)
@click.option(
    "-s",
    "--sequencer",
    default=None,
    type=click.Choice(["illumina", "aviti"], case_sensitive=True),
    help="Restrict the run to one platform ('illumina' or 'aviti'): only that platform's "
    "config keys are read and only its flowcells are watched. Required when used together "
    "with -f/--flowcellpath. Omit to watch both platforms from one config, as before.",
)
@click.option(
    "-F",
    "--forcelanesplit",
    is_flag=True,
    default=False,
    help="Force lane splitting even if specified in the sample sheet.",
)
def dissect(configfile, flowcellpath, sequencer, forcelanesplit):
    """
    define config file and start main dissect function.
    """
    print(f"This is dissectBCL version {version('dissectBCL')}")
    print(f"Loading conf from {configfile}")
    config = getConf(configfile, sequencer=sequencer)
    main(config, flowcellpath, sequencer, forcelanesplit)


def main(config, flowcellpath, platformFilter, forcelanesplit):
    """
    every hour checks for a new flow cell.
    if new flowcell:
        - initiate log
        - create flowcellClass
        - create sampleSheetClass
        - prepconvert, demux, postmux
        - QC & communication.

    platformFilter, when set, restricts every poll to that one platform's
    config keys ('illumina' or 'aviti'). It must stay fixed across loop
    iterations - it must not be overwritten by getNewFlowCell's per-call
    return value, which is None on a no-match poll.
    """

    # Set pipeline.
    while True:
        # Reload setlog
        flowcellName, flowcellDir, sequencer = getNewFlowCell(
            config, flowcellpath, platformFilter
        )

        if flowcellName:
            # Define a logfile.
            logFile = Path(
                config["Dirs"][f"flowLogDir_{sequencer}"], flowcellName + ".log"
            )
            logFile.parent.mkdir(parents=True, exist_ok=True)

            # initiate log
            logging.basicConfig(
                filename=logFile,
                level="DEBUG",
                format="%(levelname)s    %(asctime)s    %(message)s",
                filemode="a",
                force=True,
            )

            # Include log to stdout if debug mode is on
            if config["communication"]["debug_mode"]:
                # Add console handler
                console = logging.StreamHandler(sys.stdout)
                console.setLevel(logging.DEBUG)
                console.setFormatter(
                    logging.Formatter("%(levelname)s    %(asctime)s    %(message)s")
                )
                logging.getLogger().addHandler(console)

            # Set flowcellname in log.
            logging.info(
                f"Log Initiated - flowcell:{flowcellName}, filename:{logFile}, sequencer:{sequencer}"
            )

            print(f"Logfile set as {logFile}")
            # Include dissectBCL version in log
            logging.info(f"dissectBCL - version {version('dissectBCL')}")
            # Include software versions in log
            for lib in config["softwareVers"]:
                logging.debug(f"{lib} = {config['softwareVers'][lib]}")

            # Create class.
            flowcell = flowCellClass(
                name=flowcellName,
                bclPath=flowcellDir,
                logFile=logFile,
                config=config,
                sequencer=sequencer,
                forceLaneSplit=forcelanesplit,
            )
            flowcell.prepConvert()
            if sequencer == "illumina":
                # flowcell.prepConvert()
                flowcell.demux()
            else:
                flowcell.demux_aviti()
            flowcell.postmux()
            flowcell.fakenews()
            flowcell.organiseLogs()
        else:
            print("Going back to sleep for 60 minutes.")
            sleep(60 * 60)


def createFlowcell(config, fpath, sequencer, logFile=None, forceLaneSplit=False):
    config = getConf(config, sequencer=sequencer)
    flowcellName, flowcellDir, sequencer = getNewFlowCell(config, fpath, sequencer)
    if not logFile:
        logging.basicConfig(
            stream=sys.stdout,
            level="DEBUG",
            format="%(levelname)s    %(asctime)s    %(message)s",
            filemode="a",
            force=True,
        )
        logFile = "STDOUT"
    else:
        Path(logFile).parent.mkdir(parents=True, exist_ok=True)
        logging.basicConfig(
            filename=logFile,
            level="DEBUG",
            format="%(levelname)s    %(asctime)s    %(message)s",
            filemode="a",
            force=True,
        )
    return flowCellClass(
        name=flowcellName,
        bclPath=flowcellDir,
        logFile=logFile,
        config=config,
        sequencer=sequencer,
        forceLaneSplit=forceLaneSplit,
    )
