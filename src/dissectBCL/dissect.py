from dissectBCL import fakeNews, misc
from dissectBCL.classes import sampleSheetClass, flowCellClass
from dissectBCL.demux import prepConvert, demux
from dissectBCL.postmux import postmux
from dissectBCL.drHouse import initClass
from dissectBCL.fakeNews import mailHome
from rich import print, inspect
import os
import logging
from pathlib import Path
import rich_click as click
from importlib.metadata import version
from time import sleep


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
def dissect(configfile):
    '''
    define config file and start main dissect function.
    '''
    print("This is dissectBCL version {}".format(version("dissectBCL")))
    print("Loading conf from {}".format(configfile))
    config = misc.getConf(configfile)
    main(config)


def main(config):
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
        flowcellName, flowcellDir = misc.getNewFlowCell(config)
        if flowcellName:
            # set exit stats
            exitStats = {}

            # Define a logfile.
            logFile = os.path.join(
                config['Dirs']['flowLogDir'],
                flowcellName + '.log'
            )

            # initiate log
            logging.basicConfig(
                filename=logFile,
                level="DEBUG",
                format="%(levelname)s    %(asctime)s    %(message)s",
                filemode='a',
                force=True
            )

            logging.info("Log Initiated - flowcell:{}, filename:{}".format(
                flowcellName,
                logFile
            ))

            # Create classes.
            flowcell = flowCellClass(
                name=flowcellName,
                bclPath=flowcellDir,
                inBaseDir=config['Dirs']['baseDir'],
                outBaseDir=config['Dirs']['outputDir'],
                origSS=os.path.join(flowcellDir, 'SampleSheet.csv'),
                runInfo=os.path.join(flowcellDir, 'RunInfo.xml'),
                logFile=logFile,
                config=config
            )
            inspect(flowcell)

            # Parse sampleSheet information.
            sampleSheet = sampleSheetClass(
                flowcell.origSS,
                flowcell.lanes,
                config
            )
            inspect(sampleSheet)
            exitStats['premux'] = prepConvert(
                flowcell,
                sampleSheet,
            )

            # Start demultiplexing.
            exitStats['demux'] = demux(
                sampleSheet,
                flowcell,
                config
            )
            inspect(sampleSheet)

            # postmux
            exitStats['postmux'] = postmux(
                flowcell,
                sampleSheet,
                config
            )

            # transfer data
            for outLane in sampleSheet.ssDic:
                # Copy over files.
                transferTime, shipDic = fakeNews.shipFiles(
                    os.path.join(
                        flowcell.outBaseDir,
                        outLane
                    ),
                    config
                )
                exitStats[outLane] = shipDic
                # Push stats to parkour.
                exitStats[outLane]['pushParkour'] = fakeNews.pushParkour(
                    flowcell.flowcellID,
                    sampleSheet,
                    config,
                    flowcell.bclPath
                )
                # Create diagnosis + parse QC stats
                drHouse = initClass(
                    os.path.join(
                        flowcell.outBaseDir,
                        outLane
                    ),
                    flowcell.startTime,
                    sampleSheet.flowcell,
                    sampleSheet.ssDic[outLane],
                    transferTime,
                    exitStats,
                    config['Dirs']['baseDir']
                    )
                inspect(drHouse)
                # Create email.
                subject, _html = drHouse.prepMail()
                # Send it.
                mailHome(subject, _html, config)
                Path(
                        os.path.join(
                            flowcell.outBaseDir,
                            outLane,
                            'communication.done'
                        )
                    ).touch()
            # Fix logs.
            fakeNews.organiseLogs(flowcell, sampleSheet)
        else:
            print("No flowcells found. Go back to sleep.")
            sleep(60*60)
