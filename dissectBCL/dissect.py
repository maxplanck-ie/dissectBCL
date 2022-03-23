from dissectBCL import fakeNews, misc
from dissectBCL.logger import setLog
from dissectBCL.classes import sampleSheetClass, flowCellClass
from dissectBCL.demux import prepConvert, demux
from dissectBCL.postmux import postmux
from dissectBCL.drHouse import initClass
from rich import print, inspect
import os
import signal
from threading import Event
from pathlib import Path


def main():
    # Start pipeline.
    while True:
        # Set up sleeper
        HUP = Event()

        def breakSleep(signo, _frame):
            HUP.set()

        def sleep():
            HUP.wait(timeout=float(60*60))
        signal.signal(signal.SIGHUP, breakSleep)

        # Read config
        config = misc.getConf()
        flowcellName, flowcellDir = misc.getNewFlowCell(config)
        if flowcellName:
            # set exit stats
            exitStats = {}

            # Start the logs.
            logFile = os.path.join(
                config['Dirs']['flowLogDir'],
                flowcellName + '.log'
            )
            exitStats['log'] = setLog(logFile)

            # Create classes.
            flowcell = flowCellClass(
                name=flowcellName,
                bclPath=flowcellDir,
                inBaseDir=config['Dirs']['baseDir'],
                outBaseDir=config['Dirs']['outputDir'],
                origSS=os.path.join(flowcellDir, 'SampleSheet.csv'),
                runInfo=os.path.join(flowcellDir, 'RunInfo.xml'),
                logFile=logFile,
            )
            inspect(flowcell)

            # Parse sampleSheet information.
            sampleSheet = sampleSheetClass(
                flowcell.origSS,
                flowcell.lanes,
                config
            )
            exitStats['premux'] = prepConvert(
                flowcell,
                sampleSheet
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
                    config
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
                fakeNews.mailHome(subject, _html, config)
                Path(
                        os.path.join(
                            flowcell.outBaseDir,
                            outLane,
                            'communication.done'
                        )
                    ).touch()
                print()
            # Fix logs.
            fakeNews.organiseLogs(flowcell, sampleSheet)
        else:
            print("No flowcells found. I go back to sleep.")
            sleep()
            continue
