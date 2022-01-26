from dissectBCL import fakeNews, misc
from dissectBCL.logger import setLog
from dissectBCL.classes import sampleSheetClass, flowCellClass
from dissectBCL.demux import prepConvert, demux
from dissectBCL.postmux import postmux
from dissectBCL.drHouse import initClass
from rich import print, inspect
import os


def main():
    # Read config
    config = misc.getConf()
    # Search flowcell and initiate flowcell Class if found.
    flowcellName, flowcellDir = misc.getNewFlowCell(config)
    if flowcellName:
        # Start the logs.
        logFile = os.path.join(
            config['Dirs']['flowLogDir'],
            flowcellName + '.log'
        )
        setLog(logFile)

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
        sampleSheet = prepConvert(flowcell, sampleSheet)
        # Start demultiplexing.
        sampleSheet = demux(sampleSheet, flowcell, config)
        inspect(sampleSheet)
        # postmux
        postmux(flowcell, sampleSheet, config)
        # QC
        fakeNews.runSeqReports(flowcell, sampleSheet, config)
        for outLane in sampleSheet.ssDic:
            # Copy over files.
            transferTime, shipDic = fakeNews.shipFiles(
                os.path.join(
                    flowcell.outBaseDir,
                    outLane
                ),
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
                shipDic)
            inspect(drHouse)
            # Create email.
            subject, _html = drHouse.prepMail()
            # Send it.
            fakeNews.mailHome(subject, _html, config)

    else:
        print("Nothing to do. Moving on.")
