from dissectBCL import fakeNews, misc
from dissectBCL.classes import sampleSheetClass, flowCellClass
from dissectBCL.demux import prepConvert, demux
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
        fakeNews.setLog(logFile)

        # Create classes.
        flowcell = flowCellClass(
            name=flowcellName,
            bclPath=flowcellDir,
            inBaseDir=config['Dirs']['baseDir'],
            outBaseDir=config['Dirs']['outputDir'],
            origSS=os.path.join(flowcellDir, 'SampleSheet.csv'),
            runInfo=os.path.join(flowcellDir, 'RunInfo.xml'),
            logFile=logFile
        )
        inspect(flowcell)
        # Parse sampleSheet information.
        sampleSheet = sampleSheetClass(
            flowcell.origSS,
            flowcell.lanes,
            config
        )
        sampleSheet = prepConvert(flowcell, sampleSheet)
        inspect(sampleSheet)
        # Start demultiplexing.
        demux(sampleSheet, flowcell, config)

    else:
        print("Nothing to do. Moving on.")
