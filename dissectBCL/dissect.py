from dissectBCL import fakeNews, misc
from dissectBCL.fakeNews import pullParkour
from dissectBCL.classes import sampleSheetClass, flowCellClass
from dissectBCL.preFQ import prepConvert, demux
from rich import print, inspect
import logging
import os

def main():
    # Read config
    config = misc.getConf()
    # Search flowcell and initiate flowcell Class if found.
    flowcellName, flowcellDir = misc.getNewFlowCell(config)
    if flowcellName:
        # Start the logs.
        logFile = os.path.join(config['Dirs']['flowLogDir'], flowcellName + '.log')
        fakeNews.setLog(logFile)

        # Create classes.
        flowcell = flowCellClass(
            name = flowcellName,
            bclPath = flowcellDir,
            inBaseDir = config['Dirs']['baseDir'],
            outBaseDir = config['Dirs']['outputDir'],
            origSS = os.path.join(flowcellDir, 'SampleSheet.csv'),
            runInfo = os.path.join(flowcellDir, 'RunInfo.xml'),
            logFile = logFile
        )
        # Parse sampleSheet information.
        sampleSheet = sampleSheetClass( flowcell.origSS, flowcell.lanes, config )
        sampleSheet = prepConvert(flowcell, sampleSheet)
        # Start demultiplexing.
        inspect(sampleSheet)
        demux(sampleSheet, config['Dirs']['outputDir'])
        # Kick off demultiplexing.
        # Note we for loop over the different lanes.

        #for outLane in sampleSheet.ssDic:
        #    outFolder = os.path.join(config['Dirs']['outputDir'])
        #    demux(sampleSheet.ssDic[outLane]['sampleSheet'], sampleSheet.ssDic[outLane]['mask'], outFolder)
        
        

    else:
        print("Nothing to do. Moving on.")
