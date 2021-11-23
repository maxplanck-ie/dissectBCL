from dissectBCL import preFQ, FQ, fakeNews, misc
from dissectBCL.classes import sampleSheetClass, flowCellClass
from rich import print, inspect
import pandas as pd
import argparse
import logging
import os

def main():
    # Read config
    config = misc.getConf()
    # Search flowcell and initiate flowcell Class if found.
    flowcellName, flowcellDir = misc.getNewFlowCell(config)
    if flowcellName:
        # Start the logs.
        logFile = os.path.join(config['Dirs']['flowLogDir'], flowcellName)
        fakeNews.setLog(logFile)

        # Create flowcell class.
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
        sampleSheet = sampleSheetClass( flowcell.origSS, flowcell.lanes )
        

        sampleSheet = preFQ.prepConvert(flowcell, sampleSheet)
        inspect(sampleSheet)

        #################### FQ MODULE ####################
        print("FQ MODULE")
    else:
        print("Nothing to do. Moving on.")
