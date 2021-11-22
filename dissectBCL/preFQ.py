import os
import sys
import rich
import glob
from dissectBCL.classes import flowCellClass
import logging


# search for new flowcells.
def getNewFlowCell(config):
    baseDir = config['Dirs']['baseDir']
    outBaseDir = config['Dirs']['outputDir']
    flowLogDir = config['Dirs']['flowLogDir']

    # Define a dict that maps the 'illumina letters' to a sequencer.
    sequencers = {
        'A': 'NovaSeq',
        'N': 'NextSeq',
        'M': 'MiSeq'
    }
    # Get directories that are done sequencing (RTAcomplete flag.)
    flowCells = glob.glob(
        os.path.join(baseDir, '*', 'RTAComplete.txt')
        )
    # Check if the flowcell exists in the output directory.
    for flowcell in flowCells:
        flowcellName = flowcell.split('/')[-2]
        flowcellDir = flowcell.replace("/RTAComplete.txt", "")
        # Look for a folder containing the flowcellname.
        # An empty list is returned if no directory exists.
        if not glob.glob(
                os.path.join(outBaseDir, flowcellName) + "*"
        ):
            rich.print(
                "Unprocessed flowcell found: \
                    [green]{}[/green]".format(flowcellName))

            # Initiate flowcellClass
            unprocessedFlowcell = flowCellClass(
                name = flowcellName,
                bclPath = flowcellDir,
                origSS = os.path.join(flowcellDir, 'SampleSheet.csv'),
                runInfo = os.path.join(flowcellDir, 'RunInfo.xml'),
                inBaseDir = baseDir,
                outBaseDir = outBaseDir,
                logFile = os.path.join(flowLogDir, flowcellName + ".out")
            )
            return unprocessedFlowcell
    return None