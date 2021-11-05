import configparser
import os
import sys
import rich
import glob

#Define flowcell class, which will contain all necesary information
class flowCell:
    sequencer = ''
    lanes = 0
    ss = ''
    loc = ''
    def __init__(self, name):
        self.name = name 

def getConf():
    # Get userDir
    homeDir = os.path.expanduser("~")
    # Fetch ini file and stop when it's not there.
    confLoc = os.path.join(homeDir, 'dissectBCL.ini')
    if not os.path.exists(confLoc):
        sys.stderr.write("Ini file not found. Exiting.\n")
        sys.exit(1)
    else:
        # Read config and return
        config = configparser.ConfigParser()
        config.read( confLoc )
        return config

def getNewFlowCell(config):
    baseDir = config['Dirs']['baseDir']
    outDir = config['Dirs']['outputDir']
    # Define a dict that maps the 'illumina letters' to a sequencer.
    sequencers = {
        'A':'NovaSeq',
        'N':'NextSeq',
        'M':'MiSeq'
    }
    # Get directories that are done sequencing (RTAcomplete flag.)
    flowCells = glob.glob(
        os.path.join( baseDir, '*', 'RTAComplete.txt' )
        )
    # Check if the flowcell exists in the output directory.
    for flowcell in flowCells:
        flowcellName = flowcell.split('/')[-2]
        # Look for a folder containing the flowcellname.
        # An empty list is returned if no directory exists.
        if not glob.glob(
                os.path.join( outDir, flowcellName ) + "*"
            ):
            rich.print("Unprocessed flowcell found: [green]{}[/green]".format(flowcellName) )
            # Initiate flowcell class.
            unprocessedFlowcell = flowCell(flowcellName)
            # fill in the sequencer.
            unprocessedFlowcell.sequencer = sequencers[flowcellName.split('_')[1][0]]
            # fill in the location.
            unprocessedFlowcell.loc = flowcell
            # fill in the sampleSheet location.
            if not os.path.exists( os.path.join( flowcell, 'SampleSheet.csv' ) ):
                print( os.path.join( flowcell, 'SampleSheet.csv') )
                sys.stderr.write("SampleSheet not found for {}. Exiting.\n".format(flowcellName) )
                sys.exit(1)
            else:
                unprocessedFlowcell.ss = os.path.join( flowcell, 'SampleSheet.csv' )
            return unprocessedFlowcell


    

config = getConf()
a = getNewFlowCell(config)
rich.print( a() )



