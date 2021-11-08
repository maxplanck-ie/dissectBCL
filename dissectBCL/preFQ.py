import configparser
import os
import sys
import rich
import glob
import xml.etree.ElementTree as ET
import pandas as pd


# Define flowcell class, which will contain all information
class flowCell:
    def __init__(self, name):
        self.name = name


# Define config reader.
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
        config.read(confLoc)
        return config


# search for new flowcells.
def getNewFlowCell(config):
    baseDir = config['Dirs']['baseDir']
    outDir = config['Dirs']['outputDir']
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
                os.path.join(outDir, flowcellName) + "*"
        ):
            rich.print(
                "Unprocessed flowcell found: \
                    [green]{}[/green]".format(flowcellName))
            # Initiate flowcell class.
            unprocessedFlowcell = flowCell(flowcellName)
            # fill in the sequencer.
            unprocessedFlowcell.sequencer = \
                sequencers[flowcellName.split('_')[1][0]]
            # fill in the location.
            unprocessedFlowcell.bclPath = \
                flowcellDir
            # fill in the sampleSheet location.
            ssPath = os.path.join(flowcellDir, 'SampleSheet.csv')
            if not os.path.exists(ssPath):
                sys.stderr.write(
                    "SampleSheet missing: {}. Exiting.\n".format(flowcellName))
                sys.exit(1)
            else:
                unprocessedFlowcell.ssPath = ssPath
            # fill in the runInfo location.
            runinfoPath = os.path.join(flowcellDir, 'RunInfo.xml')
            if not os.path.exists(runinfoPath):
                sys.stderr.write(
                    "RunInfo.xml missing: {}. Exiting.\n".format(flowcellName))
                sys.exit(1)
            else:
                unprocessedFlowcell.runInfoPath = runinfoPath
            return unprocessedFlowcell
    return None


# Parse runInfo.xml
def parseRunInfo(flowClass):
    tree = ET.parse(flowClass.runInfoPath)
    root = tree.getroot()
    readDic = {}
    for i in root.iter():
        if i.tag == 'Read':
            if i.attrib['IsIndexedRead'] == 'Y':
                readType = 'Index'
            else:
                readType = 'Read'
            readKey = 'Read' + i.attrib['Number']
            readDic[readKey] = [i.attrib['NumCycles'], readType]
        if i.tag == 'FlowcellLayout':
            flowClass.lanes = int(i.attrib['LaneCount'])
        if i.tag == 'Instrument':
            flowClass.instrument = i.text
        if i.tag == 'Flowcell':
            flowClass.flowcellID = i.text
    flowClass.readLens = readDic
    return flowClass


# Parse sampleSheet
def parseSS(flowClass):
    ssdf = pd.read_csv(flowClass.ssPath, sep=',')
    # There is a bunch of header 'junk' that we don't want.
    # subset the df from [data] onwards.
    startIx = ssdf[ssdf.iloc[:, 0] == '[Data]'].index.values[0] + 1
    # only take those with actual sample information.
    ssdf = ssdf.iloc[startIx:, :]
    ssdf.columns = ssdf.iloc[0]
    ssdf = ssdf.drop(ssdf.index[0])
    # Reset index
    ssdf.reset_index(inplace=True)
    ssdf.index.name = None
    # Remove 'level0' column
    ssdf.drop('level_0', axis=1, inplace=True)
    ssdf = ssdf.dropna(axis=1, how='all')
    flowClass.sampleSheet = ssdf
    laneSplitStatus = True
    # Do we need lane splitting or not ?
    # If there is at least one sample in more then 1 lane, we cannot split:
    if sum(ssdf['Sample_Name'].value_counts() > 1) > 0:
        laneSplitStatus = False
    # If one project is split over multiple lanes, we also don't split:
    projects = list(ssdf['Sample_Project'].unique())
    for project in projects:
        if len(
            list(
                ssdf[ssdf['Sample_Project'] == project]['Lane'].unique()
                )) > 1:
            laneSplitStatus = False
    flowClass.laneSplitStatus = laneSplitStatus
    return flowClass
