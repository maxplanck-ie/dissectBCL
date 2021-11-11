import configparser
import os
import sys
import rich
import glob
import xml.etree.ElementTree as ET
import pandas as pd
from dissectBCL.flowCellClass import flowCell

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
    outBaseDir = config['Dirs']['outputDir']
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
            unprocessedFlowcell = flowCell(
                name = flowcellName,
                bclPath = flowcellDir,
                origSS = os.path.join(flowcellDir, 'SampleSheet.csv'),
                runInfo = os.path.join(flowcellDir, 'RunInfo.xml'),
                inBaseDir = baseDir,
                outBaseDir = outBaseDir,
            )
            # FileChecks:
            fileChecks, fileChecksDic = unprocessedFlowcell.filesExist()
            if not fileChecks:
                sys.stderr.write("Not all files exist. Exiting..")
                sys.stderr.write("{}".format(fileChecksDic))
                sys.exit(1)
            return unprocessedFlowcell
    return None


# Parse runInfo.xml
def parseRunInfo(runInfo):
    tree = ET.parse(runInfo)
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
            lanes = int(i.attrib['LaneCount'])
        if i.tag == 'Instrument':
            instrument = i.text
        if i.tag == 'Flowcell':
            flowcellID = i.text
    return {
        'readDic': readDic,
        'lanes': lanes,
        'instrument': instrument,
        'flowcellID': flowcellID
    }


# Parse sampleSheet
def parseSS(ss):
    """
    We read the sampleSheet csv, and remove the stuff above the header.
    """
    ssdf = pd.read_csv(ss, sep=',')
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
    """
    Do we need to split per lane ? 
    We like to split per lane because, less barcodes = bigger chance for allowing more mismatches.
    We can't split per lane if:
      - 1 sample is loaded on multiple lanes
    or
      - 1 project is loaded on multiple lanes
    """
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
    return {
        ssdf,
        laneSplitStatus
    }
