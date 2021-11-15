import configparser
import os
import sys
import rich
import glob
import xml.etree.ElementTree as ET
import pandas as pd
from dissectBCL.flowCellClass import flowCell
import logging


def dirCreator(flowClass):
    """
    Fills flowClass.inferredVars['outDir'] with a dictionary ({path:lane#, ...}) if lanesplitting
    else it's the absolute path (string)
    """
    # If we split lanes, we need 1 directory / lane.
    if flowClass.inferredVars['laneSplitStatus']:
        outLaneDic = {}
        for lane in range(flowClass.inferredVars['lanes']):
            outPath = flowClass.inVars['name'] + '_lanes' + str(lane + 1)
            outAbsPath = os.path.join(flowClass.inVars['outBaseDir'], outPath)
            outLaneDic[outAbsPath] = lane + 1
            if not os.path.exists(outAbsPath):
                os.mkdir(outAbsPath)
                rich.print("{} created.".format(outAbsPath))
            else:
                rich.print("{} exists. Moving on.".format(outAbsPath))
    
    # if not, we only have 1 output directory
    else:
        # we don't split, just combine the lane numbers for consistency with previous pipeline
        lanesStr = '_'.join(map(str, list(range(1,flowClass.inferredVars['lanes']+1, 1)) ))
        outPath = flowClass.inVars['name'] + '_lanes' + lanesStr
        outAbsPath = os.path.join(flowClass.inVars['outBaseDir'], outPath)
        outLaneDic = outAbsPath
        if not os.path.exists(outAbsPath):
            os.mkdir(outAbsPath)
            rich.print("{} dir created.".format(outAbsPath))
        else:
            rich.print("{} exists. Moving on.".format(outAbsPath))

    flowClass.inferredVars['outDir'] = outLaneDic
    return flowClass


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
    logDir = config['Dirs']['logDir']
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
                logOut = os.path.join(logDir, flowcellName + ".out"),
                logErr = os.path.join(logDir, flowcellName + ".err")
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