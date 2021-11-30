import os
import sys
import configparser
import xml.etree.ElementTree as ET
import glob
from dissectBCL.fakeNews import log
import pandas as pd
import numpy as np


# Define config reader.
def getConf():
    # Get userDir
    homeDir = os.path.expanduser("~")
    # Fetch ini file and stop when it's not there.
    confLoc = os.path.join(homeDir, 'dissectBCL.ini')
    if not os.path.exists(confLoc):
        log.critical(
            "[red]Ini file not found. (~/dissectBCL.ini) Exiting..[/red]"
        )
        sys.exit(1)
    else:
        # Read config and return
        config = configparser.ConfigParser()
        config.read(confLoc)
        return config


# Find new flowcells.
def getNewFlowCell(config):
    # set some config vars.
    baseDir = config['Dirs']['baseDir']
    outBaseDir = config['Dirs']['outputDir']
    # Glob over the bcl directory to get all flowcells.
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
            return flowcellName, flowcellDir
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


def hamming(s1, s2):
    # We have some basket cases (multimodal)
    # Where barcode is nan (type as float)
    # Ignore these for now.
    if type(s1) == float or type(s2) == float:
        return 0
    dist = 0
    for step in range(len(s1)):
        if s1[step] != s2[step]:
            dist += 1
    return dist


def joinLis(lis, joinStr=""):
    """
    join a list into a string (without spaces).
    elements are converted to strings.
    """
    return joinStr.join([str(i) for i in lis])


def lenMask(recipe, minl):
    """
    take length of recipe (runInfo) and length of a barcode and return a mask.
    e.g. 8bp index, 10bp sequenced, returns I8N2
    """
    if recipe-minl > 0:
        return "I{}N{}".format(minl, recipe-minl)
    else:
        return "I{}".format(minl)


def bclConvPipeLogger(PIPE):
    for line in iter(PIPE.readline()):
        log.debug('BCLConvert: {}'.format(line))


def P5Seriesret(df):
    if 'index2' in list(df.columns):
        return df['index2']
    else:
        return pd.Series()
