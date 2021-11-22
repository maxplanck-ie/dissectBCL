import os
import rich
import sys
import configparser
import xml.etree.ElementTree as ET
import pandas as pd


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
def parseSS(ssPath):
    """
    We read the sampleSheet csv, and remove the stuff above the header.
    """
    ssdf = pd.read_csv(ssPath, sep=',')
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
    return ssdf