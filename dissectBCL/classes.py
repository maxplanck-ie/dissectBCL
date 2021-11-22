import os
import xml.etree.ElementTree as ET
from fakeNews import logger
import sys
from rich import pretty

# Define flowcell class, which will contain all information
class flowCellClass:
    """This is a flowCell class, which contains:
    - prior information set in the config file and read from the flowcell directory (inVars).
    - inferred variables from the RunInfo.xml and the sampleSheet (inferredVars).
    """

    # fileChecks
    def filesExist(self):
        """
        Check if paths exist:
          - Flowcell directory with BCL files
          - the original sampleSheet
          - RunInfo.xml
          - base directory for bcl files.
          - directory where to write output.
        """
        for f in [
            self.bclPath,
            self.origSS,
            self.runInfo,
            self.inBaseDir,
            self.outBaseDir
        ]:
            if not os.path.exists(f):
                logger.critical("{} doesn't exist. Exiting".format(f) )
                sys.exit(1)

    
    # Parse runInfo
    def parseRunInfo(self):
        """
        Takes the path to runInfo.xml and parses it.
        Returns:
         - a read dictionary (containing length for each read)
         - number of lanes (int)
         - the instrument (str)
         - the flowcellID (str)
        """

        tree = ET.parse(self.runInfo)
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
        return readDic, lanes, instrument, flowcellID


    def __init__(self, name, bclPath, origSS, runInfo, inBaseDir, outBaseDir,logFile):
        sequencers = {
            'A': 'NovaSeq',
            'N': 'NextSeq',
            'M': 'MiSeq'
        }

        self.name = name
        self.sequencer = sequencers[name.split('_')[1][0]]
        self.bclPath = bclPath
        self.origSS = origSS
        self.runInfo = runInfo
        self.inBaseDir = inBaseDir
        self.outBaseDir = outBaseDir
        self.logFile = logFile
        self.readDic, self.lanes, self.instrument, self.flowcellID = self.parseRunInfo()
        logger.info("flowcell Class for {} initiated".format(self.name) )
        logger.info("flowcell Class for {} initiated".format(self.name) )

        # Check if important files exists("")
        self.filesExist()


class sampleSheetClass:
    """The sampleSheet class.
    - Contains the pandas object.
    - Many functions related to the sampleSheet (indexes, laneSplitting etc.)
    """

    def decideSplit(self):
        """
        Do we need to split per lane ? 
        We like to split per lane because, less barcodes = bigger chance for allowing more mismatches.
        We can't split per lane if:
        - 1 sample is loaded on multiple lanes
        or
        - 1 project is loaded on multiple lanes
        or
        - there are more then 1 lanes, but only 1 is specified in sampleSheet
        """
        # Do we need lane splitting or not ?
        # If there is at least one sample in more then 1 lane, we cannot split:
        if sum(self.ssdf['Sample_Name'].value_counts() > 1) > 0:
            return False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.ssdf['Sample_Project'].unique())
        for project in projects:
            if len(
                list(
                    self.ssdf[self.ssdf['Sample_Project'] == project]['Lane'].unique()
                    )) > 1:
                return False
        # Sometimes only 1 lane is listed, although there are multiple (so here we also don't split)
        if len(list(self.ssdf['Lane'].unique())) < self.runInfoLanes:
            return False
        return True


    def __init__(self, sampleSheet, lanes):
        self.ssdf = sampleSheet
        self.runInfoLanes = lanes
        self.laneSplitStatus = self.decideSplit()
        # If we split, we want to