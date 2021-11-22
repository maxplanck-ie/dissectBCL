import os
import xml.etree.ElementTree as ET
import sys
from rich import pretty
from dissectBCL.fakeNews import log


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
            log.info("Checking {}".format(f))
            if not os.path.exists(f):
                log.critical("{} doesn't exist. Exiting".format(f) )
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
        log.info("Parsing RunInfo.xml")
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
        log.info("Init flowcellClass {}".format(name))
        self.name = name
        self.sequencer = sequencers[name.split('_')[1][0]]
        self.bclPath = bclPath
        self.origSS = origSS
        self.runInfo = runInfo
        self.inBaseDir = inBaseDir
        self.outBaseDir = outBaseDir
        self.logFile = logFile
        # Run filesChecks
        self.filesExist()
        # populate runInfo vars.
        self.readDic, self.lanes, self.instrument, self.flowcellID = self.parseRunInfo()
        


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
        log.warning("Deciding lanesplit.")
        # Do we need lane splitting or not ?
        # If there is at least one sample in more then 1 lane, we cannot split:
        if sum(self.ssdf['Sample_Name'].value_counts() > 1) > 0:
            log.info("No lane splitting due to at least 1 sample in multiple lanes")
            return False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.ssdf['Sample_Project'].unique())
        for project in projects:
            if len(
                list(
                    self.ssdf[self.ssdf['Sample_Project'] == project]['Lane'].unique()
                    )) > 1:
                    log.info("No lane splitting due to 1 project over multiple lanes")
                return False
        # Sometimes only 1 lane is listed, although there are multiple (so here we also don't split)
        if len(list(self.ssdf['Lane'].unique())) < self.runInfoLanes:
            log.info("No lane splitting due to 1 lane listed, {} found.".format(self.runInfoLanes))
            return False
        log.info("Splitting up lanes.")
        return True

    # Parse sampleSheet
    def parseSS(ssPath):
        """
        We read the sampleSheet csv, and remove the stuff above the header.
        """
        log.warning("Reading sampleSheet.")
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


    def __init__(self, sampleSheet, lanes):
        self.ssdf = sampleSheet
        self.runInfoLanes = lanes
        self.laneSplitStatus = self.decideSplit()
        # If we split, we want to