import os
import xml.etree.ElementTree as ET
import sys
from rich import print
from dissectBCL.fakeNews import log
import pandas as pd


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
                print("[red]Not all necessary files found. Exiting..[/red]")
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
        readLens = []
        for i in root.iter():
            if i.tag == 'Read':
                if i.attrib['IsIndexedRead'] == 'Y':
                    readLens.append([int(i.attrib['NumCycles']), 'Index'])
                else:
                    readLens.append([int(i.attrib['NumCycles']), 'Read'])
            if i.tag == 'FlowcellLayout':
                lanes = int(i.attrib['LaneCount'])
            if i.tag == 'Instrument':
                instrument = i.text
            if i.tag == 'Flowcell':
                flowcellID = i.text
        return readLens, lanes, instrument, flowcellID


    def __init__(self, name, bclPath, origSS, runInfo, inBaseDir, outBaseDir,logFile):
        sequencers = {
            'A': 'NovaSeq',
            'N': 'NextSeq',
            'M': 'MiSeq'
        }
        log.warning("Initiating flowcellClass {}".format(name))
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
        self.readLens, self.lanes, self.instrument, self.flowcellID = self.parseRunInfo()
        


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
        log.info("Deciding lanesplit.")
        # Do we need lane splitting or not ?
        # If there is at least one sample in more then 1 lane, we cannot split:
        if sum(self.fullSS['Sample_Name'].value_counts() > 1) > 0:
            log.info("No lane splitting due to at least 1 sample in multiple lanes")
            return False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.fullSS['Sample_Project'].unique())
        for project in projects:
            if len(
                list(
                    self.fullSS[self.fullSS['Sample_Project'] == project]['Lane'].unique()
                    )) > 1:
                    log.info("No lane splitting due to 1 project over multiple lanes")
                    return False
        # Sometimes only 1 lane is listed, although there are multiple (so here we also don't split)
        if len(list(self.fullSS['Lane'].unique())) < self.runInfoLanes:
            log.info("No lane splitting due to 1 lane listed, {} found.".format(self.runInfoLanes))
            return False
        log.info("Splitting up lanes.")
        return True

    # Parse sampleSheet
    def parseSS(self):
        """
        We read the sampleSheet csv, and remove the stuff above the header.
        """
        log.info("Reading sampleSheet.")
        ssdf = pd.read_csv(self.origSs, sep=',')
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
        ssdf = ssdf.astype({'Lane': 'int32'})
        self.fullSS = ssdf
        self.laneSplitStatus = self.decideSplit()
        ssDic = {}
        # If we need to split per lane, return a dictionary with ss_lane_X as key, the dataframe as value.
        if self.laneSplitStatus:
            for lane in range(1,self.runInfoLanes + 1,1):
                key = self.flowcell + '_lanes_' + str(lane)
                ssDic[key] = {'sampleSheet':ssdf[ssdf['Lane'] == lane], 'lane' : lane}
            del self.fullSS
            return ssDic
        else:
            laneStr = '_'.join([lane for lane in range(1,self.runInfoLanes +1,1)])
            ssDic[laneStr] == {'sampleSheet':ssdf, 'lane':'all'}
        del self.fullSS
        return ssDic

    def __init__(self, sampleSheet, lanes):
        log.warning("initiating sampleSheetClass")
        self.runInfoLanes = lanes
        self.origSs = sampleSheet
        self.flowcell = sampleSheet.split('/')[-2]
        # ParseSS reads up and cleans the dataframe and will decide to split yes or no.
        # Object returned is a dictionary with dic[output folder] = {'sampleSheet':pandas df, 'lane':'all' or lane number}
        self.ssDic = self.parseSS()
        