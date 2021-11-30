import os
import xml.etree.ElementTree as ET
import sys
from rich import print
from dissectBCL.fakeNews import log, pullParkour
import pandas as pd


class flowCellClass:
    """This is a flowCell class, which contains:
    - prior information set in the config file.
    - inferred variables from the RunInfo.xml.
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
                log.critical("{} doesn't exist. Exiting".format(f))
                print("[red]Not all necessary files found. Exiting..[/red]")
                sys.exit(1)

    # Parse runInfo
    def parseRunInfo(self):
        """
        Takes the path to runInfo.xml and parses it.
        Returns:
         - a sequencing recipe dictionary. {'Read1':100, 'Index1':10, ... }
         - number of lanes (int)
         - the instrument (str)
         - the flowcellID (str)
        """
        log.info("Parsing RunInfo.xml")
        tree = ET.parse(self.runInfo)
        root = tree.getroot()
        seqRecipe = {}
        readCount = 1
        indexCount = 1
        for i in root.iter():
            if i.tag == 'Read':
                if i.attrib['IsIndexedRead'] == 'Y':
                    ixStr = 'Index' + str(indexCount)
                    seqRecipe[ixStr] = ['I', int(i.attrib['NumCycles'])]
                    indexCount += 1
                else:
                    readStr = 'Read' + str(readCount)
                    seqRecipe[readStr] = ['Y', int(i.attrib['NumCycles'])]
                    readCount += 1
            if i.tag == 'FlowcellLayout':
                lanes = int(i.attrib['LaneCount'])
            if i.tag == 'Instrument':
                instrument = i.text
            if i.tag == 'Flowcell':
                flowcellID = i.text
        return seqRecipe, lanes, instrument, flowcellID

    def __init__(
        self,
        name,
        bclPath,
        origSS,
        runInfo,
        inBaseDir,
        outBaseDir,
        logFile
    ):
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
        self.seqRecipe, \
            self.lanes, \
            self.instrument, \
            self.flowcellID = self.parseRunInfo()


class sampleSheetClass:
    """The sampleSheet class.
    - Contains the pandas object.
    - Many functions related to the sampleSheet (indexes, laneSplitting etc.)
    - At initiation stage parkour is querried.
    """

    def decideSplit(self):
        """
        Do we need to split per lane ?
       Lane splitting = preffered:
       less barcodes = allow more mismatches.
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
            log.info(
                "No lane splitting: >= 1 sample in multiple lanes."
            )
            return False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.fullSS['Sample_Project'].unique())
        for project in projects:
            if len(
                list(self.fullSS[
                    self.fullSS['Sample_Project'] == project
                ]['Lane'].unique()
                )
            ) > 1:
                log.info(
                    "No lane splitting: >= 1 project in multiple lanes."
                )
                return False
        # Don't split if 1 lane in ss, multiple in runInfo
        if len(list(self.fullSS['Lane'].unique())) < self.runInfoLanes:
            log.info(
                "No lane splitting: 1 lane listed, {} found.".format(
                    self.runInfoLanes
                )
            )
            return False
        log.info("Splitting up lanes.")
        return True

    # Parse sampleSheet
    def parseSS(self, parkourDF):
        """
        We read the sampleSheet csv, and remove the stuff above the header.
        """
        log.info("Reading sampleSheet.")
        ssdf = pd.read_csv(self.origSs, sep=',')
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
        # If lanesplit: ret dict w/ ss_lane_X:df
        if self.laneSplitStatus:
            for lane in range(1, self.runInfoLanes + 1, 1):
                key = self.flowcell + '_lanes_' + str(lane)
                # if we have a parkour dataframe, we want to merge them.
                if not parkourDF.empty:
                    mergeDF = pd.merge(
                        ssdf[ssdf['Lane'] == lane],
                        parkourDF,
                        how='left',
                        on=[
                            'Sample_ID',
                            'Sample_Name',
                            'Sample_Project',
                            'Description'
                        ]
                    )
                    ssDic[key] = {'sampleSheet': mergeDF, 'lane': lane}
                else:
                    ssDic[key] = {
                        'sampleSheet': ssdf[ssdf['Lane'] == lane],
                        'lane': lane
                    }
            del self.fullSS
            return ssDic
        else:
            laneStr = self.flowcell + '_lanes_' + '_'.join(
                [str(lane) for lane in range(1, self.runInfoLanes + 1, 1)]
            )
            if not parkourDF.empty:
                mergeDF = pd.merge(
                        ssdf,
                        parkourDF,
                        how='left',
                        on=[
                            'Sample_ID',
                            'Sample_Name',
                            'Sample_Project',
                            'Description'
                        ]
                    )
                ssDic[laneStr] = {'sampleSheet': mergeDF, 'lane': 'all'}
            else:
                ssDic[laneStr] = {'sampleSheet': ssdf, 'lane': 'all'}
        del self.fullSS
        return ssDic

    def queryParkour(self, config):
        log.info("Pulling {} with pullURL".format(self.flowcell))
        return pullParkour(self.flowcell, config)

    def __init__(self, sampleSheet, lanes, config):
        log.warning("initiating sampleSheetClass")
        self.runInfoLanes = lanes
        self.origSs = sampleSheet
        self.flowcell = sampleSheet.split('/')[-2]
        self.ssDic = self.parseSS(self.queryParkour(config))
