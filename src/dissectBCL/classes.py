import os
import xml.etree.ElementTree as ET
import sys
from dissectBCL.fakeNews import pullParkour, mailHome
import pandas as pd
import datetime
from tabulate import tabulate
from random import randint
from dominate.tags import html, div, br
import logging


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
            logging.info("Checking {}".format(f))
            if not os.path.exists(f):
                logging.critical("{} doesn't exist. Exiting".format(f))
                mailHome(
                    self.name,
                    "{} does not exist. dissectBCL crashed.".format(
                        f
                    ),
                    self.config
                )
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
        logging.info("Parsing RunInfo.xml")
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
        logFile,
        config
    ):
        sequencers = {
            'A': 'NovaSeq',
            'N': 'NextSeq',
            'M': 'MiSeq'
        }
        logging.warning("Initiating flowcellClass {}".format(name))
        self.name = name
        self.sequencer = sequencers[name.split('_')[1][0]]
        self.bclPath = bclPath
        self.origSS = origSS
        self.runInfo = runInfo
        self.inBaseDir = inBaseDir
        self.outBaseDir = outBaseDir
        self.logFile = logFile
        self.config = config
        # Run filesChecks
        self.filesExist()
        # populate runInfo vars.
        self.seqRecipe, \
            self.lanes, \
            self.instrument, \
            self.flowcellID = self.parseRunInfo()
        self.startTime = datetime.datetime.now()

    def asdict(self):
        return {
            'name': self.name,
            'sequencer': self.sequencer,
            'bclPath': self.bclPath,
            'original sampleSheet': self.origSS,
            'runInfo': self.runInfo,
            'inBaseDir': self.inBaseDir,
            'outBaseDir': self.outBaseDir,
            'dissect logFile': self.logFile,
            'seqRecipe': self.seqRecipe,
            'lanes': self.lanes,
            'instrument': self.instrument,
            'flowcellID': self.flowcellID,
            'Time initiated': self.startTime.strftime("%m/%d/%Y, %H:%M:%S"),
            'config': self.config
        }


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
        logging.info("Deciding lanesplit.")
        laneSplitStatus = True
        # Do we need lane splitting or not ?
        # If there is at least one sample in more then 1 lane, we cannot split:
        if sum(self.fullSS['Sample_Name'].value_counts() > 1) > 0:
            logging.info(
                "No lane splitting: >= 1 sample in multiple lanes."
            )
            laneSplitStatus = False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.fullSS['Sample_Project'].unique())
        for project in projects:
            if len(
                list(self.fullSS[
                    self.fullSS['Sample_Project'] == project
                ]['Lane'].unique()
                )
            ) > 1:
                logging.info(
                    "No lane splitting: >= 1 project in multiple lanes."
                )
                laneSplitStatus = False
        # Don't split if 1 lane in ss, multiple in runInfo
        if len(list(self.fullSS['Lane'].unique())) < self.runInfoLanes:
            logging.info(
                "No lane splitting: 1 lane listed, {} found.".format(
                    self.runInfoLanes
                )
            )
            laneSplitStatus = False
        # Make sure:
        # if laneSplitStatus = False:
        # No samples can clash at all!
        if 'Lane' in list(
            self.fullSS.columns
        ) and not laneSplitStatus:
            if 'index' in list(
                self.fullSS.columns
            ) and 'index2' in list(
                self.fullSS.columns
            ):
                tmpSheet = self.fullSS[['Sample_ID', 'index', 'index2']]
                # A sample can sit in multiple lanes
                # Deduplicate id - ix, ix2
                tmpSheet = tmpSheet.drop_duplicates()
                # combine index1 and index2
                testSer = tmpSheet['index'] + tmpSheet['index2']
            elif 'index' in list(self.fullSS.columns):
                tmpSheet = self.fullSS[['Sample_ID', 'index']]
                # same logic as above.
                tmpSheet = tmpSheet.drop_duplicates()
                testSer = tmpSheet['index']
            for count in testSer.value_counts().to_dict().values():
                if count > 1:
                    logging.warning(
                        "Found a sample clash even though laneSplit == False."
                    )
                    logging.info("Overriding laneSplitStatus to True!")
                    laneSplitStatus = True
        logging.info('decide_lanesplit returns {}'.format(laneSplitStatus))
        return laneSplitStatus

    # Parse sampleSheet
    def parseSS(self, parkourDF):
        """
        We read the sampleSheet csv, and remove the stuff above the header.
        """
        logging.info("Reading sampleSheet.")
        ssdf = pd.read_csv(self.origSs, sep=',')
        ssdf.columns = ssdf.iloc[0]
        ssdf = ssdf.drop(ssdf.index[0])
        # Reset index
        ssdf.reset_index(inplace=True)
        ssdf.index.name = None
        # Remove 'level0' column
        ssdf.drop('level_0', axis=1, inplace=True)

        # ssdf = ssdf.dropna(axis=1, how='all')
        # NB: don't remove NAs, as it's possible that there are no e.g.
        # indices that are specified
        # ssdf = ssdf.dropna(axis=1, how='all')

        ssdf = ssdf.astype({'Lane': 'int32'})
        # Remove spaces if we have them
        ssdf['Sample_Project'] = ssdf['Sample_Project'].str.replace(' ', '')
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
                        parkourDF.drop(columns='Description'),
                        how='left',
                        on=[
                            'Sample_ID',
                            'Sample_Name',
                            'Sample_Project',
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
            laneLis = [
                str(lane) for lane in range(1, self.runInfoLanes + 1, 1)
            ]
            laneStr = self.flowcell + '_lanes_' + '_'.join(
                laneLis
            )
            dfLaneEntry = ','.join(laneLis)
            if not parkourDF.empty:
                mergeDF = pd.merge(
                        ssdf,
                        parkourDF.drop(columns='Description'),
                        how='left',
                        on=[
                            'Sample_ID',
                            'Sample_Name',
                            'Sample_Project',
                        ]
                    )
                # Collate if one samples is split on multiple lanes.
                mergeDF['Lane'] = mergeDF['Lane'].astype(str)
                aggDic = {}
                for col in list(mergeDF.columns):
                    if col == 'Lane':
                        aggDic[col] = ','.join
                    elif col != 'Sample_ID':
                        aggDic[col] = 'first'
                mergeDF = mergeDF.groupby(
                    'Sample_ID'
                ).agg(aggDic).reset_index()
                mergeDF['Lane'] = dfLaneEntry
                ssDic[laneStr] = {'sampleSheet': mergeDF, 'lane': 'all'}
            else:
                ssdf['Lane'] = dfLaneEntry
                ssDic[laneStr] = {'sampleSheet': ssdf, 'lane': 'all'}
        del self.fullSS
        return ssDic

    def queryParkour(self, config):
        logging.info("Pulling {} with pullURL".format(self.flowcell))
        return pullParkour(self.flowcell, config)

    def __init__(self, sampleSheet, lanes, config):
        logging.warning("initiating sampleSheetClass")
        self.runInfoLanes = lanes
        self.origSs = sampleSheet
        self.flowcell = sampleSheet.split('/')[-2]
        self.ssDic = self.parseSS(self.queryParkour(config))


class drHouseClass:
    def greeter(self):
        now = datetime.datetime.now()
        morning = [
            'Guten Morgen!',
            'Good Morning!',
            'Bonjour!',
            'Buon Giorno!',
            'Buenos Dias!',
            'Sobh Bekheir!',
            'Bună Dimineața!',
            'Goeiemorgen!',
            'Dzień Dobry'

        ]
        afternoon = [
            'Guten Tag!',
            'Good Afternoon!',
            'Bonne Après-midi!',
            'Buon Pomeriggio!',
            'Buenas Tardes!',
            'Bad Az Zohr Bekheir',
            'Bună Ziua!',
            'Goede Namiddag!',
            'Dzień Dobry'
        ]
        evening = [
            'Guten Abend!',
            'Good Evening!',
            'Bonsoir!',
            'Buona Serata!',
            'Buenas Noches!',
            'Asr Bekheir!',
            'Bună Seara!',
            'Goede Avond!',
            'Dobry wieczór'
        ]
        if now.hour < 12:
            return (morning[randint(0, 6)] + '\n\n')
        elif now.hour < 18:
            return (afternoon[randint(0, 6)] + '\n\n')
        else:
            return (evening[randint(0, 6)] + '\n\n')

    def prepMail(self):
        @staticmethod
        def spaceGood(freeSpace):
            if freeSpace > 5000:
                return "All good!"
            elif freeSpace > 1000 and freeSpace < 5000:
                return "Space getting tight!"
            else:
                return "Danger zone, clear space immediately!"

        _html = html()
        # Build message
        message = self.greeter()
        message += "Flowcell: {}\n".format(self.flowcellID)
        message += "outLane: {}\n".format(self.outLane)
        message += "Runtime: {}\n".format(self.runTime)
        message += "transferTime: {}\n".format(self.transferTime)
        message += "Space Free (rapidus): {} GB - {}\n".format(
            self.spaceFree_rap[1], spaceGood(self.spaceFree_rap[1])
        )
        message += "Space Free (solexa): {} GB - {}\n".format(
            self.spaceFree_sol[1], spaceGood(self.spaceFree_sol[1])
        )
        message += "barcodeMask: {}\n".format(self.barcodeMask)
        message += self.mismatch + '\n'
        # Undetermined
        if isinstance(self.undetermined, str):
            message += "Undetermined indices: {}\n".format(self.undetermined)
        elif isinstance(self.undetermined, int):
            message += "Undetermined indices: {}% ({}M)\n".format(
                round(100*self.undetermined/self.totalReads, 2),
                round(self.undetermined/1000000, 0)
            )
        # exitStats
        for key in self.exitStats:
            if key in [
                'premux', 'demux', 'postmux'
            ]:
                message += "exit {}: {}\n".format(key, self.exitStats[key])
            elif key == self.outLane:
                for subkey in self.exitStats[key]:
                    message += "return {}: {}\n".format(
                        subkey, self.exitStats[key][subkey]
                    )
        # undetermined table
        undtableHead = ["P7", "P5", "# reads (M)", "% of und. Reads"]
        undtableCont = []
        for comb in self.topBarcodes:
            combSplit = comb.split('+')
            if len(combSplit) == 1:
                undtableCont.append(
                    [
                        combSplit[0],
                        'NA',
                        self.topBarcodes[comb][0],
                        self.topBarcodes[comb][1] * 100
                    ]
                )
            else:
                undtableCont.append(
                    [
                        combSplit[0],
                        combSplit[1],
                        self.topBarcodes[comb][0],
                        self.topBarcodes[comb][1] * 100
                    ]
                )

        # append message
        _html.add(div((i, br()) for i in message.splitlines()))
        # build the table
        tableHead = [
            "Project",
            "Sample",
            "SampleID",
            "OptDup",
            "GottenReads",
            "GotvReq",
            "%fragments",
            "kraken",
            "parkour"
        ]
        tableCont = []

        def optDupRet(optDup):
            try:
                return (
                    round(optDup/100, 2)
                )
            except TypeError:
                return (optDup)

        for optLis in self.optDup:
            # if samples come from parkour but omitted in demuxsheet
            if optLis[1] not in self.contamination:
                krakfrag = 0
                krakenOrg = 'omitted'
                parkourOrg = 'omitted'
            else:
                try:
                    krakfrag = round(
                        self.contamination[optLis[1]][0] * 100, 1
                    )
                except TypeError:  # NA / no reads retrieved
                    krakfrag = 'NA'
                krakenOrg = self.contamination[optLis[1]][1].lower()
                parkourOrg = self.contamination[optLis[1]][2].lower()
            tableCont.append(
                [
                    optLis[0],  # Project
                    optLis[2],  # Sample
                    optLis[1],  # SampleID
                    optDupRet(optLis[3]),  # OptDup,
                    "{0:.1E}".format(optLis[5]),  # gotten reads
                    optLis[4],  # got/req
                    krakfrag,  # %frags kraken
                    krakenOrg,  # krakenOrg
                    parkourOrg  # parkourOrg
                ]
            )
        msg = _html.render() +\
            '<h3>Top unknown barcodes</h3>' +\
            tabulate(undtableCont, undtableHead, tablefmt="html") +\
            '<h3>Samples</h3>' +\
            tabulate(
                tableCont, tableHead, tablefmt="html", disable_numparse=True
            )
        return (self.outLane, msg)

    def __init__(
        self,
        undetermined,
        totalReads,
        topBarcodes,
        spaceFree_rap,
        spaceFree_sol,
        runTime,
        optDup,
        flowcellID,
        outLane,
        contamination,
        barcodeMask,
        mismatch,
        transferTime,
        exitStats
    ):
        self.undetermined = undetermined
        self.totalReads = totalReads
        self.topBarcodes = topBarcodes
        self.spaceFree_rap = spaceFree_rap
        self.spaceFree_sol = spaceFree_sol
        self.runTime = runTime
        self.optDup = optDup
        self.flowcellID = flowcellID
        self.outLane = outLane
        self.contamination = contamination
        self.barcodeMask = barcodeMask
        self.mismatch = mismatch
        self.transferTime = transferTime
        self.exitStats = exitStats
