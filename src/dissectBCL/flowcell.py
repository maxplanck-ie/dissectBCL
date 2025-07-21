from dissectBCL.demux import detMask, misMatcher
from dissectBCL.demux import writeDemuxSheet, writeDemuxSheetAviti, readDemuxSheet, compareDemuxSheet
from dissectBCL.demux import evalMiSeqP5, parseStats
from dissectBCL.demux import matchingSheets
from dissectBCL.postmux import renameProject
from dissectBCL.postmux import validateFqEnds
from dissectBCL.postmux import qcs, clumper, kraken, md5_multiqc, moveOptDup
from dissectBCL.fakeNews import pullParkour, mailHome
from dissectBCL.fakeNews import shipFiles, pushParkour
from dissectBCL.fakeNews import gatherFinalMetrics
from dissectBCL.misc import umlautDestroyer, P5Seriesret

import datetime
from dominate.tags import html, div, br
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from random import randint
import ruamel.yaml
import shutil
from subprocess import Popen, PIPE
import sys
from tabulate import tabulate
import xml.etree.ElementTree as ET
import json

class flowCellClass:

    # Init - fileChecks
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
            logging.info(f"Init - Checking {f}")
            if not f.exists():
                logging.critical(f"Init - {f} doesn't exist. Exiting")
                mailHome(
                    self.name,
                    f"{f} does not exist. dissectBCL crashed.",
                    self.config
                )
                sys.exit(1)

    # Init - Parse runInfo
    def parseRunInfo(self):
        """
        Takes the path to runInfo.xml and parses it.
        Returns:
         - a sequencing recipe dictionary. {'Read1':100, 'Index1':10, ... }
         - number of lanes (int)
         - the instrument (str)
         - the flowcellID (str)
        """
        logging.info("Init - Parsing RunInfo.xml")
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

    # Init - Parse runInfo
    def parseRunInfoAviti(self):
        """
        Takes the path to RunParameters.json and parses it.
        Returns:
         - a sequencing recipe dictionary. {'Read1':100, 'Index1':10, ... }
         - number of lanes (int)
         - the instrument (str)
         - the flowcellID (str)
        """
        logging.info("Init - Parsing RunParemeters.json")
        with open (self.runInfo) as run_json:
            run_params = json.load(run_json)
        cycles_dict = run_params['Cycles']
        seqRecipe = {}
        readCount = 1
        indexCount = 1
        seqRecipe['Read1'] = ['Y',cycles_dict['R1']]
        if 'R2' in cycles_dict.keys():
            seqRecipe['Read2'] = ['Y',cycles_dict['R2']]
        if 'I1' in cycles_dict.keys():
            seqRecipe['Index1'] = ['I', cycles_dict['I1']]
        if 'I2' in cycles_dict.keys():
            seqRecipe['Index2'] = ['I', cycles_dict['I2']]
        if run_params['AnalysisLanes'] == "1+2":
            lanes = 2
        else:
            lanes = 1

        instrument = run_params["InstrumentName"]
        flowcellID = run_params["FlowcellID"]
        return seqRecipe, lanes, instrument, flowcellID

    # Init - Validate successful run.
    def validateRunCompletion(self):
        """
        validates succesfull completion status in xml.
        """
        logging.info("Init - validateRunCompletion")
        if self.sequencer == 'Miseq':
            tree = ET.parse(self.runCompletionStatus)
            root = tree.getroot()
            for i in root.iter():
                if i.tag == 'CompletionStatus':
                    _status = i.text
        else:
            # no RunCompletionStatus.xml in novaseq, assume succes.
            _status = 'SuccessfullyCompleted'
        return (_status)

    # demux - prepConvert
    def prepConvert(self):
        '''
        Determines mask, dualIx status, PE status, convertOptions and mismatches
        '''
        logging.info("Demux - prepConvert - determine masking, indices, paired ends, and other options")
        for outputFolder in self.sampleSheet.ssDic:
            ss_dict = self.sampleSheet.ssDic[outputFolder]
            ss = ss_dict['sampleSheet']

            # determine mask, dualIx, PE, convertOpts, minP5, minP7 from seqRecipe
            (ss_dict['mask'], ss_dict['dualIx'], ss_dict['PE'],
            ss_dict['convertOpts'], minP5, minP7) = detMask(
                self.seqRecipe,
                ss,
                outputFolder
            )

            # extra check to make sure all our indices are of equal size!
            for min_ix, ix_str in ((minP5, 'index'), (minP7, 'index2')):
                if min_ix and not np.isnan(min_ix):
                    ss[ix_str] = ss[ix_str].str[:min_ix]

            # determine mismatch
            ss_dict['mismatch'] = misMatcher(ss['index'], P5Seriesret(ss))
        logging.info("Demux - prepConvert - mask in sampleSheet updated.")
        self.exitStats['premux'] = 0

    # demux - demux
    def demux(self):
        # Double check for run failure
        if self.succesfullrun != 'SuccessfullyCompleted':
            logging.warning("Demux - Illumina - Run not succesfull, marking as failed.")
            for outLane in self.sampleSheet.ssDic:
                Path(self.outBaseDir, outLane).mkdir(exists_ok=True)
                Path(self.outBaseDir, outLane, 'run.failed').touch()
            mailHome(
                f"{self.name} ignored",
                "RunCompletionStatus is not successfullycompleted.\n" +
                "Marked for failure and ignored for the future.",
                self.config,
                toCore=True
            )
        else:
            logging.info("Demux - Illumina - Run succesfull, starting demux")
            for outLane in self.sampleSheet.ssDic:
                logging.info(f"Demux - {outLane}")
                _ssDic = self.sampleSheet.ssDic[outLane]
                # Set outputDirectory
                outputFolder = Path(self.outBaseDir, outLane)
                outputFolder.mkdir(exist_ok=True)
                demuxOut = outputFolder / 'demuxSheet.csv'
                # Don't remake if demuxSheet exist
                if not demuxOut.exists():
                    logging.info(f"Demux - Writing demuxSheet for {outLane}")
                    writeDemuxSheet(
                        demuxOut,
                        _ssDic,
                        self.sampleSheet.laneSplitStatus
                    )
                else:
                    logging.warning(
                        f"Demux - demuxSheet for {outLane} already exists, not changing it."
                    )
                    compareDemuxSheet(_ssDic, demuxOut)

                # Don't run bcl-convert if we have the touched flag.
                if not Path(outputFolder, 'bclconvert.done').exists():
                    # Purge pre-existing reports / log folders
                    if Path(outputFolder, 'Reports').exists():
                        shutil.rmtree(Path(outputFolder, 'Reports'))
                    if Path(outputFolder, 'Logs').exists():
                        shutil.rmtree(Path(outputFolder, 'Logs'))
                    # Run bcl-convert
                    bclOpts = [
                        self.config['software']['bclconvert'],
                        '--output-directory', outputFolder,
                        '--force',
                        '--bcl-input-directory', self.bclPath,
                        '--sample-sheet', demuxOut,
                        '--bcl-num-conversion-threads', f"{int(self.config['misc']['threads'])//2}",
                        '--bcl-num-compression-threads', f"{int(self.config['misc']['threads'])//2}",
                        "--bcl-sampleproject-subdirectories", "true",
                    ]
                    if not self.sampleSheet.laneSplitStatus:
                        bclOpts.append('--no-lane-splitting')
                        bclOpts.append('true')
                    logging.info("Demux - Starting BCLConvert")
                    logging.info(f"Demux - {bclOpts}")
                    bclRunner = Popen(bclOpts,stdout=PIPE, stderr=PIPE)
                    _stdout, _stderr = bclRunner.communicate()
                    exitcode = bclRunner.returncode
                    if exitcode == 0:
                        logging.info("Demux - bclConvert exit 0")
                        Path(outputFolder, 'bclconvert.done').touch()
                        if self.sequencer == 'MiSeq':
                            if evalMiSeqP5(
                                outputFolder,
                                _ssDic['dualIx'],
                            ):
                                logging.info("Demux - P5 RC triggered.")
                                # Purge existing reports.
                                logging.info("Demux - Purge existing Reports folder")
                                shutil.rmtree(Path(outputFolder, 'Reports'))
                                shutil.rmtree(Path(outputFolder, 'Logs'))
                                # Rerun BCLConvert
                                logging.info("Demux - Rerun BCLConvert")
                                bclRunner = Popen(bclOpts,stdout=PIPE, stderr=PIPE)
                                _stdout, _stderr = bclRunner.communicate()
                                logging.info(f"Demux - bclConvert P5fix exit {exitcode}")
                                # Update the sampleSheet with proper RC'ed indices.
                                _ssDic['sampleSheet'] = matchingSheets(
                                    _ssDic['sampleSheet'], 
                                    readDemuxSheet(demuxOut, what='df')
                                )
                                _ssDic['P5RC'] = True
                            else:
                                _ssDic['P5RC'] = False
                        else:
                            _ssDic['P5RC'] = False
                    else:
                        logging.critical(f"Demux - BCLConvert exit {exitcode}")
                        mailHome(
                            outLane,
                            f"BCL-convert exit {exitcode}. Pipeline crashed. {_stderr.decode('utf-8')}",
                            self.config,
                            toCore=True
                        )
                        sys.exit(1)

                logging.info(f"Demux - Parsing stats for {outLane}")
                _ssDic['sampleSheet'] = parseStats(outputFolder, _ssDic['sampleSheet'])
            self.exitStats['demux'] = 0

    def demux_aviti(self):
        logging.info("Demux - Aviti system.")
        if self.succesfullrun != 'SuccessfullyCompleted':
            logging.warning("Demux - Aviti - Run not succesfull, marking as failed.")
            for outLane in self.sampleSheet.ssDic:
                Path(self.outBaseDir, outLane).mkdir(exists_ok=True)
                Path(self.outBaseDir, outLane, 'run.failed').touch()
            mailHome(
                f"{self.name} ignored",
                "RunCompletionStatus is not successfullycompleted.\n" +
                "Marked for failure and ignored for the future.",
                self.config,
                toCore=True
            )
        else:
            logging.info("Demux - Aviti - Run succesfull, starting demux")
            for outLane in self.sampleSheet.ssDic:
                logging.info(f"Demux - {outLane}")
                _ssDic = self.sampleSheet.ssDic[outLane]
                # P5RC is some legacy leftover from MiSeq illumina (needed to RC).
                # To keep data structures consistent, we set it to False.
                _ssDic['P5RC'] = False
                # Set outputDirectory
                outputFolder = Path(self.outBaseDir, outLane)
                outputFolder.mkdir(exist_ok=True)
                (outputFolder / 'manifest').mkdir(exist_ok=True)
                # Ship over RunManifest.csv, only if it doesn't exist yet.
                demuxOut=Path(outputFolder, 'manifest', 'RunManifest.csv')
                if not demuxOut.exists():
                    logging.info(f"Demux - Copying RunManifest.csv to {outputFolder}")
                    #shutil.copy(self.origSS, outputFolder / 'manifest' / 'RunManifest.csv')
                    writeDemuxSheetAviti(
                        demuxOut,
                        _ssDic,
                        self.sampleSheet.laneSplitStatus
                    )
                else:
                    logging.warning(
                        f"Demux - RunManifest.csv for {outLane} already exists, not changing it."
                    )
                # Run bases2fastq
                b2fOpts = [
                    self.config['software']['bases2fastq'],
                    '--run-manifest', Path(outputFolder, 'manifest', 'RunManifest.csv'),
                    '--num-threads', f"{self.config['misc']['threads']}",
                    '--group-fastq',
                    self.bclPath,
                    Path(outputFolder),
                ]

                if not Path(outputFolder, 'bases2fastq.done').exists():
                    logging.info("Demux - Starting bases2fastq")
                    logging.info(f"Demux - {b2fOpts}")
                    b2fRunner = Popen(b2fOpts, stdout=PIPE, stderr=PIPE)
                    _stdout, _stderr = b2fRunner.communicate()
                    exitcode = b2fRunner.returncode

                    if exitcode == 0:
                        logging.info("Demux - bases2fastq exit 0")
                        Path(outputFolder, 'bases2fastq.done').touch()
                    else:
                        logging.critical(f"Demux - bases2fastq exit {exitcode}")
                        mailHome(
                            outLane,
                            f"Bases2fastq exit {exitcode}. Pipeline crashed. {_stderr.decode('utf-8')}",
                            self.config,
                            toCore=True
                        )
                        sys.exit(1)
                else:
                    logging.info("Bases2fastq.done already exists. Moving forward.")
                
                logging.info(f"Demux - Parsing stats for {outLane}")
                _ssDic['sampleSheet'] = parseStats(outputFolder, _ssDic['sampleSheet'], mode='aviti')

            self.exitStats['demux'] = 0
        
    # postmux - postmux
    def postmux(self):
        logging.info("Postmux - Demux complete, starting postmux")
        for outLane in self.sampleSheet.ssDic:
            _ssDic = self.sampleSheet.ssDic[outLane]
            laneFolder = Path(self.outBaseDir, outLane)
            renameFlag = laneFolder / 'renamed.done'
            postmuxFlag = laneFolder / 'postmux.done'
            df = _ssDic['sampleSheet']
            # A dropna is needed here, as per aviti runs phiX could be listed as a 'project'.
            projects = list(df['Sample_Project'].dropna().unique())
            
            for project in projects:
                if not renameFlag.exists():
                    logging.info("Postmux - renaming {}".format(outLane))
                    if self.sequencer == 'aviti':
                        # Aviti mode, rename project folder
                        renameProject(laneFolder / 'Samples' / project, df, self.sampleSheet.laneSplitStatus)
                    else:
                        renameProject(laneFolder / project, df, self.sampleSheet.laneSplitStatus)
                    validateFqEnds(laneFolder / project, self)
            renameFlag.touch()
            for project in projects:
                if not postmuxFlag.exists():
                    _sIDs = set(df[df['Sample_Project'] == project]['Sample_ID'])
                    # FQC
                    logging.info(f"Postmux - FastQC {outLane} - {project}")
                    qcs(project, laneFolder, _sIDs, self.config)
                    # Clump
                    logging.info(f"Postmux - Clumping {outLane} - {project}")
                    clumper(project, laneFolder, _sIDs, self.config, _ssDic['PE'], self.sequencer)
                    # kraken
                    logging.info(f"Postmux - kraken {outLane} - {project}")
                    kraken(project, laneFolder, _sIDs, self.config)
                    # multiQC
                    logging.info(f"Postmux - md5/multiqc {outLane} - {project}")
                    md5_multiqc(project, laneFolder, self)
                    # Move optical duplicates
                    moveOptDup(laneFolder)
            postmuxFlag.touch()
        self.exitStats['postmux'] = 0

    # fakenews
    def fakenews(self):
        logging.info("fakenews - Postmux complete, starting fakenews.")
        for outLane in self.sampleSheet.ssDic:
            # Ship Files
            logging.info(f"fakenews - shipFiles - {outLane}")
            self.exitStats[outLane] = shipFiles(self.outBaseDir / outLane, self.config)

            # Push parkour
            logging.info(f"fakenews - pushParkour - {outLane}")
            self.exitStats[outLane]['pushParkour'] = pushParkour(self.flowcellID, self.sampleSheet, self.config, self.bclPath, self.sequencer)

            # diagnoses / QCstats
            logging.info(f"fakenews - gatherMetrics - {outLane}")
            _h = drHouseClass(gatherFinalMetrics(outLane, self))
            logging.info(f"fakenews - prepMail - {outLane}")
            subject, _html = _h.prepMail()
            mailHome(subject, _html, self.config)
            (self.outBaseDir / outLane / 'communication.done').touch()

    # organiseLogs
    def organiseLogs(self):
        for outLane in self.sampleSheet.ssDic:
            logging.info(f"organiseLogs - Populating log dir for {outLane}")
            _logDir = self.outBaseDir / outLane / 'Logs'
            _logDir.mkdir(exist_ok=True)

            # Write out ssdf.
            outssdf = _logDir / 'sampleSheetdf.tsv'
            self.sampleSheet.ssDic[outLane]['sampleSheet'].to_csv(outssdf, sep='\t')

            # write out outLaneInfo.yaml
            dic0 = self.sampleSheet.ssDic[outLane]
            del dic0['sampleSheet']
            yaml0 = ruamel.yaml.YAML()
            yaml0.indent(mapping=2, sequence=4, offset=2)
            with open(_logDir / 'outLaneInfo.yaml', 'w') as f:
                yaml0.dump(dic0, f)

            # write out config.ini
            dic1 = self.asdict()
            with open(_logDir / 'config.ini', 'w') as f:
                dic1['config'].write(f)

            # write out flowcellInfo.yaml
            del dic1['config']
            yaml1 = ruamel.yaml.YAML()
            yaml1.indent(mapping=2, sequence=4, offset=2)
            with open(_logDir / 'flowcellInfo.yaml', 'w') as f:
                yaml1.dump(dic1, f)

    def __init__(self, name, bclPath, logFile, config, sequencer):
        sequencers = {
            'A': 'NovaSeq',
            'N': 'NextSeq',
            'M': 'MiSeq'
        }
        logging.warning("Initiating flowcellClass {}".format(name))
        self.name = name
        self.bclPath = Path(bclPath)
        self.outBaseDir = Path(config['Dirs']['outputDir'])
        self.logFile = logFile
        self.config = config

        if sequencer == 'illumina':
            # Illumina mode.
            self.inBaseDir = Path(config['Dirs']['baseDir_illumina'])
            self.sequencer = sequencers[name.split('_')[1][0]]
            self.origSS = Path(bclPath, 'SampleSheet.csv')
            self.runInfo = Path(bclPath, 'RunInfo.xml')
            self.runCompletionStatus = Path(bclPath, 'RunCompletionStatus.xml')
            # Run filesChecks
            self.filesExist()
            self.succesfullrun = self.validateRunCompletion()
            # populate runInfo vars.
            self.seqRecipe, \
                self.lanes, \
                self.instrument, \
                self.flowcellID = self.parseRunInfo()
        else:
            # Aviti mode.
            self.inBaseDir = Path(config['Dirs']['baseDir_aviti'])
            self.sequencer = sequencer
            self.origSS = Path(bclPath, 'RunManifest.csv')
            self.runInfo = Path(bclPath, 'RunParameters.json')
            self.succesfullrun = 'SuccessfullyCompleted'
            self.seqRecipe, \
            self.lanes, \
            self.instrument, \
            _fid = self.parseRunInfoAviti() # discard flowcell ID until it can be used consistently with Parkour
            self.flowcellID = self.name        
        
        self.startTime = datetime.datetime.now()
        # Create sampleSheet information
        self.sampleSheet = sampleSheetClass(
            self.origSS,
            self.lanes,
            self.sequencer,
            self.config
        )
        self.exitStats = {}
        self.transferTime = None

    def asdict(self):
        return {
            'name': self.name,
            'sequencer': self.sequencer,
            'bclPath': str(self.bclPath),
            'original sampleSheet': str(self.origSS),
            'runInfo': str(self.runInfo),
            'runCompletionStatus': str(self.runCompletionStatus) if hasattr(self, 'runCompletionStatus') and self.runCompletionStatus else "",
            'sucessfulRun': self.succesfullrun,
            'inBaseDir': str(self.inBaseDir),
            'outBaseDir': str(self.outBaseDir),
            'dissect logFile': str(self.logFile),
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

    def decideSplit(self,aviti):
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
        
        # it would actually be nicer to handle the column names on the class level
        #lane colname are identical between illumina and aviti
        lane_colname = 'Lane'
        #the other colnames are deviating -> set illumina colnames as default
        sample_colname = 'Sample_ID'
        index1_colname = 'index'
        index2_colname = 'index2'
        project_colname = 'Sample_Project'

        #deviating colnames in aviti
        if aviti:
            sample_colname = 'SampleName'
            index1_colname = "Index1"
            index2_colname = "Index2"
            project_colname = "Project"


            logging.info("Deciding lanesplit using Aviti colnames.")
        else:
            logging.info("Deciding lanesplit using Illumina colnames.")

        laneSplitStatus = True
        # Do we need lane splitting or not ?
        # If there is at least one sample in more then 1 lane, we cannot split:
        samples = list(self.fullSS[sample_colname].unique())
        if 'PhiX' in samples:
            samples.remove("PhiX")
        for _s in samples:
            if len(
                list(self.fullSS[
                    self.fullSS[sample_colname] == _s
                ][lane_colname].unique()
                )
            ) > 1:
                logging.info(
                    "No lane splitting: >= 1 sample in multiple lanes."
                )
                laneSplitStatus = False
        # If one project is split over multiple lanes, we also don't split:
        projects = list(self.fullSS[project_colname].unique())
        if '0000_PhiX_DeepSeq' in projects:
            projects.remove("0000_PhiX_DeepSeq")
        for project in projects:
            if len(
                list(self.fullSS[
                    self.fullSS[project_colname] == project
                ][lane_colname].unique()
                )
            ) > 1:
                logging.info(
                    "No lane splitting: >= 1 project in multiple lanes."
                )
                laneSplitStatus = False
        # Don't split if 1 lane in ss, multiple in runInfo
        if len(list(self.fullSS[lane_colname].unique())) < self.runInfoLanes:
            logging.info(
                "No lane splitting: 1 lane listed, {} found.".format(
                    self.runInfoLanes
                )
            )
            laneSplitStatus = False
        # Make sure:
        # if laneSplitStatus = False:
        # No samples can clash at all!
        if lane_colname in list(
            self.fullSS.columns
        ) and not laneSplitStatus:
            if index1_colname in list(
                self.fullSS.columns
            ) and index2_colname in list(
                self.fullSS.columns
            ):
                tmpSheet = self.fullSS[[sample_colname, index1_colname, index2_colname]]
                # A sample can sit in multiple lanes
                # Deduplicate id - ix, ix2
                tmpSheet = tmpSheet.drop_duplicates()
                # combine index1 and index2
                testSer = tmpSheet[index1_colname] + tmpSheet[index2_colname]
            elif index1_colname in list(self.fullSS.columns):
                tmpSheet = self.fullSS[[sample_colname, index1_colname]]
                # same logic as above.
                tmpSheet = tmpSheet.drop_duplicates()
                testSer = tmpSheet[index1_colname]
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
    def parseSS(self, parkourDF) -> dict:
        """
        We read the sampleSheet csv, and remove the stuff above the header.
        """
        logging.info("parseSS - parsing sampleSheet and parkour, Illumina mode.")
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
        # Sanitize projects names
        ssdf['Sample_Project'] = ssdf['Sample_Project'].apply(
            lambda x: umlautDestroyer(x)
        )
        # Sanitize sample names
        ssdf['Sample_Name'] = ssdf['Sample_Name'].apply(
            lambda x: umlautDestroyer(x)
        )

        self.fullSS = ssdf
        self.laneSplitStatus =self.decideSplit(aviti=False)
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
                    if '-' in self.flowcell:
                        '''
                        In case of miSeq runs,
                        assume the requested depth is 20/#samples
                        this is due to the 10M / sample
                        minimum for parkour requests
                        we do this here (and not in pullparkour)
                        since parkour returns all
                        samples, not necesarily those sequenced.
                        '''
                        newReqDepth = 20/len(
                            list(mergeDF['Sample_Name'].unique())
                        ) * 1000000
                        newReqDepth = round(newReqDepth, 0)
                        mergeDF['reqDepth'] = newReqDepth
                        logging.debug(
                            'miSEQ detected, override seqdepth: {}'.format(
                                newReqDepth
                            )
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

    def parseSS_aviti(self, parkourDF) -> dict:
        logging.info("parseSS - parsing sampleSheet and parkour, Aviti mode.")
        logging.info("Reading sampleSheet (RunManifest.csv).")
        
        # Resort to parsing RunManifest.csv line by line to get settings.
        sampleline = None
        maskstrdict = {
            'R1FastQMask': '' ,
            'R2FastQMask': '',
            'I1Mask': '',
            'I2Mask': ''
            }
        dualIx = True
        PE = True
        with open(self.origSs, 'r') as f:
            for ix, line in enumerate(f):
                if 'R1FastQMask' in line:
                    maskstrdict['R1FastQMask']=line.split(',')[1].strip()
                if 'R2FastQMask' in line:
                    maskstrdict['R2FastQMask']=line.split(',')[1].strip()
                    PE = True
                if 'I1Mask' in line:
                    maskstrdict['I1Mask']=line.split(',')[1].strip()
                if 'I2Mask' in line:
                    maskstrdict['I2Mask']=line.split(',')[1].strip()
                    dualIx = True
                if '[samples]' in line.strip().lower():
                    sampleline = ix + 1
        ssdf = pd.read_csv(self.origSs, skiprows=sampleline, sep=',')
        #don't rename
        #ssdf.rename(columns={'SampleName': 'Sample_ID'}, inplace=True)
        if 'Project' not in ssdf.columns:
            logging.critical("RunManifest.csv does not contain a 'Project' column.")
            mailHome(
                self.flowcell,
                "RunManifest.csv does not contain a 'Project' column.",
                self.config,
                toCore=False
            )
            sys.exit(1)

        #check if Description column is present, if not, add an empty one
        if 'Description' not in ssdf.columns:
            logging.info('Description column not found in RunManifest. Adding an empty column.')
            ssdf = ssdf.assign(Description='')

        self.fullSS = ssdf
        self.laneSplitStatus = self.decideSplit(aviti=True)

        if dualIx:
            mmd = {'I1MismatchThreshold': 0, 'I2MismatchThreshold': 0}
        else:
            mmd = {'I1MismatchThreshold': 0}
        
        ssDic = {}
        if self.laneSplitStatus:
            for lane in range(1, self.runInfoLanes + 1, 1):
                key = self.flowcell + '_lanes_' + str(lane)
                # if we have a parkour dataframe, we want to merge them.
                if not parkourDF.empty:
                    mergeDF = pd.merge(
                        ssdf[ssdf['Lane'] == lane],
                        parkourDF.drop(columns='Description'),
                        how='left',
                        left_on=[
                            'SampleName',
                            'Project'
                            ],
                        right_on=[
                            'Sample_ID',
                            'Sample_Project',
                        ]
                    )
                ssDic[key] = {
                    'sampleSheet': ssdf[ssdf['Lane'] == lane],
                    'lane': lane,
                    'mask': maskstrdict,
                    'dualIx': dualIx,
                    'PE': PE,
                    'convertOpts': [],
                    'mismatch': mmd
                    }
        else:
            laneLis = [
                str(lane) for lane in range(1, self.runInfoLanes + 1, 1)
            ]
            laneStr = self.flowcell + '_lanes_' + '_'.join(
                laneLis
            )
            dfLaneEntry = '+'.join(laneLis)
            if not parkourDF.empty:
                mergeDF = pd.merge(
                        ssdf,
                        parkourDF.drop(columns='Description'),
                        how='left',
                        left_on=[
                            'SampleName',
                            'Project'
                            ],
                        right_on=[
                            'Sample_ID',
                            'Sample_Project',
                        ]
                    )
                # Collate if one samples is split on multiple lanes.
                mergeDF['Lane'] = mergeDF['Lane'].astype(str)
                aggDic = {}
                for col in list(mergeDF.columns):
                    if col == 'Lane':
                        aggDic[col] = '+'.join
                    elif col != 'Sample_ID':
                        aggDic[col] = 'first'
                mergeDF = mergeDF.groupby(
                    'Sample_ID'
                ).agg(aggDic).reset_index()
                mergeDF['Lane'] = dfLaneEntry
                ssDic[laneStr] = {'sampleSheet': mergeDF,
                                  'lane': 'all',
                                  'mask': maskstrdict,
                                  'dualIx': dualIx,
                                  'PE': PE,
                                  'convertOpts': [],
                                  'mismatch': mmd}
            else:
                ssdf['Lane'] = dfLaneEntry
                ssDic[laneStr] = {'sampleSheet': ssdf,
                                  'lane': 'all',
                                  'mask': maskstrdict,
                                  'dualIx': dualIx,
                                  'PE': PE,
                                  'convertOpts': [],
                                  'mismatch': mmd}

 #       ssDic = {self.flowcell: {
 #           'sampleSheet': mergeDF,
 #           'lane': 'all',
 #           'mask': maskstr,
 #           'dualIx': dualIx,
 #           'PE': PE,
 #           'convertOpts': [],
 #           'mismatch': mmd
 #           }
 #       }
        del self.fullSS
        return ssDic

    def queryParkour(self, config, aviti=False):
        logging.info("Pulling {} with pullURL".format(self.flowcell))
        return pullParkour(self.flowcell, config, aviti)

    def __init__(self, sampleSheet, lanes, sequencer, config):
        logging.warning("initiating sampleSheetClass")
        self.origSs = sampleSheet
        self.flowcell = sampleSheet.parts[-2]
        if sequencer == 'aviti':
            self.runInfoLanes = 2
            self.ssDic = self.parseSS_aviti(self.queryParkour(config, aviti=True))
        else:
            self.runInfoLanes = lanes
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
        message += f"Flowcell: {self.flowcellID}\n"
        message += f"Sequencer: {self.sequencer}\n"
        message += f"outLane: {self.outLane}\n"
        message += f"Runtime: {self.runTime}\n"
        message += f"transferTime: {self.transferTime}\n"
        message += f"Space Free (rapidus): {self.spaceFree_rap[1]} GB - {spaceGood(self.spaceFree_rap[1])}\n"
        message += f"Space Free (solexa): {self.spaceFree_sol[1]} GB - {spaceGood(self.spaceFree_sol[1])}\n"
        message += f"barcodeMask: {self.barcodeMask}\n"
        message += self.mismatch + '\n'
        # Undetermined
        if isinstance(self.undetermined, str):
            message += f"Undetermined indices: {self.undetermined}\n"
        elif isinstance(self.undetermined, int):
            message += f"Undetermined indices: {round(100*self.undetermined/self.totalReads, 2)}% ({round(self.undetermined/1000000, 0)}M)\n"
        # exitStats
        for key in self.exitStats:
            if key in [
                'premux', 'demux', 'postmux'
            ]:
                message += f"exit {key}: {self.exitStats[key]}\n"
            elif key == self.outLane:
                for subkey in self.exitStats[key]:
                    message += f"return {subkey}: {self.exitStats[key][subkey]}\n"

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
        if not self.P5RC:
            P5RCstr = ''
        else:
            P5RCstr = '\n\n <FONT COLOR=red> '
            P5RCstr += 'Note that the P5s have been reverse complemented '
            P5RCstr += 'automatically. </strong></FONT> \n\n'
            P5RCstr += 'The multiqc report contains '
            P5RCstr += 'the index sequences '
            P5RCstr += 'as they are used for demultiplexing.'

        msg = _html.render() +\
            P5RCstr +\
            '<h3>Top unknown barcodes</h3>' +\
            tabulate(undtableCont, undtableHead, tablefmt="html") +\
            '<h3>Samples</h3>' +\
            tabulate(
                tableCont, tableHead, tablefmt="html", disable_numparse=True
            )
        return (self.outLane, msg)

    def __init__(self, qcdic):
        self.undetermined = qcdic['undetermined']
        self.totalReads = qcdic['totalReads']
        self.topBarcodes = qcdic['topBarcodes']
        self.spaceFree_rap = qcdic['spaceFree_rap']
        self.spaceFree_sol = qcdic['spaceFree_sol']
        self.runTime = qcdic['runTime']
        self.optDup = qcdic['optDup']
        self.flowcellID = qcdic['flowcellID']
        self.outLane = qcdic['outLane']
        self.contamination = qcdic['contamination']
        self.barcodeMask = qcdic['barcodeMask']
        self.mismatch = qcdic['mismatch']
        self.transferTime = qcdic['transferTime']
        self.exitStats = qcdic['exitStats']
        self.P5RC = qcdic['P5RC']
        self.sequencer = qcdic['sequencer']
