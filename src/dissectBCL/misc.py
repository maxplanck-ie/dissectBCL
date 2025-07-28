import os
import configparser
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess as sp
from importlib.metadata import version
import sys
import logging
import shutil
import re
from rich import print
from typing import Optional, Literal

def getConf(configfile, quickload=False):
    config = configparser.ConfigParser()
    logging.info("Reading configfile from {}".format(configfile))
    config.read(configfile)
    if not quickload:
        # bcl-convertVer -> Illumina demultiplexer
        p = sp.run(
            [
                config['software']['bclconvert'],
                '--version'
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        bclconvert = p.stderr.decode().splitlines()[0].split(' ')[2]
        # bases2fastq -> Aviti demultiplexer
        p = sp.run(
            [
                config['software']['bases2fastq'],
                '--version'
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        bases2fastq = p.stdout.decode().splitlines()[0].split(' ')[2].split(',')[0]
        # fastqcVer
        p = sp.run(
            [
                'fastqc',
                '--version'
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        fastqc = p.stdout.decode().splitlines()[0].split(' ')[1]
        # kraken2 version
        p = sp.run(
            [
                'kraken2',
                '--version'
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        kraken2 = p.stdout.decode().splitlines()[0].split(' ')[2]
        # clumpifyVer
        p = sp.run(
            [
                'clumpify.sh',
                '--version'
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        clumpify = p.stderr.decode().splitlines()[1].split(' ')[2]
        # Set all the versions.
        config['softwareVers'] = {}
        config['softwareVers']['bclconvert'] = bclconvert
        config['softwareVers']['bases2fastq'] = bases2fastq
        config['softwareVers']['multiqc'] = version('multiqc')
        config['softwareVers']['kraken2'] = kraken2
        config['softwareVers']['bbmap'] = clumpify
        config['softwareVers']['fastqc'] = fastqc
        for soft, ver in config['softwareVers'].items():
            print("{} = {}".format(
                soft, ver
            ))
        # Double check if fastqc_adapters is set.
        if not Path(config['software']['fastqc_adapters']).exists():
            sys.exit(f"Fastqc adapters file under {config['software']['fastqc_adapters']} not found. Please check your config file.")
    return config


def getNewFlowCell(
    config,
    fPath: Optional[str] =None,
    sequencer:Optional[Literal['aviti', 'illumina']] = None
) -> tuple[Optional[str], Optional[Path], Optional[str]]:
    # If there is a fPath set, just return that.
    outBaseDir = Path(config['Dirs']['outputDir'])
    if fPath:
        assert sequencer in ('aviti', 'illumina'), "Sequencer must be set explicitely as 'aviti' or 'illumina' when providing a direct flow cell path."
        fPath = Path(fPath)
        assert fPath.exists()
        flowcellName = fPath.name
        flowcellDir = fPath
        if not any(outBaseDir.glob(f"{flowcellName}*/communication.done")):
            return (flowcellName, flowcellDir, sequencer)
        else:
            print(f"[red]{flowcellName} exists with a communication.done flag already.[/red]")
            sys.exit()
    
    # set some config vars.
    baseDir_illumina = config['Dirs']['baseDir_illumina']
    baseDir_aviti = config['Dirs']['baseDir_aviti']

    # Illumina
    # Glob over the baseDirs to get all flowcells.
    flowCells = list(Path(baseDir_illumina).glob('*/RTAComplete.txt'))
    _patterns = ['communication.done', 'fastq.made', 'run.failed']
    # Check if the flowcell exists in the output directory.
    for flowcell in flowCells:
        flowcellName = flowcell.parent.name
        flowcellDir = flowcell.parent
        # Make sure copycomplete exists.
        if (Path(flowcellDir) / 'CopyComplete.txt').exists():
            # Look for a folder containing the flowcellname.
            # no folder with name -> start the pipeline.
            if not any(outBaseDir.glob(f"{flowcellName}*")):
                return (flowcellName, flowcellDir, 'illumina')
            # If a matching folder exists, but no flag, start the pipeline:
            elif not any(any(outBaseDir.glob(f"{flowcellName}*/{pattern}")) for pattern in _patterns):
                return (flowcellName, flowcellDir, 'illumina')
    # Aviti
    flowCells = list(Path(baseDir_aviti).glob('*/RunUploaded.json'))
    _flowcellpattern = re.compile(r"^\d{8}_[\w-]+_[\w-]+$")
    for flowcell in flowCells:
        flowcellName = flowcell.parent.name

        assert len(flowcellName.split('_')) == 3, \
            f"Aviti flow cells need to be named as 'YYYYMMDD_sequencer_runID'. Instead received: {flowcellName}"
        assert _flowcellpattern.match(flowcellName), \
            f"Aviti flow cells need to match the pattern 'YYYYMMDD_sequencer_runID'. Instead received: {flowcellName}"
        flowcellDir = flowcell.parent
        print(f"flowcellName = {flowcellName}, flowcellDir = {flowcellDir}")
        print(any(outBaseDir.glob(f"{flowcellName}*/{pattern}") for pattern in _patterns))
        # Look for a folder containing the flowcellname.
        # no folder with name -> start the pipeline.
        print("Matching name")
        if not any(outBaseDir.glob(f"{flowcellName}*")):
            return (flowcellName, flowcellDir, 'aviti')
        # If a matching folder exists, but no flag, start the pipeline:
        elif not any(any(outBaseDir.glob(f"{flowcellName}*/{pattern}")) for pattern in _patterns):
            return (flowcellName, flowcellDir, 'aviti')

    return (None, None, None)


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
    if isinstance(s1, float) or isinstance(s2, float):
        return 0
    if s1 is None or s2 is None:
        return 0
    minl1 = len(s1)
    minl2 = len(s2)
    dist = 0
    for step in range(min([minl1, minl2])):
        if s1[step] != s2[step]:
            dist += 1
    return dist


def joinLis(lis, joinStr=""):
    """
    join a list into a string (without spaces).
    elements are converted to strings.
    """
    return joinStr.join([str(i) for i in lis])


def lenMask(recipe, minl,aviti):
    """
    take length of recipe (runInfo) and length of a barcode and return a mask.
    e.g. 8bp index, 10bp sequenced, returns I8N2
    """
    if recipe-minl > 0:
        return "Y{}N{}".format(int(minl), int(recipe-minl)) if aviti else "I{}N{}".format(int(minl), int(recipe-minl))
    else:
        return "Y{}".format(int(minl)) if aviti else "I{}".format(int(minl))


def P5Seriesret(df):
    if 'index2' in list(df.columns):
        return df['index2']
    else:
        return pd.Series(dtype='float64')


def krakenfqs(IDdir):
    '''
    Returns:
    abspath to kraken report
    [--paired, R1, R2]
    or
    [R1]
    '''
    fqFiles = []
    # sort glob to ensure R1 comes before R2
    for fq in sorted(Path(IDdir).glob("*fastq.gz")):
        if fq.name.endswith('_R1.fastq.gz'):
            fqFiles.append(fq)
        elif fq.name.endswith('_R2.fastq.gz'):
            fqFiles.append(fq)
    krakRep = str(fqFiles[0]).replace(
        '_R1.fastq.gz',
        ''
    ) + '.rep'
    krakRep = krakRep.replace(
        'Project_',
        'FASTQC_Project_'
    )  # output to fastqc folder, not project.
    if len(fqFiles) == 1:
        return (
            krakRep, [str(fqFiles[0])]
        )
    elif len(fqFiles) == 2:
        return (
            krakRep, ['--paired', str(fqFiles[0]), str(fqFiles[1])]
        )


def retBCstr(ser, returnHeader=False):
    if returnHeader:
        if 'index2' in list(ser.index):
            return ("P7\tP5")
        else:
            return ("P7")
    if 'index2' in list(ser.index):
        return (
            '\t'.join(
                [str(ser['index']), str(ser['index2'])]
            )
        )
    elif 'index' in list(ser.index):
        return (str(ser['index']))
    else:
        return ("nan")


def retIxtype(ser, returnHeader=False):
    if returnHeader:
        if 'I5_Index_ID' in list(ser.index):
            return ("P7type\tP5type")
        else:
            return ("P7type")
    if 'I7_Index_ID' in list(ser.index) and 'I5_Index_ID' in list(ser.index):
        return '\t'.join(
            [str(ser['I7_Index_ID']), str(ser['I5_Index_ID'])]
        )
    elif 'I7_Index_ID' in list(ser.index):
        return str(ser['I7_Index_ID'])
    else:
        return 'NA'


def retMean_perc_Q(ser, returnHeader=False, qtype='meanQ'):
    if qtype not in ser:
        if returnHeader:
            return ('meanQ', 'NA')
        else:
            return ('NA')
    if pd.isna(ser[qtype]):
        return ('NA')
    meanQstr = str(ser[qtype])
    headers = []
    Reads = []
    for read in meanQstr.split(','):
        key = read.split(':')[0]
        val = round(float(read.split(':')[1]), 0)
        if 'I' not in key:
            headers.append('R' + str(key) + '_' + qtype)
        else:
            headers.append(str(key) + '_' + qtype)
        if qtype != 'meanQ':
            Reads.append(str(val) + '%')
        else:
            Reads.append(str(val))
    if returnHeader:
        return ('\t'.join(headers), '\t'.join(Reads))
    else:
        return ('\t'.join(Reads))


def formatSeqRecipe(seqRecipe):
    '''
    SeqRecipe is a dictionary of form:
    {key:['Y', len], ...}
    We want to return a string combining key and lens.
    with key being Read1, Read2, Index1, Index2
    '''
    if not seqRecipe:
        return "Unknown"
    retStr = ""
    for key in seqRecipe:
        retStr += "{}:{}; ".format(key, seqRecipe[key][1])
    return (retStr[:-2])


def formatMisMatches(mmDic):
    '''
    mmDic is a dictionary of form:
    {BarcodeMismatchesIndex1: int, BarcodeMismatchesIndex2: int}
    We want to return a string combining key and val.
    '''
    retStr = ""
    for key in mmDic:
        retStr += "{}:{}, ".format(key, mmDic[key])
    return (retStr[:-2])


'''
Uncomment the shutil

def fetchLatestSeqDir(PIpath, seqDir):
    # fetch sorted sequence_data directories ascending
    seqDirs = sorted(
        glob.glob(
            os.path.join(PIpath, seqDir + "*")
        )
    )
    spacedict = {}
    for _s in seqDirs:
        # total, used, free bytes
        _t, _u, _f = shutil.disk_usage(_s)
        spacedict[_s] = _u/_t
        if _u/_t < 0.9:
            return (_s)
    logging.info(spacedict)
    logging.critical(
        "No seq_data dir for {} found with space. Exiting.".format(PIpath)
    )
    sys.exit()
'''


def fetchLatestSeqDir(config, PI):
    '''
    Fetch the latest sequencing_data dir in the PI directory
    '''
    PIpath = Path(config['Dirs']['piDir'], PI)
    seqDir = config['Internals']['seqDir']
    seqDirNum = 0
    for dirs in PIpath.iterdir():
        if seqDir in dirs.name:
            seqDirStrip = dirs.name.replace('sequencing_data', '')
            if seqDirStrip != '':
                if int(seqDirStrip) > seqDirNum:
                    seqDirNum = int(seqDirStrip)
    if seqDirNum == 0:
        return Path(PIpath, 'sequencing_data')
    else:
        return Path(PIpath) / 'sequencing_data{}'.format(str(seqDirNum))


def umlautDestroyer(germanWord):
    '''
    Destroy umlauts.
    Illumina destroys: Förtsch -> Fortsch.
    We do too.
    Only exception is ß, which goes to ss.
    Add in a replacement for spaces as well.
    '''

    _u = 'ü'.encode()
    _ef = 'é'.encode()
    _er = 'è'.encode()
    _af = 'á'.encode()
    _ar = 'à'.encode()
    _U = 'Ü'.encode()
    _a = 'ä'.encode()
    _A = 'Ä'.encode()
    _o = 'ö'.encode()
    _O = 'Ö'.encode()
    _ss = 'ß'.encode()
    _apstrf = "'".encode()

    _string = germanWord.encode()
    _string = _string.replace(_u, b'u')
    _string = _string.replace(_ef, b'e')
    _string = _string.replace(_er, b'e')
    _string = _string.replace(_af, b'a')
    _string = _string.replace(_ar, b'a')
    _string = _string.replace(_U, b'U')
    _string = _string.replace(_a, b'a')
    _string = _string.replace(_A, b'A')
    _string = _string.replace(_o, b'o')
    _string = _string.replace(_O, b'O')
    _string = _string.replace(_ss, b'ss')
    _string = _string.replace(_apstrf, b'')
    return (_string.decode('utf-8').replace(' ', ''))


def multiQC_yaml(flowcell, project, laneFolder):
    '''
    This function creates:
     - config yaml, containing appropriate header information
     - data string adding gen stats
     - data string containing our old seqreport statistics.
    Keep in mind we delete these after running mqc
    '''
    logging.info("Postmux - multiqc yaml creation")
    ssDic = flowcell.sampleSheet.ssDic[laneFolder.name]
    ssdf = ssDic['sampleSheet'][ssDic['sampleSheet']['Sample_Project'] == project]
    # data string genstats
    mqcData = "# format: 'tsv'\n"
    mqcData += "# plot_type: 'generalstats'\n"
    mqcData += "# pconfig: \n"
    mqcData += "Sample_Name\tSample_ID\tRequested\n"
    reqDict = {}
    for sample in list(ssdf['Sample_Name'].unique()):
        sampleID = ssdf[ssdf['Sample_Name'] == sample]['Sample_ID'].values[0]
        if ssdf[ssdf['Sample_Name'] == sample]['reqDepth'].values[0] == 'NA':
            reqDepth = 'NA'
        else:
            reqDepth = float(
                ssdf[ssdf['Sample_Name'] == sample]['reqDepth'].values[0]
            )
            reqDepth = round(reqDepth/1000000, 2)
        # 
        sampleLis = sorted(laneFolder.glob(f"*/Sample_{sampleID}/*[IR][12].fastq.gz"))
        for fqFile in sampleLis:
            # basename, with 'fastq.gz'
            reqDict[fqFile.with_suffix('').with_suffix('').name] = [sampleID, reqDepth]
    for sample in reqDict:
        mqcData += f"{sample}\t{reqDict[sample][0]}\t{reqDict[sample][1]}\n"
    # seqreport Data
    seqreportData = ""
    for index, row in ssdf.iterrows():
        if seqreportData == "":
            meanQ_headers, Meanq = retMean_perc_Q(row, returnHeader=True)
            percq30_headers, perc30 = retMean_perc_Q(
                row, returnHeader=True, qtype='percQ30'
            )
            seqreportData += \
                "\tSample ID\tLane\t{}\t{}\n".format(
                    meanQ_headers,
                    percq30_headers
                )
            seqreportData += "{}\t{}\t{}\t{}\t{}\n".format(
                row['Sample_Name'],
                row['Sample_ID'],
                "L" + str(row['Lane']),
                Meanq,
                perc30
            )
        else:
            Meanq = retMean_perc_Q(row)
            perc30 = retMean_perc_Q(row, qtype='percQ30')
            seqreportData += "{}\t{}\t{}\t{}\t{}\n".format(
                row['Sample_Name'],
                row['Sample_ID'],
                "L" + str(row['Lane']),
                Meanq,
                perc30
            )

    # Index stats.
    indexreportData = ""
    # indexreportData = "\tSample ID\tBarcodes\tBarcode types\n"
    for index, row in ssdf.iterrows():
        if indexreportData == "":
            indexreportData += "\tSample ID\t{}\t{}\n".format(
                retBCstr(row, returnHeader=True),
                retIxtype(row, returnHeader=True)
            )
        indexreportData += "{}\t{}\t{}\t{}\n".format(
            row['Sample_Name'],
            row['Sample_ID'],
            retBCstr(row),
            retIxtype(row)
        )

    # config yaml
    # libraryTypes
    libTypes = ', '.join(list(
        ssdf['Library_Type'].fillna('None').unique()
    ))
    # indexTypes
    ixTypes = ', '.join(list(
        ssdf["indexType"].fillna('None').unique()
    ))
    # Protocols
    protTypes = ', '.join(list(
        ssdf["Description"].fillna('None').unique()
    ))
    # Organisms
    orgs = ', '.join(list(
        ssdf["Organism"].str[0].fillna('None').unique()
    ))
    # Resequencing runs are screwed up (e.g. don't contain the samples)
    # switch total requested to NA
    try:
        sumReqRound = str(
            round((ssdf['reqDepth'].sum())/1000000, 0)
        )
    except TypeError:
        sumReqRound = 'NA'

    if flowcell.sequencer == 'aviti':
        _demuxver = {'bases2fastq': flowcell.config['softwareVers']['bases2fastq']}
    else:
        _demuxver = {'bclconvert': flowcell.config['softwareVers']['bclconvert']}
    mqcyml = {
        "title": project,
        "custom_logo": flowcell.config["misc"]["mpiImg"],
        "custom_logo_url": "https://www.ie-freiburg.mpg.de/",
        "custom_logo_title": "MPI-IE",
        "show_analysis_paths": False,
        "show_analysis_time": False,
        "fastqscreen_simpleplot": False,
        "log_filesize_limit": 2000000000,
        "report_header_info": [
            {"Contact E-mail": flowcell.config["communication"]["bioinfoCore"]},
            {"Flowcell": flowcell.name},
            {"Sequencer Type": flowcell.sequencer},
            {"Read Lengths": formatSeqRecipe(flowcell.seqRecipe)},
            {"Demux. Mask": ssDic["mask"]},
            {"Mismatches": formatMisMatches(ssDic["mismatch"])},
            {"dissectBCL version": f'{version("dissectBCL")}'},
            _demuxver,
            {"Library Type": libTypes},
            {"Library Protocol": protTypes},
            {"Index Type": ixTypes},
            {"Organism": orgs},
            {"Requested reads": sumReqRound},
            {"Received reads": str(
                round(
                    (ssdf['gotDepth'].replace(
                        'NA', np.nan
                    ).dropna().sum())/1000000,
                    0
                )
                )}
        ],
        "section_comments": {
            "kraken": flowcell.config["misc"]['krakenExpl']
        }

    }
    return (mqcyml, mqcData, seqreportData, indexreportData)


def stripRights(enduserBase):
    for root in Path(enduserBase).rglob('*'):
        if root.is_dir() or root.is_file():
            os.chmod(root, 0o700)

def getDiskSpace(outputDir):
    '''
    Return space free in GB
    '''
    total, used, free = shutil.disk_usage(outputDir)
    return (total // (2**30), free // (2**30))

def fexUpload(outLane, project, fromA, opas):
    '''
    outLane = 240619_M01358_0047_000000000-LKGP2_lanes_1
    project = Project_xxxx_user_PI
    fromA = from sender (comes from config)
    opas = (path/to/project_xxx_user_PI, path/to/FASTQC_project_xx_user_PI)
    '''
    replaceStatus ='Uploaded'
    tarBall = outLane + '_' + project + '.tar'
    fexList = sp.check_output(
        ['fexsend', '-l', fromA]
    ).decode("utf-8").replace("\n", " ").split(' ')
    if tarBall in fexList:
        logging.info("fakenews - {project} found in fex. Replacing.")
        fexRm = ['fexsend', '-d', tarBall, fromA]
        fexdel = sp.Popen(fexRm)
        fexdel.wait()
        replaceStatus = 'Replaced'
    fexsend = f"tar cf - {opas[0]} {opas[1]} | fexsend -s {tarBall} {fromA}"
    os.system(fexsend)
    return replaceStatus

def sendMqcReports(outPath, tdirs):
    '''
    Ship mqc reports to seqfacdir and bioinfocoredir.
    outPath = /path/to/240619_M01358_0047_000000000-LKGP2_lanes_1
    tdirs = Dirs part from config, contains bioinfoCoreDir and seqFacDir
    '''
    outLane = outPath.name
    yrstr = '20' + outLane[:2]
    BioInfoCoreDir = Path(tdirs['bioinfoCoreDir'])
    seqFacDir = Path(tdirs['seqFacDir']) / f"Sequence_Quality_{yrstr}" / f"Illumina_{yrstr}" / outLane
    seqFacDir.mkdir(parents=True, exist_ok=True)
    for _mq in outPath.glob("*/*multiqc_report.html"):
        sout = seqFacDir / _mq.name
        bout = BioInfoCoreDir / f"{outLane}_{_mq.name}"
        logging.info(f"fakenews - sedMqcReports - seqfac: {sout}")
        logging.info(f"fakenews - sedMqcReports - bioinfo: {bout}")
        shutil.copyfile(_mq, sout)
        shutil.copyfile(_mq, bout)

def matchOptdupsReqs(optDups, ssdf):
    '''
    Takes a nested list (optDups) with:
    [
        project,
        sampleID,
        sampleName,
        optical_dups,
    ]
    Matches sampleID with ssdf, and gets gotten / req reads.
    returns a nested list including req/got & got
    this list is sorted by sampleID.
    '''
    _optDups = []
    for lis in optDups:
        sampleID = lis[1]
        sampleName = lis[2]
        req = ssdf[
            ssdf['Sample_ID'] == sampleID
        ]['reqDepth'].values
        got = ssdf[
            ssdf['Sample_ID'] == sampleID
        ]['gotDepth'].values
        # At this stage, 'PhiX' projects in aviti yield a list
        if isinstance(got, np.ndarray) and len(got) > 1:
            if len(set(got)) != 1:
                _failmsg = f"Received multiple depths for single sample {sampleID}, {sampleName}"
                logging.critical(_failmsg)
                sys.exit(_failmsg)
            got = got[0]
            if len(req) == 1:
                _failmsg = f"Multiple depths received for {sampleID}, {sampleName}, but only one reqDepth. Is this listly wrongly in parkour ?"
                logging.critical(_failmsg)
                sys.exit(_failmsg)
            req = 2000000 # Assume 2 million phiX reads ~= 2% for 800M flow cell
        reqvgot = float(got/req)
        # isnull if sample is omitted from demuxsheet but in parkour.
        if pd.isnull(got):
            _optDups.append(
                [
                    lis[0],
                    sampleID,
                    sampleName,
                    lis[3],
                    0,
                    0
                ]  # fill in zeroes
            )
        else:
            _optDups.append(
                [
                    lis[0],
                    sampleID,
                    sampleName,
                    lis[3],
                    round(reqvgot, 2),
                    int(got)
                ]
            )
    return (sorted(_optDups, key=lambda x: x[1]))
