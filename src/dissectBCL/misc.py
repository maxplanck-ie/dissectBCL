import os
import configparser
import xml.etree.ElementTree as ET
import glob
from dissectBCL.logger import log
import pandas as pd


def getConf(configfile):
    config = configparser.ConfigParser()
    log.info("Writing configfile from {}".format(configfile))
    config.read(configfile)
    return config


def getNewFlowCell(config, fPath=None):
    # If there is a fPath set, just return that.
    if fPath:
        flowcellName = fPath.split('/')[-1]
        flowcellDir = fPath
        return (flowcellName, flowcellDir)
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
        # If there is no folder matching the flowcell name, start the pipeline.
        if not glob.glob(
            os.path.join(outBaseDir, flowcellName) + "*"
        ):
            return flowcellName, flowcellDir
        # If a matching folder exists, but no flag, start the pipeline:
        elif not glob.glob(
            os.path.join(outBaseDir, flowcellName + '*', 'communication.done')
        ) and not glob.glob(
            os.path.join(outBaseDir, flowcellName + '*', 'fastq.made')  # bfq
        ):
            return (flowcellName, flowcellDir)
    return (None, None)


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


def lenMask(recipe, minl):
    """
    take length of recipe (runInfo) and length of a barcode and return a mask.
    e.g. 8bp index, 10bp sequenced, returns I8N2
    """
    print("recip:{} minl:{}".format(recipe, minl))
    if recipe-minl > 0:
        return "I{}N{}".format(int(minl), int(recipe-minl))
    else:
        return "I{}".format(int(minl))


def P5Seriesret(df):
    if 'index2' in list(df.columns):
        return df['index2']
    else:
        return pd.Series(dtype='float64')


def screenFqFetcher(IDdir):
    """
    Return what fastq file should be used in the fastq screen.
    Prioritize R2 > R1
    R3 no longer needed since we don't produce it anymore.
    """
    fqFiles = glob.glob(
        os.path.join(
            IDdir,
            "*fastq.gz"
        )
    )
    for substr in ["_R2.fastq.gz", "_R1.fastq.gz"]:
        hit = [s for s in fqFiles if substr in s and 'optical' not in s]
        if hit:
            return hit[0]


def moveOptDup(laneFolder):
    for txt in glob.glob(
        os.path.join(
            laneFolder,
            '*',
            '*',
            '*duplicate.txt'
        )
    ):
        # Field -3 == project folder
        # escape those already in a fastqc folder (reruns)
        if 'FASTQC' not in txt:
            pathLis = txt.split('/')
            pathLis[-3] = 'FASTQC_' + pathLis[-3]
            ofile = "/".join(pathLis)
            ofile.replace('duplicate.txt', 'opticalduplicates.txt')
            os.rename(txt, ofile)


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
    if str(ser[qtype]) == 'NA':
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


def fetchLatestSeqDir(PIpath, seqDir):
    seqDirNum = 0
    for dir in os.listdir(PIpath):
        if seqDir in dir and dir.replace(seqDir, ''):
            seqDirNum = int(dir[-1])
    if seqDirNum == 0:
        return (os.path.join(PIpath, seqDir))
    else:
        return (os.path.join(PIpath, seqDir + str(seqDirNum)))


def umlautDestroyer(germanWord):
    '''
    Destroy umlauts.
    Illumina destroys: Förtsch -> Fortsch.
    We do too.
    Only exception is ß, which goes to ss.
    Add in a replacement for spaces as well.
    '''

    _u = 'ü'.encode()
    _U = 'Ü'.encode()
    _a = 'ä'.encode()
    _A = 'Ä'.encode()
    _o = 'ö'.encode()
    _O = 'Ö'.encode()
    _ss = 'ß'.encode()

    _string = germanWord.encode()
    _string = _string.replace(_u, b'u')
    _string = _string.replace(_U, b'U')
    _string = _string.replace(_a, b'a')
    _string = _string.replace(_A, b'A')
    _string = _string.replace(_o, b'o')
    _string = _string.replace(_O, b'O')
    _string = _string.replace(_ss, b'ss')
    return (_string.decode('utf-8').replace(' ', ''))
