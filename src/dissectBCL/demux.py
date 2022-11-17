from dissectBCL.logger import log
from dissectBCL.misc import joinLis, hamming, lenMask
from dissectBCL.misc import P5Seriesret, matchingSheets
from dissectBCL.fakeNews import mailHome
from itertools import combinations
import os
from subprocess import Popen, PIPE
import sys
from pathlib import Path
import pandas as pd
import numpy as np


def hamming2Mismatch(minVal):
    if minVal > 2 and minVal <= 4:
        return 1
    elif minVal > 4:
        return 2
    else:
        return 0


def misMatcher(P7s, P5s):
    """
    return the number of mismatches allowed in demux.
    [0, 1 or 2]
    """
    mmDic = {}
    hammings = []
    for comb in combinations(P7s, 2):
        hammings.append(hamming(comb[0], comb[1]))
    mmDic['BarcodeMismatchesIndex1'] = hamming2Mismatch(min(hammings))
    if not P5s.empty and not P5s.isnull().all():
        hammings = []
        for comb in combinations(P5s, 2):
            hammings.append(hamming(comb[0], comb[1]))
        mmDic['BarcodeMismatchesIndex2'] = hamming2Mismatch(min(hammings))
    return mmDic


def detMask(seqRecipe, sampleSheetDF, outputFolder):
    """
    Determine the actual masking.
    Based on:
        - seqRecipe (RunInfo.xml)
        - sampleSheet
        - parkour info
    """
    log.info("determine masking for {}".format(outputFolder))
    log.info("masking for seqRecipe {}".format(seqRecipe))
    mask = []
    dualIx = False
    PE = False
    P5seq = False
    # Find out the actual index size and how much was sequenced.
    minP7 = sampleSheetDF['index'].str.len().min()
    recipeP7 = seqRecipe['Index1'][1]
    if 'index2' in list(sampleSheetDF.columns):
        dualIx = True
        minP5 = sampleSheetDF['index2'].str.len().min()
    else:
        minP5 = np.NAN
    if 'Index2' in seqRecipe:
        P5seq = True
        recipeP5 = seqRecipe['Index2'][1]
    if 'Read2' in seqRecipe:
        PE = True
    # Capture NuGEN Ovation Solo or scATAC
    if 'indexType' in list(sampleSheetDF.columns):
        log.info("indexType column found.")
        log.info("indexType series:")
        # Nugen Ovation SOLO
        if any(
            sampleSheetDF['indexType'].dropna().str.contains(
                "NuGEN Ovation SoLo RNA-Seq System"
            )
        ):
            log.info("NuGEN Ovation SoLo found for {}.".format(outputFolder))
            # Read 1
            mask.append(joinLis(seqRecipe['Read1']))
            # Index1 (index1 = 8bp index, 8bp UMI)
            if recipeP7-minP7 > 0:
                mask.append("I{}U{}".format(minP7, recipeP7-minP7))
            else:
                log.warning(
                    "NuGEN Ovation solo Index read length == P7!"
                )
                mask.append("I{}".format(minP7))
            # Index 2
            if dualIx:
                mask.append(lenMask(recipeP5, minP5))
            # Read 2
            if PE:
                mask.append(joinLis(seqRecipe['Read2']))
            # set UMI ops.
            convertOpts = ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
            return ";".join(mask), dualIx, PE, convertOpts, None, None
        # scATAC
        elif any(sampleSheetDF['Description'].dropna().str.contains(
            "scATAC-Seq 10xGenomics"
            )
        ):
            log.info("scATAC seq found for {}".format(outputFolder))
            # Read 1
            mask.append(joinLis(seqRecipe['Read1']))
            # Index 1 (sample barcode)
            mask.append(lenMask(recipeP7, minP7))
            # Index2 (10x barcode)
            if 'Index2' in seqRecipe:
                mask.append("U{}".format(recipeP5))
            if dualIx:
                log.warning(
                    "P5 detected, Mixed modalities not processed by default."
                )
                dualIx = False
            # Read 2
            if PE:
                mask.append(joinLis(seqRecipe['Read2']))
            else:
                log.warning("Single end sequencing.")
            convertOpts = ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
            return ";".join(mask), dualIx, PE, convertOpts, None, None
        else:  # general case.
            log.info(
                "{} is not special. Setting default mask".format(outputFolder)
            )
            # Read 1
            mask.append(joinLis(seqRecipe['Read1']))
            # Index 1 (sample barcode)
            mask.append(lenMask(recipeP7, minP7))
            # Index2
            if np.isnan(minP5):
                log.info("P5 is sequenced, but libs in lane are P7 only!")
                dualIx = False
            if dualIx:
                mask.append(lenMask(recipeP5, minP5))
            if not dualIx and P5seq:
                mask.append("N{}".format(recipeP5))
            # Read2
            if PE:
                mask.append(joinLis(seqRecipe['Read2']))
            convertOpts = []
            return ";".join(mask), dualIx, PE, convertOpts, minP5, minP7
    else:
        log.info("parkour failure probably, revert back to what we can.")


def prepConvert(flowcell, sampleSheet, config):
    log.warning("PreFQ module")
    log.info("determine masking")
    for outputFolder in sampleSheet.ssDic:
        # add check for bclconvert also here in addition to demux
        if os.path.exists(
            os.path.join(
                config['Dirs']['outputDir'], 
                outputFolder, 
                'bclconvert.done'
            )
        ):
            continue
        sampleSheet.ssDic[outputFolder]['mask'], \
            sampleSheet.ssDic[outputFolder]['dualIx'], \
            sampleSheet.ssDic[outputFolder]['PE'], \
            sampleSheet.ssDic[outputFolder]['convertOpts'], \
            minP5, \
            minP7 = detMask(
                flowcell.seqRecipe,
                sampleSheet.ssDic[outputFolder]['sampleSheet'],
                outputFolder,
            )
        # extra check to make sure all our indices are of equal size!
        if minP7:
            sampleSheet.ssDic[
                outputFolder
            ]['sampleSheet']['index'] = sampleSheet.ssDic[
                outputFolder
            ]['sampleSheet']['index'].str[:minP7]
        if minP5:
            if not np.isnan(minP5):
                sampleSheet.ssDic[
                    outputFolder
                ]['sampleSheet']['index2'] = sampleSheet.ssDic[
                    outputFolder
                ]['sampleSheet']['index2'].str[:minP5]
        sampleSheet.ssDic[outputFolder]['mismatch'] = misMatcher(
            sampleSheet.ssDic[outputFolder]['sampleSheet']['index'],
            P5Seriesret(
                sampleSheet.ssDic[outputFolder]['sampleSheet']
            )
        )
    log.info("mask in sampleSheet updated.")
    return (0)


def writeDemuxSheet(demuxOut, ssDic, laneSplitStatus):
    demuxSheetLines = []
    # Header first.
    demuxSheetLines.append("[Header],,,")
    demuxSheetLines.append("FileFormatVersion,2,,")
    demuxSheetLines.append(",,,")
    demuxSheetLines.append("[BCLConvert_Settings],,,")
    demuxSheetLines.append(
        "BarcodeMismatchesIndex1,{},,".format(
            ssDic['mismatch']['BarcodeMismatchesIndex1']
        )
    )
    if ssDic['dualIx']:
        # crash when single + dual index mixed but still want to write ss.
        if 'BarcodeMismatchesIndex2' in ssDic['mismatch']:
            demuxSheetLines.append(
                "BarcodeMismatchesIndex2,{},,".format(
                    ssDic['mismatch']['BarcodeMismatchesIndex2']
                )
            )
        else:
            log.warning("dualIx set, but no mismatch returned. Overriding.")
            ssDic['dualIx'] = False
    demuxSheetLines.append("OverrideCycles,{},,".format(ssDic['mask']))
    if len(ssDic['convertOpts']) > 0:
        for opts in ssDic['convertOpts']:
            demuxSheetLines.append(opts)
    demuxSheetLines.append(",,,")
    if ssDic['dualIx']:
        demuxSheetLines.append("[BCLConvert_Data],,,,,,")
        if laneSplitStatus:
            demuxSheetLines.append(
                "Lane,Sample_ID,index,index2,Sample_Project"
            )
            for index, row in ssDic['sampleSheet'].iterrows():
                demuxSheetLines.append(
                    joinLis(
                        [
                            row['Lane'],
                            row['Sample_ID'],
                            row['index'],
                            row['index2'],
                            row['Sample_Project']
                        ],
                        joinStr=","
                    )
                )
        else:
            demuxSheetLines.append(
                "Sample_ID,index,index2,Sample_Project"
            )
            for index, row in ssDic['sampleSheet'].iterrows():
                demuxSheetLines.append(
                    joinLis(
                        [
                            row['Sample_ID'],
                            row['index'],
                            row['index2'],
                            row['Sample_Project']
                        ],
                        joinStr=","
                    )
                )
    else:
        demuxSheetLines.append("[BCLConvert_Data],,,,,")
        if laneSplitStatus:
            demuxSheetLines.append(
                "Lane,Sample_ID,index,Sample_Project"
            )
            for index, row in ssDic['sampleSheet'].iterrows():
                demuxSheetLines.append(
                    joinLis(
                        [
                            row['Lane'],
                            row['Sample_ID'],
                            row['index'],
                            row['Sample_Project']
                        ],
                        joinStr=","
                    )
                )
        else:
            demuxSheetLines.append(
                "Sample_ID,index,Sample_Project"
            )
            for index, row in ssDic['sampleSheet'].iterrows():
                demuxSheetLines.append(
                    joinLis(
                        [
                            row['Sample_ID'],
                            row['index'],
                            row['Sample_Project']
                        ],
                        joinStr=","
                    )
                )
    with open(demuxOut, 'w') as f:
        for line in demuxSheetLines:
            f.write('{}\n'.format(line))


def readDemuxSheet(demuxSheet):
    '''
    In case of manual intervention.
    We want to have the correct info in reports / emails.
    Thus, we reparse the demuxSheet just to be sure.
    We check for:
     - 'mask' (overridecycles)
     - indices used.
    '''
    with open(demuxSheet) as f:
        sampleStatus = False
        nesLis = []
        for line in f:
            line = line.strip()
            if line.startswith('OverrideCycles'):
                mask = line.replace(
                    'OverrideCycles', ''
                ).replace(
                    ',', ''
                )
            if sampleStatus:
                nesLis.append(line.split(','))
            if line.startswith('[BCLConvert_Data]') or \
                line.startswith('[Data]'):
                sampleStatus = True
    df = pd.DataFrame(
        nesLis[1:],
        columns=nesLis[0]
    )
    if 'index2' in list(df.columns):
        dualIx = True
    else:
        dualIx = False
    
    # test if mask has been defined
    try: 
        mask
    except NameError: 
        mask = None

    return (mask, df, dualIx)


def parseStats(outputFolder, ssdf):
    QCmetFile = os.path.join(
        outputFolder,
        'Reports',
        'Quality_Metrics.csv'
    )
    DemuxmetFile = os.path.join(
        outputFolder,
        'Reports',
        'Demultiplex_Stats.csv'
    )
    QCdf = pd.read_csv(QCmetFile, sep=',', header=0)
    muxdf = pd.read_csv(DemuxmetFile, sep=',', header=0)
    QCmetDic = {}
    for index, row in QCdf.iterrows():
        sampleID = row['SampleID']
        readnum = str(row['ReadNumber'])
        if sampleID not in QCmetDic:
            QCmetDic[sampleID] = {}
        if readnum not in QCmetDic[sampleID]:
            QCmetDic[sampleID][readnum] = [
                float(row['Mean Quality Score (PF)']),
                int(float(row['% Q30']) * 100)
            ]
        else:
            new_QS = (
                QCmetDic[sampleID][readnum][0] +
                float(row['Mean Quality Score (PF)'])
            )/2
            new_Q30 = (
                QCmetDic[sampleID][readnum][1] + float(row['% Q30']) * 100
            )/2
            QCmetDic[sampleID][readnum] = [new_QS, new_Q30]
    muxDic = {}
    for index, row in muxdf.iterrows():
        sampleID = row['SampleID']
        if sampleID not in muxDic:
            muxDic[sampleID] = int(row['# Reads'])
        else:
            muxDic[sampleID] += int(row['# Reads'])
    MetrixDic = {}
    for ID in QCmetDic:
        QCstr = ""
        Perc30str = ""
        for read in QCmetDic[ID]:
            QCstr += "{}:{},".format(read, round(QCmetDic[ID][read][0], 2))
            Perc30str += "{}:{},".format(read, round(QCmetDic[ID][read][1], 2))
        QCstr = QCstr[:-1]
        Perc30str = Perc30str[:-1]
        MetrixDic[ID] = {
            'meanQ': QCstr,
            'percQ30': Perc30str,
            'gotDepth': int(round(muxDic[ID]/1000000, 0))
        }
    MetrixDF = pd.DataFrame(MetrixDic).T
    MetrixDF['Sample_ID'] = MetrixDF.index
    # left join to only get samples already present
    newDF = pd.merge(ssdf, MetrixDF, on='Sample_ID', how='left')
    newDF = newDF[newDF['Sample_ID'] != 'Undetermined']
    return (newDF)


def demux(sampleSheet, flowcell, config):
    log.warning("Demux module")
    for outLane in sampleSheet.ssDic:
        log.info("Demuxing {}".format(outLane))
        # Check outDir
        outputFolder = os.path.join(flowcell.outBaseDir, outLane)
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
            log.info("{} created.".format(outputFolder))
        else:
            log.info("{} already exists. Moving on.".format(outputFolder))
        # Write the demuxSheet in the outputfolder
        demuxOut = os.path.join(outputFolder, "demuxSheet.csv")
        # Don't remake if demuxSheet exist
        if not os.path.exists(demuxOut):
            log.info("Writing demuxSheet for {}".format(outLane))
            writeDemuxSheet(
                demuxOut,
                sampleSheet.ssDic[outLane],
                sampleSheet.laneSplitStatus
            )
        else:
            log.warning(
                "demuxSheet for {} already exists!".format(outLane)
            )
            manual_mask, manual_df, manual_dualIx = readDemuxSheet(demuxOut)
            # if mask is changed, update:
            # Mask
            if 'mask' in sampleSheet.ssDic[outLane] and \
                manual_mask != sampleSheet.ssDic[outLane]['mask']:
                log.info("Mask is changed from {} into {}.".format(
                    sampleSheet.ssDic[outLane]['mask'],
                    manual_mask
                ))
                sampleSheet.ssDic[outLane]['mask'] = manual_mask
            # dualIx status
            if manual_dualIx != sampleSheet.ssDic[outLane]['dualIx']:
                log.info("dualIx is changed from {} into {}.".format(
                    sampleSheet.ssDic[outLane]['dualIx'],
                    manual_dualIx
                ))
            # sampleSheet
            sampleSheet.ssDic[outLane]['sampleSheet'] = matchingSheets(
                sampleSheet.ssDic[outLane]['sampleSheet'],
                manual_df
            )
        # Don't run bcl-convert if we have the touched flag.
        if not os.path.exists(
            os.path.join(outputFolder, 'bclconvert.done')
        ):
            # Run bcl-convert
            bclOpts = [
                config['software']['bclconvert'],
                '--output-directory', outputFolder,
                '--force',
                '--bcl-input-directory', flowcell.bclPath,
                '--sample-sheet', demuxOut,
                '--bcl-num-conversion-threads', "20",
                '--bcl-num-compression-threads', "20",
                "--bcl-sampleproject-subdirectories", "true",
            ]
            if not sampleSheet.laneSplitStatus:
                bclOpts.append('--no-lane-splitting')
                bclOpts.append('true')
            log.info("Starting BCLConvert")
            log.info(" ".join(bclOpts))
            bclRunner = Popen(
                bclOpts,
                stdout=PIPE
            )
            exitcode = bclRunner.wait()
            if exitcode == 0:
                log.info("bclConvert exit {}".format(exitcode))
                Path(
                    os.path.join(outputFolder, 'bclconvert.done')
                ).touch()
            elif exitcode == -6:
                # known bug async, related to I/O lag over network (?).
                # Illumina tech support wasn't sure, ticket still open.
                bclRunner = Popen(
                    bclOpts,
                    stdout=PIPE
                )
                exitcode = bclRunner.wait()
                if exitcode == 0:
                    log.info(
                        "bclConvert exit {} after second try.".format(
                            exitcode
                        )
                    )
                    Path(
                        os.path.join(outputFolder, 'bclconvert.done')
                    ).touch()
                else:
                    log.critical("bclConvert exit {}".format(exitcode))
                    mailHome(
                        outLane,
                        'BCL-convert exit {}. Investigate.'.format(
                            exitcode
                        ),
                        config,
                        toCore=True
                    )
                    sys.exit(1)
            else:
                log.critical("bclConvert  exit {}".format(exitcode))
                mailHome(
                        outLane,
                        'BCL-convert exit {}. Investigate.'.format(
                            exitcode
                        ),
                        config,
                        toCore=True
                    )
                sys.exit(1)
        log.info("Parsing stats for {}".format(outLane))
        sampleSheet.ssDic[outLane]['sampleSheet'] = parseStats(
                    outputFolder,
                    sampleSheet.ssDic[outLane]['sampleSheet']
        )
    return (0)
