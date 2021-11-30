from dissectBCL.fakeNews import log
from dissectBCL.misc import joinLis, hamming, lenMask
from dissectBCL.misc import bclConvPipeLogger, P5Seriesret
from itertools import combinations
import os
from subprocess import Popen, PIPE
import sys


def misMatcher(P7s, P5s):
    """
    return the number of mismatches allowed in demux.
    [0, 1 or 2]
    """
    mmDic = {}
    hammings = []
    for comb in combinations(P7s, 2):
        hammings.append(hamming(comb[0], comb[1]))
    if min(hammings) > 2 and min(hammings) <= 4:
        mmDic['BarcodeMismatchesIndex1'] = 1
    elif min(hammings) > 4:
        mmDic['BarcodeMismatchesIndex1'] = 2
    else:
        mmDic['BarcodeMismatchesIndex1'] = 0
    if not P5s.empty:
        hammings = []
        for comb in combinations(P5s, 2):
            hammings.append(hamming(comb[0], comb[1]))
        if min(hammings) == 3:
            mmDic['BarcodeMismatchesIndex2'] = 1
        elif min(hammings) > 3:
            mmDic['BarcodeMismatchesIndex2'] = 2
        else:
            mmDic['BarcodeMismatchesIndex2'] = 0
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
    if 'Index2' in seqRecipe:
        P5seq = True
        recipeP5 = seqRecipe['Index2'][1]
    if 'Read2' in seqRecipe:
        PE = True
    # Capture NuGEN Ovation Solo or scATAC
    if 'indexType' in list(sampleSheetDF.columns):
        # Nugen Ovation SOLO
        if any(
            sampleSheetDF['indexType'].str.contains(
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
            return ";".join(mask), dualIx, PE, convertOpts
        # scATAC
        elif any(sampleSheetDF['indexType'].str.contains(
            "ATAC-Seq single cell"
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
            # Read 2
            if PE:
                mask.append(joinLis(seqRecipe['Read2']))
            else:
                log.warning("Single end sequencing.")
            convertOpts = []
            return ";".join(mask), dualIx, PE, convertOpts
        else:  # general case.
            log.info(
                "{} is not special. Setting default mask".format(outputFolder)
            )
            # Read 1
            mask.append(joinLis(seqRecipe['Read1']))
            # Index 1 (sample barcode)
            mask.append(lenMask(recipeP7, minP7))
            # Index2
            if dualIx:
                mask.append(lenMask(recipeP5, minP5))
            if not dualIx and P5seq:
                mask.append("N{}".format(recipeP5))
            # Read2
            if PE:
                mask.append(joinLis(seqRecipe['Read2']))
            convertOpts = []
            return ";".join(mask), dualIx, PE, convertOpts


def prepConvert(flowcell, sampleSheet):
    log.warning("PreFQ module")
    log.info("determine masking")
    for outputFolder in sampleSheet.ssDic:
        sampleSheet.ssDic[outputFolder]['mask'], \
            sampleSheet.ssDic[outputFolder]['dualIx'], \
            sampleSheet.ssDic[outputFolder]['PE'], \
            sampleSheet.ssDic[outputFolder]['convertOpts'] = detMask(
                flowcell.seqRecipe,
                sampleSheet.ssDic[outputFolder]['sampleSheet'],
                outputFolder,
            )
        sampleSheet.ssDic[outputFolder]['mismatch'] = misMatcher(
            sampleSheet.ssDic[outputFolder]['sampleSheet']['index'],
            P5Seriesret(
                sampleSheet.ssDic[outputFolder]['sampleSheet']
            )
        )
    return sampleSheet


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
        demuxSheetLines.append(
            "BarcodeMismatchesIndex2,{},,".format(
                ssDic['mismatch']['BarcodeMismatchesIndex2']
            )
        )
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
                "demuxSheet for {} already exists.".format(outLane)
            )
        # Run bcl-convert
        bclOpts = [
            config['software']['bclconvert'],
            '--output-directory', outputFolder,
            '--force',
            '--bcl-input-directory', flowcell.bclPath,
            '--sample-sheet', demuxOut,
            '--bcl-num-conversion-threads', "20",
            '--bcl-num-compression-threads', "20",
            "--bcl-sampleproject-subdirectories", "true"
        ]
        log.info("Starting BCLConvert")
        log.info(" ".join(bclOpts))
        bclRunner = Popen(
            bclOpts,
            stdout=PIPE
        )
        with bclRunner.stdout:
            bclConvPipeLogger(bclRunner.stdout)
        exitcode = bclRunner.wait()
        if exitcode == 0:
            log.info("bclConvert exit {}".format(exitcode))
        else:
            log.critical("bclConvert exit {}".format(exitcode))
            sys.exit(1)