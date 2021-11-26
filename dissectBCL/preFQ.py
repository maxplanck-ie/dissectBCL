from dissectBCL.fakeNews import log
from dissectBCL.misc import joinLis, hamming
from itertools import combinations
import os
import pandas as pd


'''
Cycle, masking, mismatching and sampleSheetClass changes.
'''

def misMatcher(P7s, P5s=None):
    """
    return the number of mismatches allowed in demux.
    [0, 1 or 2]
    """
    mmDic = {}
    mm = 0
    hammings = []
    for comb in combinations(P7s, 2):
        hammings.append(hamming(comb[0], comb[1]))
    if min(hammings) == 3:
        mmDic['BarcodeMismatchesIndex1'] = 1
    elif min(hammings) > 3:
        mmDic['BarcodeMismatchesIndex1'] = 2
    else:
        mmDic['BarcodeMismatchesIndex1'] = 0
    if P5s:
        mm = 0
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
    Determine the actual masking based on the runInfo stats ('seqRecipe') and the information in the 
    """
    log.info("determine masking for {}".format(outputFolder))
    mask = []
    dualIx = False
    ## Capture both cases with either NuGEN Ovation Solo (P7 UMI), or scATAC (P5 UMI)
    if 'indexType' in list(sampleSheetDF.columns):
        if any(sampleSheetDF['indexType'].str.contains("NuGEN Ovation SoLo RNA-Seq System")):
            log.info("NuGEN Ovation SoLo found for {}.".format(outputFolder))
            ########## Add read 1 ##########
            mask.append(joinLis(seqRecipe['Read1']))
            ########## Add index 1 ##########
            # Index1 contains the UMI after the barcode.
            # Smallest length of index in this lane/flowcell:
            minP7 = sampleSheetDF['index'].str.len().min()
            recipeP7 = seqRecipe['Index1'][1]
            if recipeP7-minP7 > 0:
                mask.append("I{}U{}".format(minP7, recipeP7-minP7))
            else:
                log.warning("NuGEN Ovation solo detected but Index read length == P7, no UMI will be written!")
                mask.append("I{}".format(minP7))
            if "index_2" in sampleSheetDF: #If no index2 they were all NaN in sampleSheet.
                dualIx = True
                minP5 = sampleSheetDF['index_2'].str.len().min()
                recipeP5 = seqRecipe['Index2'][1]
                if recipeP5-minP5 > 0:
                    mask.append("I{}N{}".format(minP5, recipeP5-minP5))
                else:
                    mask.append("I{}".format(minP5, recipeP5-minP5))
            if 'Read2' in seqRecipe:
                mask.append(joinLis(seqRecipe['Read2']))
            return ";".join(mask), dualIx


def prepConvert(flowcell, sampleSheet):
    log.warning("PreFQ module")
    log.info("determine masking")
    for outputFolder in sampleSheet.ssDic:
        sampleSheet.ssDic[outputFolder]['mask'], sampleSheet.ssDic[outputFolder]['dualIx'] = detMask(
                flowcell.seqRecipe,
                sampleSheet.ssDic[outputFolder]['sampleSheet'],
                outputFolder,
            )
        sampleSheet.ssDic[outputFolder]['mismatch'] = misMatcher(
            sampleSheet.ssDic[outputFolder]['sampleSheet']['index'],
            sampleSheet.ssDic[outputFolder]['sampleSheet']['index2'] if 'index2' in sampleSheet.ssDic[outputFolder]['sampleSheet'] else None
        )
    return sampleSheet

#demux(sampleSheet.ssDic[outLane]['sampleSheet'], sampleSheet.ssDic[outLane]['mask'], outFolder)
def writeDemuxSheet(demuxOut, ssDic):
    demuxSheetLines = []
    # Header first.
    demuxSheetLines.append("[Header],,,")
    demuxSheetLines.append("FileFormatVersion,2,,")
    demuxSheetLines.append(",,,")
    demuxSheetLines.append("[BCLConvert_Settings],,,")
    demuxSheetLines.append("BarcodeMismatchesIndex1,{},,".format(ssDic['mismatch']['BarcodeMismatchesIndex1']))
    if ssDic['dualIx']:
        demuxSheetLines.append("BarcodeMismatchesIndex2,{},,".format(ssDic['mismatch']['BarcodeMismatchesIndex2']))
    demuxSheetLines.append("OverrideCycles,{},,".format(ssDic['mask']))
    demuxSheetLines.append(",,,")
    if ssDic['dualIx']:
        demuxSheetLines.append("[BCLConvert_Data],,,,,,")
        demuxSheetLines.append("Lane,Sample_ID,index,index2,Sample_Project")
        for index, row in ssDic['sampleSheet'].iterrows():
            demuxSheetLines.append(joinLis([row['Lane'], row['Sample_ID'], row['index'], row['index2'], row['Sample_Project']], joinStr = ","))
    else:
        demuxSheetLines.append("[BCLConvert_Data],,,,,")
        demuxSheetLines.append("Lane,Sample_ID,index,Sample_Project")
        for index, row in ssDic['sampleSheet'].iterrows():
            demuxSheetLines.append(joinLis([row['Lane'], row['Sample_ID'], row['index'], row['Sample_Project']], joinStr=","))
    log.info(demuxSheetLines)
    if not os.path.exists(demuxOut):
        log.info("No demuxSheet found, writing it out")
        with open(demuxOut, 'w') as f:
            for line in demuxSheetLines:
                f.write('{}\n'.format(line))
    else:
        log.info("demuxSheet already exists, not overwriting it.")


def demux(sampleSheet, outBase):
    log.warning("Demux module")
    for outLane in sampleSheet.ssDic:
        log.info("Demuxing {}".format(outLane))
        # Check outDir
        outputFolder = os.path.join(outBase, outLane)
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
            log.info("{} created.".format(outputFolder))
        else:
            log.info("{} already exists. Moving on.".format(outputFolder))
        # Write the demuxSheet in the outputfolder
        demuxOut = os.path.join(outputFolder, "demuxSheet.csv")
        writeDemuxSheet(demuxOut, sampleSheet.ssDic[outLane])
    

    


