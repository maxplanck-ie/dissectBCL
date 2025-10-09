from dissectBCL.misc import joinLis
from dissectBCL.misc import hamming
from dissectBCL.misc import lenMask
from dissectBCL.misc import P5Seriesret
from Bio.Seq import Seq
from itertools import combinations
import logging
import numpy as np
import json
import pandas as pd
from pathlib import Path
import shutil


def misMatcher(P7s, P5s, sequencer):
    """
    return the number of mismatches allowed in demux.
    [0, 1 or 2]

    if P7s and P5s are both empty, return an empty dictionary.
    """
    mmDic = {}
    for i, ix_list in enumerate((P7s, P5s)):
        hammings = []
        if not ix_list.empty and not ix_list.isnull().all():
            for comb in combinations(ix_list, 2):
                hammings.append(hamming(comb[0], comb[1]))
                if sequencer == 'aviti':
                    barcode_mm = 'I{}MismatchThreshold'.format(i + 1)
                else:
                    barcode_mm = 'BarcodeMismatchesIndex{}'.format(i + 1)
                mmDic[barcode_mm] = hamming2Mismatch(min(hammings))
    return mmDic


def hamming2Mismatch(minVal):
    if minVal > 2 and minVal <= 4:
        return 1
    elif minVal > 4:
        return 2
    else:
        return 0


def detMask(seqRecipe, sampleSheetDF, outputFolder, sequencer):
    """
    Determine the actual masking.
    Based on:
        - seqRecipe (RunInfo.xml)
        - sampleSheet
        - parkour info

    Special masking care for
     - scATAC (UMI/CB in P5)
     - nugen ovationsolo (P7 = index+UMI)
    """
    logging.info("determine masking for {}".format(outputFolder))
    logging.info("masking for seqRecipe {}".format(seqRecipe))

    scATACl = [
        "scATAC-Seq 10xGenomics",
        "NextGEM_Multiome_ATAC",
        "Next GEM Single Cell ATAC",
        "MUXscATAC-seq v3"
    ]
    # initialize variables
    if sequencer == 'aviti':
        mask = {
            'R1FastQMask': '' ,
            'R2FastQMask': '',
            'I1Mask': '',
            'I2Mask': ''
            }
    else:    
        mask = []
    dualIx = False
    PE = False
    P5seq = False
    convertOpts = []

    # set initial values
    if 'Index2' in seqRecipe and seqRecipe['Index2'][1] > 0:
        P5seq = True
        recipeP5 = seqRecipe['Index2'][1]
    else:
        P5seq = False
    

    if 'Read2' in seqRecipe and seqRecipe['Read2'][1] > 0:
        PE = True
    
    # if there is no index in the reads, then set the return values manually
    # TODO a lot of this is code duplication, refactor me!
    if not any(sr.startswith('Index') for sr in seqRecipe.keys()):
        if sequencer == 'aviti':
            mask['R1FastQMask'] = joinLis(seqRecipe['Read1'])
            if not dualIx and P5seq:
                mask['I2Mask'] = "N{}".format(recipeP5)
            if PE:
                mask['R2FastQMask'] = joinLis(seqRecipe['Read2'])
            return mask, dualIx, PE, convertOpts, None, None
        else:
            mask.append(joinLis(seqRecipe['Read1']))

            # Index 1 (sample barcode)
            if not dualIx and P5seq:
                mask.append("N{}".format(recipeP5))

            if PE:
                mask.append(joinLis(seqRecipe['Read2']))

            # minP5 and minP7 are None here
            return ";".join(mask), dualIx, PE, convertOpts, None, None

    # Find out the actual index size and how much was sequenced.
    index1_colname = 'index'
    if sequencer == 'aviti':
        index1_colname = "Index1"
    minP7 = sampleSheetDF[index1_colname].str.len().min()
    recipeP7 = seqRecipe['Index1'][1]
    
    # Since we don't strip the df for nans anymore, get minLength of P5
    # This will equal out to nan if any P5 == nan ?
    index2_colname = "index2"
    if sequencer == 'aviti':
        index2_colname = "Index2"
    if sampleSheetDF[index2_colname].isna().all():
        minP5 = np.nan
    else:
        minP5 = sampleSheetDF[index2_colname].str.len().min()

    # Capture NuGEN Ovation Solo or scATAC
    if 'indexType' in list(sampleSheetDF.columns):
        logging.info("indexType column found.")
        logging.info("indexType series:")
        # Nugen Ovation SOLO
        if any(
            sampleSheetDF['indexType'].dropna().str.contains(
                "NuGEN Ovation SoLo RNA-Seq System"
            )
        ):
            logging.info(
                "NuGEN Ovation SoLo found for {}.".format(outputFolder)
            )
            # Read 1
            if sequencer == 'aviti':
                mask['R1FastQMask'] = joinLis(seqRecipe['Read1'])
            else:
                mask.append(joinLis(seqRecipe['Read1']))
            # Index1 (index1 = 8bp index, 8bp UMI)
            if recipeP7-minP7 > 0:
                if sequencer == 'aviti':
                    mask['I1Mask'] = "Y{}N{}".format(minP7, recipeP7-minP7)
                else:
                    mask.append("I{}U{}".format(minP7, recipeP7-minP7))
            else:
                logging.warning(
                    "NuGEN Ovation solo Index read length == P7!"
                )
                if sequencer == 'aviti':
                    mask['I1Mask'] = "Y{}".format(minP7)
                else:
                    mask.append("I{}".format(minP7))
            # Index 2
            if dualIx:
                if sequencer == 'aviti':
                    mask['I2Mask'] = lenMask(recipeP5, minP5,aviti=True)
                else:
                    mask.append(lenMask(recipeP5, minP5,aviti=False))
            # Read 2
            if PE:
                if sequencer == 'aviti':
                    mask['R2FastqMask'] = joinLis(seqRecipe['Read2'])
                else:
                    mask.append(joinLis(seqRecipe['Read2']))
            # set UMI ops.
            if sequencer == 'aviti':
                convertOpts = ['']
                return mask, dualIx, PE, convertOpts, None, None
            else:
                convertOpts = ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
                return ";".join(mask), dualIx, PE, convertOpts, None, None
        # scATAC
        elif any(
            sampleSheetDF['Description'].dropna().str.contains(
                '|'.join(scATACl)
            )
        ):
            logging.info("scATAC seq found for {}".format(outputFolder))
            # Read 1
            if sequencer == 'aviti':
                mask['R1FastQMask'] = joinLis(seqRecipe['Read1'])
            else:
                mask.append(joinLis(seqRecipe['Read1']))
            # Index 1 (sample barcode)
            if sequencer == 'aviti':
                mask['I1Mask'] = lenMask(recipeP7, minP7,aviti=True)
            else:
                mask.append(lenMask(recipeP7, minP7,aviti=False))
            # Index2 (10x barcode)
            if 'Index2' in seqRecipe:
                if sequencer == 'aviti':
                    if np.isnan(recipeP5):
                        recipeP5 = seqRecipe['Index2'][1]
                    mask['I2Mask'] = "I2:N{}".format(recipeP5)
                    mask['UmiMask'] = "I2:Y{}".format(recipeP5)
                    mask['I1FastQ'] = True
                    mask['I2FastQ'] = True
                    mask['UmiFastQ'] = True
                else:
                    mask.append("U{}".format(recipeP5))
            if dualIx:
                logging.warning(
                    "P5 detected, Mixed modalities not processed by default."
                )
                dualIx = False
            # Read 2
            if PE:
                if sequencer == 'aviti':
                    mask['R2FastQMask'] = joinLis(seqRecipe['Read2'])
                else:
                    mask.append(joinLis(seqRecipe['Read2']))
            else:
                logging.warning("Single end sequencing.")
            if sequencer == 'aviti':
                convertOpts = ['']
                return mask, dualIx, PE, convertOpts, None, None
            else:
                convertOpts = ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
                return ";".join(mask), dualIx, PE, convertOpts, None, None
        else:  # general case.
            logging.info(
                "{} is not special. Setting default mask".format(outputFolder)
            )
            # Read 1
            if sequencer == 'aviti':
                mask['R1FastQMask'] = joinLis(seqRecipe['Read1'])
            else:
                mask.append(joinLis(seqRecipe['Read1']))
            # Index 1 (sample barcode)
            if sequencer == 'aviti':
                mask['I1Mask'] = lenMask(recipeP7, minP7,aviti=True)
            else:
                mask.append(lenMask(recipeP7, minP7,aviti=False))
            # Index2
            if np.isnan(minP5):
                logging.info("P5 is sequenced, but libs in lane are P7 only!")
                dualIx = False
            else:
                # we have P5s in our sampleSheet, dualIx must be true.
                dualIx = True
            if dualIx:
                if sequencer == 'aviti':
                    mask['I2Mask'] = lenMask(recipeP5, minP5,aviti=True)
                else:
                    mask.append(lenMask(recipeP5, minP5,aviti=False))
            if not dualIx and P5seq:
                if sequencer == 'aviti':
                    mask['I2Mask'] = "N{}".format(recipeP5)
                else:
                    mask.append("N{}".format(recipeP5))
            # Read2
            if PE:
                if sequencer == 'aviti':
                    mask['R2FastQMask'] = joinLis(seqRecipe['Read2'])
                else:
                    mask.append(joinLis(seqRecipe['Read2']))
            if sequencer == 'aviti':
                return mask, dualIx, PE, convertOpts, minP5, minP7
            else:
                return ";".join(mask), dualIx, PE, convertOpts, minP5, minP7
    else:
        logging.info("parkour failure probably, revert back to what we can.")


def prepConvert(flowcell, sampleSheet):
    logging.warning("PreFQ module")
    logging.info("determine masking, indices, paired ends, and other options")
    for outputFolder in sampleSheet.ssDic:
        # assign variables for brevity
        ss_dict = sampleSheet.ssDic[outputFolder]
        ss = ss_dict['sampleSheet']

        # determine mask, dualIx, PE, convertOpts, minP5, minP7 from seqRecipe
        (ss_dict['mask'], ss_dict['dualIx'], ss_dict['PE'],
         ss_dict['convertOpts'], minP5, minP7) = detMask(
            flowcell.seqRecipe,
            ss,
            outputFolder,
            flowcell.sequencer

        )

        # extra check to make sure all our indices are of equal size!
        index1_colname = "index"
        index2_colname = "index2"
        if flowcell.sequencer == 'aviti':
            index1_colname = "Index1"
            index2_colname = "Index2"
        for min_ix, ix_str in ((minP5, index1_colname), (minP7, index2_colname)):  #is this correct? isn't index1 P7 and index2 P5 ?
            if min_ix and not np.isnan(min_ix):
                ss[ix_str] = ss[ix_str].str[:min_ix]

        # determine mismatch
        ss_dict['mismatch'] = misMatcher(ss[index1_colname], P5Seriesret(ss),flowcell.sequencer)
    
    logging.info("mask in sampleSheet updated.")
    return (0)


def writeDemuxSheet(demuxOut, ssDic, laneSplitStatus):
    demuxSheetLines = []
    # Header first.
    demuxSheetLines.append("[Header],,,")
    demuxSheetLines.append("FileFormatVersion,2,,")
    demuxSheetLines.append(",,,")
    demuxSheetLines.append("[BCLConvert_Settings],,,")
    if 'mismatch' in ssDic:
        for i in (1, 2):
            bc_str = 'BarcodeMismatchesIndex{}'.format(i)
            if i == 2 and not ssDic['dualIx']:
                continue
            if bc_str in ssDic['mismatch']:
                demuxSheetLines.append(
                    "{},{},,".format(bc_str, ssDic['mismatch'][bc_str])
                )

            # TODO do we want this behavior?
            # log.warning("dualIx set, but no mismatch returned. Overriding.")
            # ssDic['dualIx'] = False

    demuxSheetLines.append("OverrideCycles,{},,".format(ssDic['mask']))
    if len(ssDic['convertOpts']) > 0:
        for opts in ssDic['convertOpts']:
            demuxSheetLines.append(opts)
    demuxSheetLines.append(",,,")

    # replace nans with empty string
    ssdf_towrite = ssDic['sampleSheet'].fillna('')

    # Todo - deduplicate this mess.
    if ssDic['dualIx']:
        demuxSheetLines.append("[BCLConvert_Data],,,,,,")
        if laneSplitStatus:
            demuxSheetLines.append(
                "Lane,Sample_ID,index,index2,Sample_Project"
            )
            for index, row in ssdf_towrite.iterrows():
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
            for index, row in ssdf_towrite.iterrows():
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
            for index, row in ssdf_towrite.iterrows():
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
            for index, row in ssdf_towrite.iterrows():
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


def writeDemuxSheetAviti(demuxOut, ssDic, laneSplitStatus):
    demuxSheetLines = []
    # Header first.
    demuxSheetLines.append("[SETTINGS],,,,")
    demuxSheetLines.append("SettingName,Value,Lane,,")
    if 'mismatch' in ssDic:
        for i in (1, 2):
            bc_str = 'I{}MismatchThreshold'.format(i)
            if i == 2 and not ssDic['dualIx']:
                continue
            if bc_str in ssDic['mismatch']:
                demuxSheetLines.append(
                    "{},{},,,".format(bc_str, ssDic['mismatch'][bc_str])
                )

            # TODO do we want this behavior?
            # log.warning("dualIx set, but no mismatch returned. Overriding.")
            # ssDic['dualIx'] = False
    
    for k,v in ssDic['mask'].items():

        demuxSheetLines.append("{},{},,,".format(k,v))
    
    if len(ssDic['convertOpts']) > 0:
        for opts in ssDic['convertOpts']:
            demuxSheetLines.append(opts)

    # replace nans with empty string
    ssdf_towrite = ssDic['sampleSheet'].fillna('')

    # Todo - deduplicate this mess.
    if ssDic['dualIx']:
        demuxSheetLines.append("[SAMPLES],,,,")
        demuxSheetLines.append(
            "SampleName,Index1,Index2,Lane,Project"
        )
        for index, row in ssdf_towrite.iterrows():
            demuxSheetLines.append(
                joinLis(
                    [
                        row['SampleName'],
                        row['Index1'],
                        row['Index2'],
                        row['Lane'],
                        row['Project']
                    ],
                    joinStr=","
                )
            )

    else:
        demuxSheetLines.append("[SAMPLES],,,,")
        demuxSheetLines.append(
            "SampleName,Index1,Lane,Project"
        )
        for index, row in ssdf_towrite.iterrows():
            demuxSheetLines.append(
                joinLis(
                    [
                        row['SampleName'],
                        row['Index1'],
                        row['Lane'],
                        row['Project']
                    ],
                    joinStr=","
                )
            )
    
    with open(demuxOut, 'w') as f:
        for line in demuxSheetLines:
            f.write('{}\n'.format(line))


def readDemuxSheet(demuxSheet, what='all'):
    '''
    In case of manual intervention.
    We want to have the correct info in reports / emails.
    Thus, we reparse the demuxSheet just to be sure.
    We check for:
     - 'mask' (overridecycles)
     - indices used.
     - mismatches definition
    '''
    with open(demuxSheet) as f:
        sampleStatus = False
        nesLis = []
        mmdic = {}
        for line in f:
            line = line.strip()
            if line.startswith('BarcodeMismatchesIndex1'):
                mmdic['BarcodeMismatchesIndex1'] = int(line.split(',')[1])
            if line.startswith('BarcodeMismatchesIndex2'):
                mmdic['BarcodeMismatchesIndex2'] = int(line.split(',')[1])
            if line.startswith('OverrideCycles'):
                mask = line.replace(
                    'OverrideCycles', ''
                ).replace(
                    ',', ''
                )
            if sampleStatus:
                nesLis.append(line.split(','))

            if (
                line.startswith('[BCLConvert_Data]')
                or line.startswith('[Data]')
            ):
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
    if what == 'all':
        return (mask, df, dualIx, mmdic)
    elif what == 'df':
        return (df)


def parseStats(outputFolder, ssdf, mode='illumina') -> pd.DataFrame:
    if mode == 'illumina':
        QCmetFile = outputFolder / 'Reports' / 'Quality_Metrics.csv'
        DemuxmetFile = outputFolder / 'Reports' / 'Demultiplex_Stats.csv'
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
                'gotDepth': int(muxDic[ID])
            }
        MetrixDF = pd.DataFrame(MetrixDic).T
        MetrixDF['Sample_ID'] = MetrixDF.index
        # left join to only get samples already present
        newDF = pd.merge(ssdf, MetrixDF, on='Sample_ID', how='left')
        newDF = newDF[newDF['Sample_ID'] != 'Undetermined']
        return newDF
    elif mode =='aviti':
        # After bases2fq, samples are organized under Samples/project/sampleID
        # In that directory is a stats.json file we can use to get statistics.
        for samplestat in (outputFolder).rglob("*_stats.json"):
            with open(samplestat) as f:
                stats = json.load(f)
            sampleID = samplestat.name.replace('_stats.json', '')
            
            if sampleID not in ssdf['Sample_ID'].values:
                logging.warning(f"Sample {sampleID} not found in sampleSheet.")
                continue
            _qsm = stats['QualityScoreMean']
            _q30 = stats['PercentQ30']
            _depth = stats['NumPolonies']
            match len(stats['Reads']):
                case 1:
                    # Single-end
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'meanQ'] = f"1:{_qsm}"
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'percQ30'] = f"1:{_q30}"
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'gotDepth'] = int(_depth)
                case 2:
                    # paired-end
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'meanQ'] = f"1:{_qsm},2:{_qsm}"
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'percQ30'] = f"1:{_q30},2:{_q30}"
                    ssdf.loc[ssdf['Sample_ID'] == sampleID, 'gotDepth'] = int(_depth)
        return ssdf
    else:
        logging.error(f"parseStats - mode not supported: {mode}")
        

def compareDemuxSheet(ssDic, demuxSheet):
    '''
    a demuxSheet already exists, but we inferred demux parameters from the original data.
    This happens in case we want to override some settings.
    This function updates these parameters if necessary
    Parameters:
      - ssDic: Dictionary with parameters from a specific outLane, comes from sampleSheet class.
      - demuxSheet: Path to the existing demuxSheet.csv
    '''
    man_mask, man_df, man_dualIx, man_mmdic = readDemuxSheet(
        demuxSheet
    )
    # mmdic
    if ssDic['mismatch'] != man_mmdic:
        logging.warning(f"Demux - mmdic set as {man_mmdic}")
        ssDic['mismatch'] = man_mmdic
    # mask
    if 'mask' in ssDic and man_mask != ssDic['mask']:
        logging.warning(f"Demux - mask set as {man_mask}")
        ssDic['mask'] = man_mask
    # dualIx
    if 'dualIx' in ssDic and man_dualIx != ssDic['dualIx']:
        logging.warning(f"Demux - dualIx set as {man_dualIx}")
        ssDic['dualIx'] = man_dualIx
    # sampleSheet
    ssDic['sampleSheet'] = matchingSheets(
        ssDic['sampleSheet'],
        man_df
    )
    # P5RC
    ssDic['P5RC'] = False
    if demuxSheet.with_suffix('.bak').exists():
        logging.warning("Demux - P5RC detected")
        ssDic['P5RC'] = True


def matchingSheets(autodf, mandf):
    '''
    if demuxSheet is overwritten:
        update indices from the manually-edited dataframe to the
        given automatically generated DataFrame.

          :param autodf: Automatically generated sampleSheet DataFrame
          from SampleSheet.ssDic[lane]['samplesheet'].
          :param mandf: Manually edited samplesheet (from demuxSheet.csv in
          output lane directory)
          :type autodf: pandas.DataFrame
          :type mandf: pandas.DataFrame
          :returns: updated autodf with indices
    '''
    # if there are no indices, just return the same autodf
    if 'index' not in autodf.columns:
        return autodf

    if len(autodf.index) != len(mandf.index):
        logging.warning(
            "Demux - number of samples changed in overwritten demuxSheet !"
        )

    dualIx = 'index2' in list(mandf.columns)

    for index, row in mandf.iterrows():
        sample_ID = row['Sample_ID']
        index = row['index']
        if dualIx:
            index2 = row['index2']
        # grab the index in the autodf.
        pdIx = autodf[autodf['Sample_ID'] == sample_ID].index
        if dualIx:
            if autodf.loc[pdIx, 'index'].values != index:
                oriP7 = autodf.loc[pdIx, 'index'].values[0]
                logging.debug(f"Demux - Changing P7 {oriP7} to {index} for {sample_ID}")
                autodf.loc[pdIx, 'index'] = index
                autodf.loc[pdIx, 'I7_Index_ID'] = np.nan
            if autodf.loc[pdIx, 'index2'].values != index2:
                oriP5 = autodf.loc[pdIx, 'index2'].values[0]
                logging.debug(f"Demux - Changing P5 {oriP5} to {index2} for {sample_ID}")
                autodf.loc[pdIx, 'index2'] = index2
                autodf.loc[pdIx, 'I5_Index_ID'] = np.nan
        else:
            # check index1, set index2 to na
            if autodf.loc[pdIx, 'index'].values != index:
                oriP7 = autodf.loc[pdIx, 'index'].values[0]
                logging.debug(f"Demux - Changing P7 {oriP7} to {index} for {sample_ID}")
                autodf.loc[pdIx, 'index'] = index
                # change type as well!
                autodf.loc[pdIx, 'I7_Index_ID'] = np.nan
                # it's not dualIx, so set index2/I5_Index_ID to nan.
                if 'index2' in list(autodf.columns):
                    autodf.loc[pdIx, 'index2'] = np.nan
                if 'I5_Index_ID' in list(autodf.columns):
                    autodf.loc[pdIx, 'I5_Index_ID'] = np.nan
    return (autodf)


def evalMiSeqP5(outPath, dualIx):
    '''
    Evaluates MiSeq runs if P5's need to be RC'ed.
    Takes the path for an outlane,
    find out if all samples are empty
    if that is the case, and the run is dualindexed,
    rerun bclConvert with all P5s RC'ed.
    '''
    # Known barcodes
    kbcDF = pd.read_csv( outPath / 'Reports' / 'Demultiplex_Stats.csv' )
    # Test if > 90% of samples are virtually empty.
    numLowreadSamples = len(kbcDF[kbcDF['# Reads'] < 1000])
    totalSamples = len(kbcDF[kbcDF['SampleID'] != 'Undetermined'])
    if not numLowreadSamples/totalSamples == 1:
        return False
    logging.warning(
        'Demux - EvalP5: More then 90% samples empty. Attempting to salvage by RC the P5.'
    )
    if not dualIx:  # Only RC P5 operations for now.
        return False

    # Read demuxSheet
    demuxSheetPath = Path(outPath, 'demuxSheet.csv')
    demuxSheet = []
    with open(demuxSheetPath) as f:
        headStatus = True
        for line in f:
            if 'Sample_ID' in line.strip():
                headStatus = False
                colnames = line.strip().split(',')
                demuxSheet.append(colnames)
            if headStatus:
                demuxSheet.append(line.strip().split(','))
            else:
                if 'Sample_ID' not in line.strip():
                    demuxSheetLine = line.strip().split(',')
                    ixPos = colnames.index('index2')
                    oldIx = demuxSheetLine[ixPos]
                    newIx = str(Seq(oldIx).reverse_complement())
                    demuxSheetLine[ixPos] = newIx
                    demuxSheet.append(demuxSheetLine)
    shutil.move(
        demuxSheetPath,
        demuxSheetPath.with_suffix('.bak')
    )
    with open(demuxSheetPath, 'w') as f:
        for _l in demuxSheet:
            f.write(','.join(_l) + '\n')
    return True
