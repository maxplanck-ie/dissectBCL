from dissectBCL.classes import drHouseClass
from dissectBCL.misc import joinLis
import pandas as pd
import os
import shutil
import glob
import datetime
import logging


def getDiskSpace(outputDir):
    '''
    Return space free in GB
    '''
    total, used, free = shutil.disk_usage(outputDir)
    return (total // (2**30), free // (2**30))


def matchOptdupsReqs(optDups, ssdf):
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
        reqvgot = float(got/req)
        _optDups.append(
            [lis[0], sampleID, sampleName, lis[3], round(reqvgot, 2), int(got)]
        )
    return (sorted(_optDups, key=lambda x: x[1]))


def initClass(
    outPath, initTime, flowcellID, ssDic, transferTime, exitStats, solPath
        ):
    logging.info("init drHouse for {}".format(outPath))
    ssdf = ssDic['sampleSheet']
    barcodeMask = ssDic['mask']
    mismatch = " ".join(
        [i + ': ' + str(j) for i, j in ssDic['mismatch'].items()]
    )
    # Get undetermined
    muxPath = os.path.join(
        outPath,
        'Reports',
        'Demultiplex_Stats.csv'
    )
    muxDF = pd.read_csv(muxPath)
    totalReads = int(muxDF['# Reads'].sum())
    if len(muxDF[muxDF['SampleID'] == 'Undetermined']) == 1:
        undReads = int(muxDF[muxDF['SampleID'] == 'Undetermined']['# Reads'])
    else:
        undDic = dict(
            muxDF[
                muxDF['SampleID'] == 'Undetermined'
            ][['Lane', '# Reads']].values
        )
        undStr = ""
        for lane in undDic:
            undStr += "Lane {}: {}% {}M, ".format(
                lane,
                round(100*undDic[lane]/totalReads, 2),
                round(undDic[lane]/1000000, 2)
            )
            undReads = undStr[:-2]
    # topBarcodes
    BCPath = os.path.join(
        outPath,
        'Reports',
        'Top_Unknown_Barcodes.csv'
    )
    bcDF = pd.read_csv(BCPath)
    bcDF = bcDF.head(5)
    BCs = [
        joinLis(
            list(x), joinStr='+'
        ) for x in bcDF.filter(like='index', axis=1).values
    ]
    BCReads = list(bcDF['# Reads'])
    BCReadsPerc = list(bcDF['% of Unknown Barcodes'])
    BCDic = {}
    for entry in list(
        zip(BCs, BCReads, BCReadsPerc)
    ):
        BCDic[entry[0]] = [round(float(entry[1])/1000000, 2), entry[2]]
    # runTime
    runTime = datetime.datetime.now() - initTime
    # optDups
    optDups = []
    for opt in glob.glob(
        os.path.join(
            outPath,
            '*/*/*duplicate.txt'
        )
    ):
        project = opt.split('/')[-3].replace("FASTQC_", "")
        sample = opt.split('/')[-1].replace(".duplicate.txt", "")
        sampleID = opt.split('/')[-2].replace("Sample_", "")
        with open(opt) as f:
            dups = f.read()
            dups = dups.strip().split()
            if float(dups[1]) != 0:
                optDups.append(
                    [
                        project,
                        sampleID,
                        sample,
                        round(100*float(dups[0])/float(dups[1]), 2)
                    ]
                )
            else:
                optDups.append(
                    [
                        project,
                        sampleID,
                        sample,
                        "NA"
                    ]
                )
    IDprojectDic = pd.Series(
        ssdf['Sample_Project'].values,
        index=ssdf['Sample_ID']
    ).to_dict()
    nameIDDic = pd.Series(
        ssdf['Sample_Name'].values,
        index=ssdf['Sample_ID']
    ).to_dict()
    for sampleID in nameIDDic:
        if not any(sampleID in sl for sl in optDups):
            optDups.append(
                [
                    IDprojectDic[sampleID],
                    sampleID,
                    nameIDDic[sampleID],
                    'NA'
                ]
            )
    optDups = matchOptdupsReqs(optDups, ssdf)
    # optDups = matchIDtoName(optDups, ssdf)
    # Fetch organism and fastqScreen
    sampleDiv = {}
    for screen in glob.glob(
        os.path.join(
            outPath,
            '*/*/*.rep'
        )
    ):
        sampleID = screen.split('/')[-2].replace("Sample_", "")
        sample = screen.split('/')[-1].replace('.rep', '')

        # samples with 0 reads still make an empty report.
        # hence the try / except.
        parkourOrg = str(  # To string since NA is a float
                ssdf[ssdf["Sample_ID"] == sampleID]['Organism'].values[0]
            )
        try:
            screenDF = pd.read_csv(
                screen, sep='\t', header=None
            )
            # tophit == max in column 2.
            # ParkourOrganism
            krakenOrg = screenDF.iloc[
                screenDF[2].idxmax()
            ][5].replace(' ', '')
            fraction = round(
                screenDF[2].max()/screenDF[2].sum(),
                2
            )
            sampleDiv[sampleID] = [fraction, krakenOrg, parkourOrg]
        except pd.errors.EmptyDataError:
            sampleDiv[sampleID] = ['NA', 'None', parkourOrg]

    return (drHouseClass(
        undetermined=undReads,
        totalReads=totalReads,
        topBarcodes=BCDic,
        spaceFree_rap=getDiskSpace(outPath),
        spaceFree_sol=getDiskSpace(solPath),
        runTime=runTime,
        optDup=optDups,
        flowcellID=flowcellID,
        outLane=outPath.split('/')[-1],
        contamination=sampleDiv,
        mismatch=mismatch,
        barcodeMask=barcodeMask,
        transferTime=transferTime,
        exitStats=exitStats
    ))
