from dissectBCL.classes import drHouseClass
from dissectBCL.logger import log
from dissectBCL.misc import joinLis
import pandas as pd
import os
import shutil
import glob
import datetime
import re


def getDiskSpace(outputDir):
    total, used, free = shutil.disk_usage(outputDir)
    return(total // (2**30), free // (2**30))


def matchOptdupsReqs(optDups, ssdf):
    _optDups = []
    for lis in optDups:
        sample = lis[1]
        req = ssdf[
            ssdf['Sample_Name'] == sample
        ]['reqDepth'].values
        got = ssdf[
            ssdf['Sample_Name'] == sample
        ]['gotDepth'].values
        reqvgot = float(got/req)
        _optDups.append(
            [lis[0], sample, lis[2], round(reqvgot, 2)]
        )
    return(_optDups)

def matchIDtoName(optDups, ssdf):
    _optDups = []
    for lis in optDups:
        sample = lis[1]
        sid = ssdf[
            ssdf['Sample_Name'] == sample
        ]['Sample_ID'].values[0]
        _optDups.append(
            [lis[0], sample, lis[2], lis[3], sid]
        )
    return(sorted(_optDups, key=lambda x: x[4]))


def initClass(
    outPath, initTime, flowcellID, ssDic, transferTime, exitStats, solPath
        ):
    log.info("init drHouse for {}".format(outPath))
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
    print(bcDF)
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
        BCDic[entry[0]] = [round(entry[1]/1000000, 0), entry[2]]
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
        with open(opt) as f:
            dups = f.read()
            dups = dups.strip().split()
            if float(dups[1]) != 0:
                optDups.append(
                    [
                        project,
                        sample,
                        round(100*float(dups[0])/float(dups[1]), 2)
                    ]
                )
            else:
                optDups.append(
                    [
                        project,
                        sample,
                        "NA"
                    ]
                )
    projSamDic = pd.Series(
        ssdf['Sample_Project'].values,
        index=ssdf['Sample_Name']
    ).to_dict()
    for sample in projSamDic:
        if not any(sample in sl for sl in optDups):
            optDups.append(
                [
                    projSamDic[sample],
                    sample,
                    'NA'
                ]
            )
    optDups = matchOptdupsReqs(optDups, ssdf)
    optDups = matchIDtoName(optDups, ssdf)
    # Fetch organism and fastqScreen
    sampleDiv = {}
    for screen in glob.glob(
        os.path.join(
            outPath,
            '*/*/*screen.txt'
        )
    ):
        screenDF = pd.read_csv(
            screen, sep='\t', skip_blank_lines=True, header=0, skiprows=[0]
        )
        screenDF = screenDF.dropna()
        sample = re.sub('_R[123]_screen.txt', '', screen.split('/')[-1])
        # ParkourOrganism
        parkourOrg = str(  # To string since NA is a float
            ssdf[ssdf["Sample_Name"] == sample]['Organism'].values[0]
        )
        # Top_oneonone
        if not screenDF['#One_hit_one_genome'].sum() == 0:
            maxHit = screenDF["%One_hit_one_genome"].max()
            fqscreenOrg = screenDF[
                screenDF['%One_hit_one_genome'] == maxHit
            ]['Genome'].values[0]
            sampleDiv[sample] = [maxHit, fqscreenOrg, parkourOrg]
        else:
            sampleDiv[sample] = ['NA', 'NA', parkourOrg]
    return(drHouseClass(
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
