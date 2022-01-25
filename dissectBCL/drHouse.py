from dissectBCL.classes import drHouseClass
import pandas as pd
import os
import shutil
import glob
import datetime
import numpy as np
import re

def getDiskSpace(outputDir):
    total, used, free = shutil.disk_usage(outputDir)
    return(total // (2**30), free // (2**30))


def initClass(outPath, initTime, flowcellID, ssdf):
    # Get undetermined
    muxPath = os.path.join(
        outPath,
        'Reports',
        'Demultiplex_Stats.csv'
    )
    muxDF = pd.read_csv(muxPath)
    totalReads = int(muxDF['# Reads'].sum())
    undReads = int(muxDF[muxDF['SampleID'] == 'Undetermined']['# Reads'])
    # topBarcodes
    BCPath = os.path.join(
        outPath,
        'Reports',
        'Top_Unknown_Barcodes.csv'
    )
    bcDF = pd.read_csv(BCPath)
    bcDF = bcDF.head(5)
    BCs = [
        '+'.join(list(x)) for x in bcDF.filter(like='index', axis=1).values
    ]
    BCReads = list(bcDF['# Reads'])
    BCDic = dict(zip(
        BCs, BCReads
    ))
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
        print(opt)
        project = opt.split('/')[-3].replace("FASTQC_","")
        sample = opt.split('/')[-1].replace(".duplicate.txt", "")
        with open(opt) as f:
            dups = f.read()
            dups = dups.strip().split()
            optDups.append(
                [
                    project,
                    sample,
                    round(100*float(dups[0])/float(dups[1]), 2)
                ]
            )
    projSamDic = pd.Series(
        ssdf['Sample_Project'].values,
        index = ssdf['Sample_Name']
    ).to_dict()
    for sample in projSamDic:
        if not any(sample in sl for sl in optDups):
            optDups.append(
                [projSamDic[sample],
                sample,
                'NA']
            )
    # Fetch organism and fastqScreen
    sampleDiv = {}
    if 'Organism' in list(ssdf.columns):
        for screen in glob.glob(
            os.path.join(
                outPath,
                '*/*/*screen.txt'
            )
        ):
            screenDF = pd.read_csv(screen, sep='\t', skip_blank_lines=True, header=0, skiprows=[0])
            screenDF = screenDF.dropna()
            sample = re.sub('_R[123]_screen.txt', '', screen.split('/')[-1])
            # Simpson diversity.
            simpson = sum(
                [(i/screenDF['#One_hit_one_genome'].sum())**2 for i in screenDF['#One_hit_one_genome']]
            )
            sampleDiv[sample] = round(simpson,2)
    return(drHouseClass(
        undetermined = undReads,
        totalReads = totalReads,
        topBarcodes = BCDic,
        spaceFree = getDiskSpace(outPath),
        runTime = runTime,
        optDup = optDups,
        flowcellID = flowcellID,
        outLane = outPath.split('/')[-1],
        simpson = sampleDiv
    ))