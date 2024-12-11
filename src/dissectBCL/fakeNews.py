from dissectBCL.misc import fetchLatestSeqDir
from dissectBCL.misc import fexUpload
from dissectBCL.misc import getDiskSpace
from dissectBCL.misc import joinLis
from dissectBCL.misc import matchOptdupsReqs
from dissectBCL.misc import sendMqcReports
from dissectBCL.misc import stripRights
from dissectBCL.misc import umlautDestroyer
import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from importlib.metadata import version
import interop
import json
import logging
import os
import pandas as pd
from pathlib import Path
import requests
import ruamel.yaml
import shutil
import smtplib
import sys



def pullParkour(flowcellID, config):
    """
    Look for the flowcell/lane in parkour for the library type.
    The flowcell ID is of form:
     - 210608_A00931_0309_BHCCMWDRXY
     - 211105_M01358_0001_000000000-JTYPH
    we need "HCCMWDRXY" or "JTYPH" for the request.
    """
    FID = flowcellID.split('_')[3][1::]
    if '-' in FID:
        FID = FID.split('-')[1]
    logging.info(
        "Pulling parkour for with flowcell {} using FID {}".format(
            flowcellID,
            FID
        )
    )
    d = {'flowcell_id': FID}
    res = requests.get(
        config['parkour']['URL'] + '/api/analysis_list/analysis_list/',
        auth=(
            config['parkour']['user'],
            config['parkour']['password']
        ),
        params=d,
        verify=config['parkour']['cert']
    )
    if res.status_code == 200:
        logging.info("parkour API code 200")
        """
        res.json() returns nested dict
        {
            project:
                sampleID:
                [
                    name,
                    libType,
                    protocol,
                    genome,
                    indexType,
                    depth
                ]
        }
        We flatten it and return.
        """
        flatLis = []
        for project in res.json():
            for sample in res.json()[project]:
                flatLis.append(
                    [
                        umlautDestroyer(project), sample
                    ] + res.json()[project][sample]
                )
        parkourDF = pd.DataFrame(flatLis)
        parkourDF.columns = [
                'Sample_Project',
                'Sample_ID',
                'Sample_Name',
                'Library_Type',
                'Description',
                'Organism',
                'indexType',
                'reqDepth'
            ]
        # Sanitize sample names.
        parkourDF['Sample_Name'] = parkourDF['Sample_Name'].apply(
            lambda x: umlautDestroyer(x)
        )
        # parkour lists requested in millions.
        parkourDF['reqDepth'] = parkourDF['reqDepth']*1000000
        # Some exceptions where there is a ' in the description..
        parkourDF['Description'] = parkourDF[
            'Description'
        ].str.replace(r"[\’,]", '', regex=True)
        return parkourDF
    logging.warning("parkour API not 200!")
    mailHome(
        flowcellID,
        "Parkour pull failed [{}]".format(
            res.status_code
        ),
        config
    )
    sys.exit("Parkour pull failed.")


def pushParkour(flowcellID, sampleSheet, config, flowcellBase):
    # pushing out the 'Run statistics in parkour'.
    '''
    we need:
     - R1 > Q30 %
     - R2 > Q30 % (if available)
     - clusterPF (%)
     - name (== laneStr)
     - undetermined_indices (%)
     - reads_pf (M)
    '''

    # Parse interop.
    iop_df = pd.DataFrame(
        interop.summary(
            interop.read(
                str(flowcellBase)
            ),
            'Lane'
        )
    )

    FID = flowcellID
    if '-' in FID:
        FID = FID.split('-')[1]
    d = {}
    d['flowcell_id'] = FID
    laneDict = {}
    for outLane in sampleSheet.ssDic:
        # Quality_Metrics.csv contains all the info we need.
        qMetPath = Path(config['Dirs']['outputDir'], outLane, 'Reports', 'Quality_Metrics.csv')
        qdf = pd.read_csv(qMetPath)
        # If a flowcell is split, qMetPath contains only Lane 1 e.g.
        # If not, all lanes sit in qdf
        # Anyhow, iterating and filling will capture all we need.
        for lane in list(qdf['Lane'].unique()):
            subdf = qdf[qdf['Lane'] == lane]
            laneStr = 'Lane {}'.format(lane)
            laneDict[laneStr] = {}
            # reads PF.
            readsPF = iop_df[
                (iop_df['ReadNumber'] == 1) & (iop_df['Lane'] == lane)
            ]['Reads Pf'].values
            logging.info('lane {}, reads PF = {}'.format(lane, float(readsPF)))
            laneDict[laneStr]['reads_pf'] = float(readsPF)
            # Und indices.
            laneDict[laneStr]["undetermined_indices"] = \
                round(
                    subdf[
                        subdf["SampleID"] == "Undetermined"
                    ]["YieldQ30"].sum() / subdf['YieldQ30'].sum() * 100,
                    2
                )
            Q30Dic = subdf[subdf['SampleID'] != 'Undetermined'].groupby("ReadNumber")['% Q30'].mean().to_dict()
            for read in Q30Dic:
                if 'I' not in str(read):
                    readStr = 'read_{}'.format(read)
                    laneDict[laneStr][readStr] = round(Q30Dic[read]*100, 2)
            laneDict[laneStr]["cluster_pf"] = round(
                subdf[subdf['SampleID'] != 'Undetermined']["YieldQ30"].sum()/subdf[subdf['SampleID'] != 'Undetermined']["Yield"].sum() * 100,
                2
            )
            laneDict[laneStr]["name"] = laneStr
    d['matrix'] = json.dumps(list(laneDict.values()))
    logging.info(f"fakenews - pushParkour - Pushing FID with dic {FID} {d}")
    pushParkStat = requests.post(
        config.get("parkour", "URL") + '/api/run_statistics/upload/',
        auth=(
            config.get("parkour", "user"),
            config.get("parkour", "password")
        ),
        data=d,
        verify=config['parkour']['cert']
    )
    logging.info("fakenews - ParkourPush - return {}".format(pushParkStat))
    return pushParkStat


def mailHome(subject, _html, config, toCore=False):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = f"[{config['communication']['subject']}] [{version('dissectBCL')}] " + subject
    mailer['From'] = config['communication']['fromAddress']
    if toCore:
        mailer['To'] = config['communication']['bioinfoCore']
    else:
        mailer['To'] = config['communication']['finishedTo']
    email = MIMEText(_html, 'html')
    mailer.attach(email)
    s = smtplib.SMTP(config['communication']['host'])
    if toCore:
        s.sendmail(
            config['communication']['fromAddress'],
            config['communication']['bioinfoCore'],
            mailer.as_string()
            )
    else:
        s.sendmail(
            config['communication']['fromAddress'],
            config['communication']['finishedTo'].split(', '),
            mailer.as_string()
            )
    s.quit()


def shipFiles(outPath, config):
    transferStart = datetime.datetime.now()
    shipDic = {}
    outLane = outPath.name
    # Get directories from outPath.
    for projectPath in outPath.glob('Project*'):
        project = projectPath.name
        shipDic[project] = 'No'
        logging.info("fakenews - Shipping {}".format(project))
        PI = project.split('_')[-1].lower().replace(
            "cabezas-wallscheid", "cabezas"
        )
        fqcPath = Path(str(projectPath).replace("Project_", "FASTQC_Project_"))
        if PI in config['Internals']['PIs']:
            # Shipping
            fqc = fqcPath.name
            enduserBase = fetchLatestSeqDir(config, PI) / outLane
            logging.info(f"fakenews - Found {PI}. Shipping internally to {enduserBase}.")
            enduserBase.mkdir(mode=0o750, exist_ok=True)
            replaceStatus = 'Copied'
            if (enduserBase / fqc).exists():
                shutil.rmtree(enduserBase / fqc)
                replaceStatus = 'Replaced'
            try:
                shutil.copytree(fqcPath, enduserBase / fqc)
            except OSError:
                logging.critical(f"Copying {fqcPath} into {enduserBase} failed.")
                mailHome(
                    outPath.name,
                    f"{fqcPath} copying into {enduserBase} failed.",
                    config
                )
                sys.exit()
            if (enduserBase / project).exists():
                shutil.rmtree(enduserBase / project)
                replaceStatus = 'Replaced'
            try:
                shutil.copytree(projectPath, enduserBase / project)
            except OSError:
                logging.critical(f"Copying {projectPath} into {enduserBase} failed.")
                mailHome(
                    outPath.name,
                    f"{projectPath} copying into {enduserBase} failed.",
                    config
                )
                sys.exit()
            # Strip rights
            stripRights(enduserBase)
            shipDic[project] = [replaceStatus, f"{getDiskSpace(enduserBase)[1]}GB free"]
        else:
            if not config['Internals'].getboolean('fex'):
                shipDic[project] = "Ignored( by config)"
                logging.info(f"fakenews - {project} not fex uploaded by config.")
            else:
                shipDic[project] = fexUpload(
                    outLane, project, config['communication']['fromAddress'],
                    (projectPath, fqcPath)
                )
    sendMqcReports(outPath, config['Dirs'])
    transferStop = datetime.datetime.now()
    transferTime = transferStop - transferStart
    return {'transfertime': transferTime, 'shipDic': shipDic}


def organiseLogs(flowcell, sampleSheet):
    for outLane in sampleSheet.ssDic:
        logging.info("Populating log dir for {}".format(outLane))
        _logDir = os.path.join(
            flowcell.outBaseDir,
            outLane,
            'Logs'
        )
        _logBCLDir = os.path.join(
            _logDir,
            'BCLConvert'
        )
        # move bclConvert logFiles.
        if not os.path.exists(_logBCLDir):
            os.mkdir(_logBCLDir)
            bclConvertFiles = [
                'Errors.log',
                'FastqComplete.txt',
                'Info.log',
                'Warnings.log'
            ]
            for mvFile in bclConvertFiles:
                fileIn = os.path.join(
                    _logDir,
                    mvFile
                )
                fileOut = os.path.join(
                    _logBCLDir,
                    mvFile
                )
                shutil.move(fileIn, fileOut)

        # Write out ssdf.
        outssdf = os.path.join(_logDir, 'sampleSheetdf.tsv')
        sampleSheet.ssDic[outLane]['sampleSheet'].to_csv(outssdf, sep='\t')

        # write out outLaneInfo.yaml
        dic0 = sampleSheet.ssDic[outLane]
        del dic0['sampleSheet']
        yaml0 = ruamel.yaml.YAML()
        yaml0.indent(mapping=2, sequence=4, offset=2)
        outLaneInfo = os.path.join(_logDir, 'outLaneInfo.yaml')
        with open(outLaneInfo, 'w') as f:
            yaml0.dump(dic0, f)

        # write out config.ini
        dic1 = flowcell.asdict()
        flowcellConfig = os.path.join(_logDir, 'config.ini')
        with open(flowcellConfig, 'w') as f:
            dic1['config'].write(f)

        # write out flowcellInfo.yaml
        del dic1['config']
        yaml1 = ruamel.yaml.YAML()
        yaml1.indent(mapping=2, sequence=4, offset=2)
        flowcellInfo = os.path.join(_logDir, 'flowcellInfo.yaml')
        with open(flowcellInfo, 'w') as f:
            yaml1.dump(dic1, f)

# outPath, initTime, flowcellID, ssDic, transferTime, exitStats, solPath
def gatherFinalMetrics(outLane, flowcell):
    logging.info(f"fakenews - gatherFinalMetrics - {outLane}")
    outPath = flowcell.outBaseDir / outLane
    ssDic = flowcell.sampleSheet.ssDic[outLane] 
    ssdf = ssDic['sampleSheet']
    barcodeMask = ssDic['mask']
    mismatch = " ".join(
        [i + ': ' + str(j) for i, j in ssDic['mismatch'].items()]
    )
    # Get undetermined
    muxDF = pd.read_csv(outPath / 'Reports' / 'Demultiplex_Stats.csv')
    totalReads = int(muxDF['# Reads'].sum())
    if len(muxDF[muxDF['SampleID'] == 'Undetermined']) == 1:
        undReads = int(
            muxDF[
                muxDF['SampleID'] == 'Undetermined'
            ]['# Reads'].iloc[0]
        )
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
    bcDF = pd.read_csv(outPath / 'Reports' / 'Top_Unknown_Barcodes.csv')
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
    runTime = datetime.datetime.now() - flowcell.startTime
    # optDups
    optDups = []
    for opt in outPath.glob("*/*/*duplicate.txt"):
        project = opt.parts[-3].replace("FASTQC_", "")
        sample = opt.name.replace(".duplicate.txt", "")
        sampleID = opt.parts[-2].replace("Sample_", "")
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
    # Fetch organism and kraken reports
    sampleDiv = {}
    for screen in outPath.glob("*/*/*.rep"):
        sampleID = screen.parts[-2].replace("Sample_", "")
        sample = screen.name.replace('.rep', '')

        # samples with 0 reads still make an empty report.
        # hence the try / except.
        # 'mouse (GRCm39)' -> 'mouse'
        parkourOrg = str(
            ssdf[ssdf["Sample_ID"] == sampleID]['Organism'].values[0][1]
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

    return {
        'undetermined':undReads,
        'totalReads':totalReads,
        'topBarcodes':BCDic,
        'spaceFree_rap':getDiskSpace(outPath),
        'spaceFree_sol':getDiskSpace(flowcell.bclPath),
        'runTime':runTime,
        'optDup':optDups,
        'flowcellID':flowcell.flowcellID,
        'outLane':outLane,
        'contamination':sampleDiv,
        'mismatch':mismatch,
        'barcodeMask':barcodeMask,
        'transferTime': flowcell.transferTime,
        'exitStats': flowcell.exitStats,
        'P5RC':ssDic['P5RC']
    }
