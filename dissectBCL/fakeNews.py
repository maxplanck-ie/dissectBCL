import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import requests
import pandas as pd
from dissectBCL.logger import log
from dissectBCL.misc import retBCstr, retIxtype, retMean_perc_Q
from dissectBCL.misc import fetchLatestSeqDir, formatSeqRecipe
from dissectBCL.misc import umlautDestroyer, formatMisMatches
import os
import shutil
import smtplib
import glob
import ruamel.yaml
import json
from subprocess import check_output, Popen
import sys
import numpy as np


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
    log.info(
        "Pulling parkour for with flowcell {} using FID {}".format(
            flowcellID,
            FID
        )
    )
    d = {'flowcell_id': FID}
    res = requests.get(
        config['parkour']['pullURL'],
        auth=(
            config['parkour']['user'],
            config['parkour']['password']
        ),
        params=d
    )
    if res.status_code == 200:
        log.info("parkour API code 200")
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
        return parkourDF
    log.warning("parkour API not 200!")
    sys.exit("Parkour pull failed.")


def pushParkour(flowcellID, sampleSheet, config):
    # pushing out the 'Run statistics in parkour'.
    '''
    we need:
     - R1 > Q30 % - done
     - R2 > Q30 % (if available) - done
     - clusterPF (%) - done
     - name (== laneStr) - done
     - undetermined_indices (%) - done
     - reads_pf (M) - no longer available (only Bases)

    '''
    d = {}
    d['flowcell_id'] = flowcellID
    laneDict = {}
    for outLane in sampleSheet.ssDic:
        # Quality_Metrics.csv contains all the info we need.
        qMetPath = os.path.join(
            config['Dirs']['outputDir'],
            outLane,
            'Reports',
            'Quality_Metrics.csv'
        )
        qdf = pd.read_csv(qMetPath)
        # If a flowcell is split, qMetPath contains only Lane 1 e.g.
        # If not, all lanes sit in qdf
        # Anyhow, iterating and filling will capture all we need.
        for lane in list(qdf['Lane'].unique()):
            subdf = qdf[qdf['Lane'] == lane]
            laneStr = 'Lane {}'.format(lane)
            laneDict[laneStr] = {}
            # Und indices.
            laneDict[laneStr]["undetermined_indices"] = \
                round(
                    subdf[
                        subdf["SampleID"] == "Undetermined"
                    ]["YieldQ30"].sum() / subdf['YieldQ30'].sum() * 100,
                    2
                )
            Q30Dic = subdf.groupby("ReadNumber").mean()['% Q30'].to_dict()
            for read in Q30Dic:
                readStr = 'read_{}'.format(read)
                laneDict[laneStr][readStr] = round(Q30Dic[read]*100, 2)
            laneDict[laneStr]["cluster_pf"] = round(
                subdf["YieldQ30"].sum()/subdf["Yield"].sum() * 100,
                2
            )
            laneDict[laneStr]["name"] = laneStr
    d['matrix'] = json.dumps(list(laneDict.values()))
    log.info("Pushing FID with dic {} {}".format(flowcellID, d))
    pushParkStat = requests.post(
        config.get("parkour", "pushURL"),
        auth=(
            config.get("parkour", "user"),
            config.get("parkour", "password")
        ),
        data=d
    )
    log.info("ParkourPush return {}".format(pushParkStat))
    return pushParkStat


def multiQC_yaml(config, flowcell, ssDic, project, laneFolder):
    '''
    This function creates:
     - config yaml, containing appropriate header information
     - data string adding gen stats
     - data string containing our old seqreport statistics.
    Keep in mind we delete these after running mqc
    '''
    ssdf = ssDic['sampleSheet'][
        ssDic['sampleSheet']['Sample_Project'] == project
    ].fillna('NA')

    # data string genstats
    mqcData = "# format: 'tsv'\n"
    mqcData += "# plot_type: 'generalstats'\n"
    mqcData += "# pconfig: \n"
    mqcData += "Sample_Name\tSample_ID\tRequested\n"
    reqDict = {}
    reqsMax = 0
    for sample in list(ssdf['Sample_Name'].unique()):
        sampleID = ssdf[ssdf['Sample_Name'] == sample]['Sample_ID'].values[0]
        if ssdf[ssdf['Sample_Name'] == sample]['reqDepth'].values[0] == 'NA':
            reqDepth = 'NA'
        else:
            reqDepth = float(
                ssdf[ssdf['Sample_Name'] == sample]['reqDepth'].values[0]
            )
        if reqDepth != 'NA':
            if reqDepth > reqsMax:
                reqsMax = reqDepth
        sampleLis = glob.glob(
            os.path.join(
                laneFolder,
                '*/*/' + sample + '*fastq.gz'
            )
        )
        purgeSampleLis = []
        for i in sampleLis:
            if 'optical' not in i:
                purgeSampleLis.append(i)
        for fullfqFile in purgeSampleLis:
            fqFile = fullfqFile.split('/')[-1]
            sampleBase = fqFile.replace(".fastq.gz", "")
            reqDict[sampleBase] = [sampleID, reqDepth]
    for sample in reqDict:
        mqcData += "{}\t{}\t{}\n".format(
            sample, reqDict[sample][0], reqDict[sample][1]
        )
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
        ssdf['Library_Type'].unique()
    ))
    # indexTypes
    ixTypes = ', '.join(list(
        ssdf["indexType"].unique()
    ))
    # Protocols
    protTypes = ', '.join(list(
        ssdf["Description"].unique()
    ))
    # Organisms
    orgs = ', '.join(list(
        ssdf["Organism"].unique()
    ))
    # Resequencing runs are screwed up (e.g. don't contain the samples)
    # switch total requested to NA
    try:
        sumReqRound = str(
            round(ssdf['reqDepth'].sum(), 0)
        )
    except TypeError:
        sumReqRound = 'NA'
    sumReqRound 
    mqcyml = {
        "title": project,
        "custom_logo": config["misc"]["mpiImg"],
        "custom_logo_url": "https://www.ie-freiburg.mpg.de/",
        "custom_logo_title": "MPI-IE",
        "show_analysis_paths": False,
        "show_analysis_time": False,
        "fastqscreen_simpleplot": False,
        "log_filesize_limit": 2000000000,
        "report_header_info": [
            {"Contact E-mail": config["communication"]["bioinfoCore"]},
            {"Flowcell": flowcell.name},
            {"Sequencer Type": flowcell.sequencer},
            {"Read Lengths": formatSeqRecipe(flowcell.seqRecipe)},
            {"Demux. Mask": ssDic["mask"]},
            {"Mismatches": formatMisMatches(ssDic["mismatch"])},
            {"dissectBCL version": "0.0.1"},
            {"bcl-convert version": config["softwareVers"]["bclconvertVer"]},
            {"Library Type": libTypes},
            {"Library Protocol": protTypes},
            {"Index Type": ixTypes},
            {"Organism": orgs},
            {"Requested reads": sumReqRound},
            {"Received reads": str(
                round(ssdf['gotDepth'].replace('NA', np.nan).dropna().sum(), 0)
                )}
        ]
    }
    return(mqcyml, mqcData, seqreportData, indexreportData)


def mailHome(subject, _html, config, toCore=False):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = '[dissectBCL] [v0.0.1] ' + subject
    mailer['From'] = config['communication']['fromAddress']
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
    outLane = outPath.split('/')[-1]
    # Get directories from outPath.
    for projectPath in glob.glob(
        os.path.join(
            outPath,
            'Project*'
        )
    ):
        project = projectPath.split('/')[-1]
        shipDic[project] = 'No'
        log.info("Shipping {}".format(project))
        PI = project.split('_')[-1].lower().replace(
            "cabezas-wallscheid", "cabezas"
        )
        fqcPath = projectPath.replace("Project_", "FASTQC_Project_")
        if PI in config['Internals']['PIs']:
            log.info("Found {}. Shipping internally.".format(PI))
            fqc = fqcPath.split('/')[-1]
            enduserBase = os.path.join(
                fetchLatestSeqDir(
                    os.path.join(
                        config['Dirs']['piDir'],
                        PI
                    ),
                    config['Internals']['seqDir']
                ),
                outLane
            )
            if not os.path.exists(enduserBase):
                os.mkdir(enduserBase, mode=0o750)
            copyStat_FQC = ""
            copyStat_Project = ""
            if os.path.exists(
                os.path.join(enduserBase, fqc)
            ):
                shutil.rmtree(
                    os.path.join(
                        enduserBase, fqc
                    )
                )
                copyStat_FQC = "Replaced"
            shutil.copytree(
                fqcPath,
                os.path.join(
                    enduserBase,
                    fqc
                )
            )
            if copyStat_FQC != "Replaced":
                copyStat_FQC = "Copied"
            if os.path.exists(
                os.path.join(enduserBase, project)
            ):
                shutil.rmtree(
                    os.path.join(
                        enduserBase,
                        project
                    )
                )
                copyStat_Project = "Replaced"
            shutil.copytree(
                projectPath,
                os.path.join(
                    enduserBase,
                    project
                )
            )
            log.info("Stripping group rights for {}".format(enduserBase))
            for r, dirs, files in os.walk(enduserBase):
                for d in dirs:
                    os.chmod(os.path.join(r, d), 0o700)
                for f in files:
                    os.chmod(os.path.join(r, f), 0o700)
            if copyStat_Project != "Replaced":
                copyStat_Project = "Copied"
            if 'Replaced' in [copyStat_Project, copyStat_FQC]:
                shipDic[project] = 'Transferred'
            else:
                shipDic[project] = 'Copied'
        else:
            shipDicStat = "Uploaded"
            laneStr = fqcPath.split('/')[-2]
            # If the same tarball is already present, replace it.
            fexList = check_output(
                [
                    'fexsend',
                    '-l',
                    config['communication']['fromAddress']
                ]
            ).decode("utf-8").replace("\n", "").split(' ')
            tarBall = laneStr + '_' + project + '.tar'
            if tarBall in fexList:
                fexRm = [
                    'fexsend',
                    '-d',
                    tarBall,
                    config['communication']['fromAddress']
                ]
                fexdel = Popen(fexRm)
                fexdel.wait()
                shipDicStat = "Replaced"
            fexer = "tar cf - {} {} | fexsend -s {}.tar {}".format(
                projectPath,
                fqcPath,
                laneStr + '_' + project,
                config['communication']['fromAddress']
            )
            os.system(fexer)
            shipDic[project] = shipDicStat
    # Ship multiQC reports.
    seqFacDir = os.path.join(
        config['Dirs']['seqFacDir'],
        outLane
    )
    if not os.path.exists(seqFacDir):
        os.mkdir(seqFacDir)
    for qcRepo in glob.glob(
        os.path.join(outPath, 'Project_*', 'multiqc_report.html')
    ):
        # to seqfacdir
        outqcRepo = os.path.join(
            seqFacDir, qcRepo.split('/')[-2] + '_multiqcreport.html'
        )
        shutil.copyfile(qcRepo, outqcRepo)
        # to bioinfoCoredir
        outqcBioinfo = os.path.join(
            config['Dirs']['bioinfoCoreDir'],
            qcRepo.split('/')[-2] + '_multiqcreport.html'
        )
        shutil.copyfile(qcRepo, outqcBioinfo)
    transferStop = datetime.datetime.now()
    transferTime = transferStop - transferStart
    return(transferTime, shipDic)


def organiseLogs(flowcell, sampleSheet):
    for outLane in sampleSheet.ssDic:
        log.info("Populating log dir for {}".format(outLane))
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
        # Write out the yaml files.
        yaml = ruamel.yaml.YAML()
        yaml.indent(mapping=2, sequence=4, offset=2)

        # write out outLaneInfo.yaml
        outLaneInfo = os.path.join(_logDir, 'outLaneInfo.yaml')
        dic = sampleSheet.ssDic[outLane]
        del dic['sampleSheet']
        with open(outLaneInfo, 'w') as f:
            ruamel.yaml.dump(dic, f)
        # write out flowcellInfo.yaml
        flowcellInfo = os.path.join(_logDir, 'flowcellInfo.yaml')
        dic = flowcell.asdict()
        with open(flowcellInfo, 'w') as f:
            ruamel.yaml.dump(dic, f)
