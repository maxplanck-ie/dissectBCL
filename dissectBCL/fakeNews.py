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
    return pd.DataFrame()


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
    if ssDic['PE']:
        ssdf['reqDepth/2'] = ssdf['reqDepth']/2

    # data string genstats
    mqcData = "# format: 'tsv'\n"
    mqcData += "# plot_type: 'generalstats'\n"
    mqcData += "# pconfig: \n"
    mqcData += "Sample ID\tRequested\n"
    reqDict = {}
    reqsMax = 0
    for sample in list(ssdf['Sample_Name'].unique()):
        sampleID = ssdf[ssdf['Sample_Name'] == sample]['Sample_ID'].values[0]
        if ssDic['PE']:
            reqDepth = float(
                ssdf[ssdf['Sample_Name'] == sample]['reqDepth/2'].values[0]
            )
        else:
            reqDepth = float(
                ssdf[ssdf['Sample_Name'] == sample]['reqDepth'].values[0]
            )
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
    seqreportData = ""
    for index, row in ssdf.iterrows():
        if seqreportData == "":
            meanQ_headers, Meanq = retMean_perc_Q(row, returnHeader=True)
            percq30_headers, perc30 = retMean_perc_Q(
                row, returnHeader=True, qtype='percQ30'
            )
            seqreportData += \
                "\tSample ID\tBarcodes\tindexTypes\tLane\t{}\t{}\n".format(
                    meanQ_headers,
                    percq30_headers
                )
            seqreportData += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row['Sample_Name'],
                row['Sample_ID'],
                retBCstr(row),
                retIxtype(row),
                "L" + str(row['Lane']),
                Meanq,
                perc30
            )
        else:
            Meanq = retMean_perc_Q(row)
            perc30 = retMean_perc_Q(row, qtype='percQ30')
            seqreportData += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                row['Sample_Name'],
                row['Sample_ID'],
                retBCstr(row),
                retIxtype(row),
                "L" + str(row['Lane']),
                Meanq,
                perc30
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
    mqcyml = {
        "title": project,
        "intro_text": "This is a placeholder.",
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
            {"Organism": orgs}
        ]
    }
    return(mqcyml, mqcData, seqreportData)


def mailHome(subject, _html, config):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = '[dissectBCL] [v0.0.1] ' + subject
    mailer['From'] = config['communication']['fromAddress']
    mailer['To'] = config['communication']['bioinfoCore']
    email = MIMEText(_html, 'html')
    mailer.attach(email)
    s = smtplib.SMTP(config['communication']['host'])
    s.sendmail(
        config['communication']['fromAddress'],
        config['communication']['bioinfoCore'],
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
        if PI in config['Internals']['PIs']:
            log.info("Found {}. Shipping internally.".format(PI))
            fqcPath = projectPath.replace("Project_", "FASTQC_Project_")
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
            copyStat_FQC = False
            copyStat_Project = False
            if not os.path.exists(
                os.path.join(enduserBase, fqc)
            ):
                shutil.copytree(
                    fqcPath,
                    os.path.join(
                        enduserBase,
                        fqc
                    )
                )
                copyStat_FQC = True
            if not os.path.exists(
                os.path.join(enduserBase, project)
            ):
                shutil.copytree(
                    projectPath,
                    os.path.join(
                        enduserBase,
                        project
                    )
                )
            if copyStat_FQC and copyStat_Project:
                shipDic[project] = 'Transferred'
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
