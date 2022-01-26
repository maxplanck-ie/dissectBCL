import datetime
import requests
import pandas as pd
from dissectBCL.logger import log
from dissectBCL.misc import joinLis, retBCstr, retIxtype
from dissectBCL.misc import TexformatQual, TexformatDepFrac
from dissectBCL.misc import ReportDFSlicer, truncStr
from dissectBCL.misc import fetchLatestSeqDir
from subprocess import Popen, DEVNULL
import os
import shutil
from random import randint
import smtplib
import datetime
import glob


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
                flatLis.append([project, sample] + res.json()[project][sample])
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


def buildTexTable(PEstatus, df):
    if PEstatus:
        headPath = os.path.join(
            os.path.dirname(__file__),
            'templates',
            'PEhead.tex'
        )
        tablePath = os.path.join(
            os.path.dirname(__file__),
            'templates',
            'PEtable.tex'
        )
        # less then 25 samples = 1 page.
        if len(df.index) < 25:
            with open(headPath) as f:
                txTable = f.read()
            with open(tablePath) as f:
                txTemplate = f.read()
            for index, row in df.iterrows():
                txTable += txTemplate % {
                    'samID': row['Sample_ID'],
                    'samName': truncStr(row['Sample_Name'].replace('_', r'\_')),
                    'BC': retBCstr(row),
                    'BCID': retIxtype(row),
                    'lane': row['Lane'],
                    'reads': row['gotDepth'],
                    'readvreq': TexformatDepFrac(
                        row['gotDepth']/row['reqDepth']
                    ),
                    'meanQ': TexformatQual(row['meanQ']),
                    'perc30': TexformatQual(row['percQ30']),
                }
            txTable += r'''
            \end{tabular}
            \end{center}
            '''
            return txTable
        elif len(df.index) > 25:
            slices = ReportDFSlicer(len(df.index))
            multiTable = r''
            for slice in slices:
                with open(headPath) as f:
                    txTable = f.read()
                with open(tablePath) as f:
                    txTemplate = f.read()
                for index, row in df.iloc[slice[0]:slice[1]].iterrows():
                    txTable += txTemplate % {
                        'samID': row['Sample_ID'],
                        'samName': truncStr(row['Sample_Name'].replace('_', r'\_')),
                        'BC': retBCstr(row),
                        'BCID': retIxtype(row),
                        'lane': row['Lane'],
                        'reads': row['gotDepth'],
                        'readvreq': TexformatDepFrac(
                            row['gotDepth']/row['reqDepth']
                        ),
                        'meanQ': TexformatQual(row['meanQ']),
                        'perc30': TexformatQual(row['percQ30']),
                    }
                txTable += r'''
                \end{tabular}
                \end{center}
                \newpage
                '''
                multiTable += txTable
            return multiTable
                


def buildSeqReport(project, ssdf, config, flowcell, outLane, sampleSheet):
    absOutTex = os.path.join(
        flowcell.outBaseDir,
        outLane,
        'Project_' + project,
        'SequencingReport.tex'
    )
    absOutPdf = absOutTex.replace("tex", 'pdf')
    outDir = os.path.join(
        flowcell.outBaseDir,
        outLane,
        'Project_' + project
    )
    # Always rebuild report.
    if os.path.exists(absOutPdf):
        os.remove(absOutPdf)
    ss = ssdf[ssdf['Sample_Project'] == project]
    libTypes = ','.join(
        [str(x) for x in list(ss['Library_Type'].unique())]
    ).replace('_', r'\_')
    Protocol = ','.join(
        [str(x) for x in list(ss['Description'].unique())]
    ).replace('_', r'\_')
    indexType = ','.join(
        [str(x) for x in list(ss['indexType'].unique())]
    ).replace('_', r'\_')
    # Read up the tex template.
    templatePath = os.path.join(
        os.path.dirname(__file__),
        'templates',
        'SequencingReport.tex'
    )
    with open(templatePath) as f:
        txTemp = f.read()
    txTemp = txTemp % {
        'project': project.replace('_', r'\_'),
        'date': str(datetime.datetime.now().replace(microsecond=0)),
        'flowcellname': flowcell.name.replace('_', r'\_'),
        'flowcellsequencer': flowcell.sequencer,
        'readlen': ';'.join(
            [str(x[-1]) for x in list(flowcell.seqRecipe.values())]
        ),
        'mask': sampleSheet.ssDic[outLane]['mask'],
        'vers': '0.0.1',
        'convvers': config['softwareVers']['bclconvertVer'],
        'mismatch': joinLis(
            list(sampleSheet.ssDic[outLane]['mismatch'].values()),
            joinStr=", "
        ),
        'libtyp': libTypes,
        'ixtyp': indexType,
        'prot': Protocol,
        'texTable': buildTexTable(sampleSheet.ssDic[outLane]['PE'], ssdf)
    }
    with open(absOutTex, 'w') as f:
        f.write(txTemp)
    pdfProc = Popen(['tectonic', absOutTex, '--outdir', outDir], stdout=DEVNULL, stderr=DEVNULL)
    pdfProc.wait()
    #os.remove(absOutTex)
    log.info("Attempting copy")
    shutil.copy(
        absOutPdf,
        os.path.join(
            config['Dirs']['bioinfoCoreDir'],
            project + '_seqrep.pdf'
        )
    )


def runSeqReports(flowcell, sampleSheet, config):
    log.info("Building sequencing Reports")
    for outLane in sampleSheet.ssDic:
        ssdf = sampleSheet.ssDic[outLane]['sampleSheet']
        projects = list(
            ssdf['Sample_Project'].unique()
        )
        for project in projects:
            buildSeqReport(
                project,
                ssdf[ssdf['Sample_Project'] == project],
                config,
                flowcell,
                outLane,
                sampleSheet
            )


def mailHome(subject, msg, config):
    msg['Subject'] = subject
    msg['From'] = config['communication']['fromAddress']
    msg['To'] = config['communication']['bioinfoCore']
    s = smtplib.SMTP(config['communication']['host'])
    s.send_message(msg)
    s.quit()


def shipFiles(outPath, config):
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
        shipDic[projectPath] = []
        log.info("Shipping {}".format(project))
        PI = project.split('_')[-1].lower().replace("cabezas-wallscheid", 'cabezas')
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





