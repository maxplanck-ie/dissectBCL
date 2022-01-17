import datetime
import requests
import pandas as pd
from dissectBCL.logger import log
from dissectBCL.misc import joinLis, retBCstr, retIxtype
from subprocess import Popen
import os
import shutil


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


def greeter():
    now = datetime.datetime.now()
    if now.hour < 12:
        return "Good Morning!"
    elif now.hour < 18:
        return "Good Afternoon!"
    else:
        return "Good Evening!"


def buildTexTable(PEstatus, df):
    if PEstatus:
        texTable = r'''
\begin{center}
\begin{tabular}{|c | c | c | c | c | c | c | c | c |}
\hline
SampleID & Sample Name & Barcode(s) & Barcode ID(s) & Lane(s) & Fragments (M) & Frag/Req & mean Q PF & perc Q30\\
\hline\hline'''
        for index, row in df.iterrows():
            texTable = texTable + r'''
%(samID)s & %(samName)s & %(BC)s & %(BCID)s & %(lane)s & %(reads)s & %(readvreq)s & %(meanQ)s & %(perc30)s\\
\hline''' % {
                'samID': row['Sample_ID'],
                'samName': row['Sample_Name'].replace('_', '\_'),
                'BC': retBCstr(row),
                'BCID': retIxtype(row),
                'lane': int(row['Lane']),
                'reads': row['gotDepth'],
                'readvreq': row['gotDepth']/row['reqDepth'],
                'meanQ': row['meanQ'],
                'perc30': row['percQ30'],
            }
        texTable += r'''
        \end{tabular}
        \end{center}
        '''
        return texTable

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
    if os.path.exists(absOutPdf):
        os.remove(absOutPdf)
    if not os.path.exists(absOutPdf):
        ss = ssdf[ssdf['Sample_Project'] == project]
        libTypes = ','.join(
            list(ss['Library_Type'].unique())
        )
        Protocol = ','.join(
            list(ss['Description'].unique())
        )
        indexType = ','.join(
            list(ss['indexType'].unique())
        )
        # Read up the tex template.
        templatePath = os.path.join(
            os.path.dirname(__file__),
            'templates',
            'SequencingReport.tex'
        )
        with open(templatePath) as f:
            txTemp = f.read()
        txTemp = txTemp % {
            'project': project.replace('_', '\_'),
            'date': str(datetime.datetime.now().replace(microsecond=0)),
            'flowcellname': flowcell.name.replace('_', '\_'),
            'flowcellsequencer': flowcell.sequencer,
            'readlen': ';'.join([str(x[-1]) for x in list(flowcell.seqRecipe.values())]),
            'mask': sampleSheet.ssDic[outLane]['mask'],
            'vers': '0.0.1',
            'convvers' : config['softwareVers']['bclconvertVer'],
            'mismatch' : joinLis(list(sampleSheet.ssDic[outLane]['mismatch'].values()), joinStr=", "),
            'libtyp' : libTypes,
            'ixtyp' : indexType,
            'prot': Protocol,
            'texTable': buildTexTable(sampleSheet.ssDic[outLane]['PE'], ssdf)
        }
        with open(absOutTex, 'w') as f:
            f.write(txTemp)
        pdfProc = Popen(['tectonic', absOutTex,'--keep-logs', '--outdir', outDir])
        pdfProc.wait()
        #os.remove(absOutTex)
        shutil.copy(absOutPdf, '/data/manke/group/deboutte/qc/SEQREP.pdf')


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
                ssdf,
                config,
                flowcell,
                outLane,
                sampleSheet
            )
