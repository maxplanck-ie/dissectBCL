import os
import glob
from dissectBCL.fakeNews import multiQC_yaml, mailHome
from dissectBCL.misc import krakenfqs, moveOptDup, validateFqEnds
from pathlib import Path
import re
import shutil
import sys
from multiprocessing import Pool
from pandas import isna
from subprocess import Popen, DEVNULL
import ruamel.yaml
import logging


def matchIDtoName(ID, ssdf):

    if (ID not in set(ssdf['Sample_ID'])):
        # can happen if filename is not legit
        # e.g. if demuxSheet did not match SampleSheet
        logging.critical(
            "ID {} is not defined in SampleSheet.".format(ID)
        )
        sys.exit(1)

    name = ssdf[ssdf['Sample_ID'] == ID]['Sample_Name'].values

    if len(name) > 1:
        # It can happen one sample sits in 2 lanes.
        if len(set(name)) > 1:
            logging.cricital(
                "SampleID {} has multiple names {}, exiting.".format(
                    ID, name
                )
            )
            sys.exit(1)

        # can happen if ID is not listed in SampleSheet --> no Sample_Name
        elif isna(name[0]):
            logging.critical(
                "Sample_Name is not defined for ID {} .".format(ID)
            )
            sys.exit(1)
        else:
            return name[0]
    else:
        return name[0]


def renamefq(fqFile, projectFolder, ssdf, laneSplitStatus):
    oldName = fqFile.split('/')[-1]
    sampleID = oldName.split('_')[0]
    sampleName = matchIDtoName(sampleID, ssdf)
    sampleIDPath = os.path.join(
        projectFolder,
        "Sample_" + sampleID
    )
    if not os.path.exists(sampleIDPath):
        os.mkdir(sampleIDPath)
    # Create new name
    if laneSplitStatus:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+L[0-9][0-9][0-9]_+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(
            regstr,
            r'_\1',
            newName
        )
        newName.replace(sampleID, sampleName)
    else:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(
            regstr,
            r'_\1',
            newName
        )
    return os.path.join(sampleIDPath, newName)


def renameProject(projectFolder, ssdf, laneSplitStatus):
    """
    rename and move files under sample_ID folders.
    rename project folder from e.g.
    1906_Hein_B03_Hein -> Project_1906_Hein_B03_Hein
    """
    logging.info("Renaming {}".format(projectFolder))
    for fq in glob.glob(
        os.path.join(
            projectFolder,
            '*fastq.gz'
        )
    ):
        newName = renamefq(fq, projectFolder, ssdf, laneSplitStatus)
        shutil.move(fq, newName)
    # Finally rename the project folder.
    projID = projectFolder.split('/')[-1]
    shutil.move(
        projectFolder,
        projectFolder.replace(projID, "Project_" + projID)
    )


def fqcRunner(cmd):
    cmds = cmd.split(" ")
    qcRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = qcRun.wait()
    return (exitcode)


def qcs(project, laneFolder, sampleIDs, config):
    # make fastqc folder.
    fqcFolder = os.path.join(
        laneFolder,
        "FASTQC_Project_" + project
    )
    if not os.path.exists(fqcFolder):
        os.mkdir(fqcFolder)
    fastqcCmds = []
    for ID in sampleIDs:
        # Colliding samples are omitted, and don't have a folder.
        if not os.path.exists(
            os.path.join(
                laneFolder,
                'Project_' + project,
                'Sample_' + ID
            )
        ):
            continue
        else:
            IDfolder = os.path.join(
                laneFolder,
                "FASTQC_Project_" + project,
                "Sample_" + ID
            )
            if not os.path.exists(IDfolder):
                os.mkdir(IDfolder)
            fqFiles = glob.glob(
                os.path.join(
                    laneFolder,
                    "Project_" + project,
                    "Sample_" + ID,
                    '*fastq.gz'
                )
            )
            # Don't do double work.
            if len(glob.glob(
                os.path.join(IDfolder, "*zip")
            )) == 0:
                fastqcCmds.append(
                    " ".join([
                        'fastqc',
                        '-a',
                        config['software']['fastqc_adapters'],
                        '-q',
                        '-t',
                        str(len(fqFiles)),
                        '-o',
                        IDfolder
                    ] + fqFiles)
                )
    if fastqcCmds:
        logging.info(
            "fastqc example for {} - {}".format(
                project, fastqcCmds[0]
            )
        )
        with Pool(20) as p:
            fqcReturns = p.map(fqcRunner, fastqcCmds)
            if fqcReturns.count(0) == len(fqcReturns):
                logging.info("FastQC done for {}.".format(project))
            else:
                logging.critical(
                    "FastQC runs failed for {}. exiting.".format(project)
                )
                mailHome(
                    laneFolder,
                    "FastQC runs failed. Investigate.",
                    config,
                    toCore=True
                )
                sys.exit(1)
    else:
        logging.info("FastQCs already done for {}".format(project))


def clmpRunner(cmd):
    cmds = cmd.split(" ")
    baseName = cmds.pop(-1)
    PE = str(cmds.pop(-1))
    samplePath = cmds.pop(-1)
    os.chdir(samplePath)
    clumpRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = clumpRun.wait()
    splitCmd = ['splitFastq', 'tmp.fq.gz', PE, baseName, '10']
    splitFq = Popen(splitCmd, stdout=DEVNULL, stderr=DEVNULL)
    exitcode_split = splitFq.wait()
    os.remove('tmp.fq.gz')
    return (
        (exitcode, exitcode_split)
    )


def clumper(project, laneFolder, sampleIDs, config, PE, sequencer):
    logging.info("Clump for {}".format(project))
    clmpOpts = {
        'general': [
            'out=tmp.fq.gz',
            'dupesubs=0',
            'qin=33',
            'markduplicates=t',
            'optical=t',
            '-Xmx400G',
            'threads=15',
            'tmpdir={}'.format(config['Dirs']['tempDir'])
        ],
        'NextSeq': [
            'spany=t',
            'adjacent=t',
            'dupedist=40'
        ],
        'NovaSeq': ['dupedist=12000']
    }
    clmpCmds = []
    if sequencer != 'MiSeq':
        for ID in sampleIDs:
            sampleDir = os.path.join(
                laneFolder,
                "Project_" + project,
                "Sample_" + ID
            )
            if os.path.exists(sampleDir):
                if len(glob.glob(
                    os.path.join(sampleDir, '*optical_duplicates.fastq.gz')
                )) == 0:
                    fqFiles = glob.glob(
                        os.path.join(
                            sampleDir,
                            "*fastq.gz"
                        )
                    )
                    if len(fqFiles) < 3:
                        if PE and len(fqFiles) == 2:
                            for i in fqFiles:
                                if '_R1.fastq.gz' in i:
                                    in1 = "in=" + i
                                    baseName = i.split('/')[-1].replace(
                                        "_R1.fastq.gz",
                                        ""
                                    )
                                elif '_R2.fastq.gz' in i:
                                    in2 = "in2=" + i
                            clmpCmds.append(
                                'clumpify.sh' + " " +
                                in1 + " " +
                                in2 + " " +
                                " ".join(clmpOpts['general']) + " " +
                                " ".join(clmpOpts[sequencer]) + " " +
                                sampleDir + " " +
                                "1" + " " +
                                baseName
                            )
                        elif not PE and len(fqFiles) == 1:
                            if '_R1.fastq.gz' in fqFiles[0]:
                                in1 = "in=" + fqFiles[0]
                                baseName = fqFiles[0].split('/')[-1].replace(
                                    "_R1.fastq.gz",
                                    ""
                                )
                                clmpCmds.append(
                                    'clumpify.sh' + " " +
                                    in1 + " " +
                                    " ".join(clmpOpts['general']) + " " +
                                    " ".join(clmpOpts[sequencer]) + " " +
                                    sampleDir + " " +
                                    "0" + " " +
                                    baseName
                                )
                            else:
                                logging.info("Not clumping {}".format(ID))
        if clmpCmds:
            logging.info(
                "clump example for {} - {}".format(
                    project,
                    clmpCmds[0]
                )
            )
            with Pool(5) as p:
                clmpReturns = p.map(clmpRunner, clmpCmds)
                if clmpReturns.count((0, 0)) == len(clmpReturns):
                    logging.info("Clumping done for {}.".format(project))
                else:
                    logging.critical(
                        "Clumping failed for {}. exiting.".format(project)
                    )
                    mailHome(
                        laneFolder,
                        "Clump runs failed. Investigate.",
                        config,
                        toCore=True
                    )
                    sys.exit(1)
        else:
            logging.info("No clump run for {}".format(project))
    else:
        logging.info("No clump for MiSeq.")


def krakRunner(cmd):
    cmds = cmd.split(" ")
    krakRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = krakRun.wait()
    return exitcode


def kraken(project, laneFolder, sampleIDs, config):
    logging.info("Kraken for {}".format(project))
    krakenCmds = []
    for ID in sampleIDs:
        IDfolder = os.path.join(
            laneFolder,
            "FASTQC_Project_" + project,
            "Sample_" + ID
        )
        if os.path.exists(IDfolder):
            if len(glob.glob(
                os.path.join(
                    IDfolder,
                    '*.rep'
                )
            )) == 0:
                sampleFolder = os.path.join(
                    laneFolder,
                    "Project_" + project,
                    "Sample_" + ID
                )
                reportname, fqs = krakenfqs(sampleFolder)
                krakenCmds.append(
                    ' '.join([
                        'kraken2',
                        '--db',
                        config['software']['kraken2db'],
                        '--out',
                        '-',
                        '--threads',
                        '4',
                        '--report',
                        reportname
                    ] + fqs)
                )
    if krakenCmds:
        logging.info(
            "kraken example for {} - {}".format(
                project,
                krakenCmds[0]
            )
        )
        with Pool(10) as p:
            screenReturns = p.map(krakRunner, krakenCmds)
            if screenReturns.count(0) == len(screenReturns):
                logging.info("kraken ran {}.".format(project))
            else:
                logging.critical(
                    "kraken failed for {}. exiting.".format(project)
                )
                mailHome(
                    laneFolder,
                    "kraken runs failed. Investigate.",
                    config,
                    toCore=True
                )
                sys.exit(1)
    else:
        logging.info(
            "kraken files already present. Skipping {}".format(project)
        )


def multiqc(project, laneFolder, config, flowcell, sampleSheet):
    logging.info("multiqc for {}".format(project))
    outLane = laneFolder.split('/')[-1]
    QCFolder = os.path.join(
        laneFolder,
        'FASTQC_Project_' + project
    )
    projectFolder = os.path.join(
        laneFolder,
        "Project_" + project
    )
    # md5sums
    md5CmdStr = \
        r"md5sum {} | {} | {} | {} > {}".format(
            projectFolder + '/*/*fastq.gz',
            r"sed 's/  //g'",
            r"cut -d '/' -f1,8",
            r"sed 's/\//\t/g'",
            os.path.join(projectFolder, 'md5sums.txt')
            )
    if not os.path.exists(
        os.path.join(projectFolder, 'md5sums.txt')
    ):
        os.system(md5CmdStr)
    # Always overwrite the multiQC reports. RunTimes are marginal anyway.
    mqcConf, mqcData, seqrepData, indexreportData = multiQC_yaml(
        config,
        flowcell,
        sampleSheet.ssDic[outLane],
        project,
        laneFolder
    )
    yaml = ruamel.yaml.YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    confOut = os.path.join(
        projectFolder,
        'mqc.yaml'
    )
    dataOut = os.path.join(
        QCFolder,
        'parkour_mqc.tsv'
    )
    seqrepOut = os.path.join(
        QCFolder,
        'Sequencing_Report_mqc.tsv'
    )
    indexrepOut = os.path.join(
        QCFolder,
        'Index_Info_mqc.tsv'
    )
    with open(confOut, 'w') as f:
        yaml.dump(mqcConf, f)
    with open(seqrepOut, 'w') as f:
        f.write(seqrepData)
    with open(dataOut, 'w') as f:
        f.write(mqcData)
    with open(indexrepOut, 'w') as f:
        f.write(indexreportData)
    multiqcCmd = [
        'multiqc',
        '--quiet',
        '--no-data-dir',
        '-f',
        '-o',
        projectFolder,
        '-c',
        confOut,
        QCFolder
    ]
    multiqcRun = Popen(multiqcCmd, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = multiqcRun.wait()
    if exitcode == 0:
        logging.info('multiqc ran for {}'.format(project))
        os.remove(confOut)
        os.remove(dataOut)
        os.remove(seqrepOut)
        os.remove(indexrepOut)
    else:
        logging.critical("multiqc failed for {}".format(project))
        mailHome(
            laneFolder,
            "multiQC runs failed. Investigate.",
            config,
            toCore=True
        )
        sys.exit(1)
    return exitcode


def postmux(flowcell, sampleSheet, config):
    logging.warning("Postmux module")
    for outLane in sampleSheet.ssDic:
        laneFolder = os.path.join(flowcell.outBaseDir, outLane)
        # Don't rename if renamed.done flag is there.
        renameFlag = os.path.join(
            laneFolder,
            'renamed.done'
        )
        postmuxFlag = os.path.join(
            laneFolder,
            'postmux.done'
        )
        df = sampleSheet.ssDic[outLane]['sampleSheet']
        projects = list(df['Sample_Project'].unique())
        # Renaming module
        if not os.path.exists(renameFlag):
            logging.info("Renaming {}".format(outLane))
            for project in projects:
                renameProject(
                    os.path.join(
                        flowcell.outBaseDir,
                        outLane,
                        project
                    ),
                    df,
                    sampleSheet.laneSplitStatus
                )
            # Sanity check for file endings.
            fqEnds = validateFqEnds(laneFolder)
            if not fqEnds:
                logging.info("All fastq files have proper extension.")
            else:
                logging.critical(
                    "some fastq files aren't properly renamed: {}".format(
                        fqEnds
                    )
                )
                sys.exit(
                    "some fastq files aren't properly renamed: {}".format(
                        fqEnds
                    )
                )
            Path(
                os.path.join(laneFolder, 'renamed.done')
            ).touch()
        # postMux module
        if not os.path.exists(postmuxFlag):
            logging.info("FastQC pool {}".format(outLane))
            for project in projects:
                # FastQC
                qcs(
                    project,
                    laneFolder,
                    set(
                        df[df['Sample_Project'] == project]['Sample_ID']
                    ),
                    config
                )
                # clump
                clumper(
                    project,
                    laneFolder,
                    set(
                        df[df['Sample_Project'] == project]['Sample_ID']
                    ),
                    config,
                    sampleSheet.ssDic[outLane]['PE'],
                    flowcell.sequencer
                )
                # kraken
                kraken(
                    project,
                    laneFolder,
                    set(
                        df[df['Sample_Project'] == project]['Sample_ID']
                    ),
                    config
                )
                # multiQC
                multiqc(
                    project,
                    laneFolder,
                    config,
                    flowcell,
                    sampleSheet
                )
            logging.info("Moving optical dup txt into FASTQC folder")
        moveOptDup(laneFolder)
        Path(
                os.path.join(laneFolder, 'postmux.done')
        ).touch()
    return (0)
