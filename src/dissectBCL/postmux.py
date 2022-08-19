import os
import glob
from dissectBCL.fakeNews import multiQC_yaml, mailHome
from dissectBCL.logger import log
from dissectBCL.misc import screenFqFetcher, moveOptDup
from pathlib import Path
import re
import shutil
import sys
from multiprocessing import Pool
from pandas import isna
from subprocess import Popen, DEVNULL
import ruamel.yaml


def matchIDtoName(ID, ssdf):

    if (ID not in set(ssdf['Sample_ID'])):
        # can happen if filename is not legit
        # e.g. if demuxSheet did not match SampleSheet
        log.critical(
            "ID {} is not defined in SampleSheet.".format(ID)
        )
        sys.exit(1)

    name = ssdf[ssdf['Sample_ID'] == ID]['Sample_Name'].values

    if len(name) > 1:
        # It can happen one sample sits in 2 lanes.
        if len(set(name)) > 1:
            log.cricital(
                "SampleID {} has multiple names {}, exiting.".format(
                    ID, name
                )
            )
            sys.exit(1)

        # can happen if ID is not listed in SampleSheet --> no Sample_Name
        elif isna(name[0]):
            log.critical(
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
        newName = re.sub(
            r"_S[0-9]?[0-9]_+L[0-9][0-9][0-9]_+([IR][123])+_[0-9][0-9][0-9]",
            r'_\1',
            newName
        )
        newName.replace(sampleID, sampleName)
    else:
        newName = oldName.replace(sampleID, sampleName)
        newName = re.sub(
            r"_S[0-9]?[0-9]_+([IR][123])+_[0-9][0-9][0-9]",
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
    log.info("Renaming {}".format(projectFolder))
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
                    config['software']['fastqc'],
                    '-q',
                    '-t',
                    str(len(fqFiles)),
                    '-o',
                    IDfolder
                ] + fqFiles)
            )
    if fastqcCmds:
        with Pool(20) as p:
            fqcReturns = p.map(fqcRunner, fastqcCmds)
            if fqcReturns.count(0) == len(fqcReturns):
                log.info("FastQC done for {}.".format(project))
            else:
                log.critical(
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
        log.info("FastQCs already done for {}".format(project))


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
    log.info("Clump for {}".format(project))
    clmpOpts = {
        'general': [
            'out=tmp.fq.gz',
            'dupesubs=0',
            'qin=33',
            'markduplicates=t',
            'optical=t',
            '-Xmx400G',
            'threads=15',
            'tmpdir=/scratch/local'
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
                            config['software']['clumpify'] + " " +
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
                                config['software']['clumpify'] + " " +
                                in1 + " " +
                                " ".join(clmpOpts['general']) + " " +
                                " ".join(clmpOpts[sequencer]) + " " +
                                sampleDir + " " +
                                "0" + " " +
                                baseName
                            )
                        else:
                            log.info("Not clumping {}".format(ID))
        if clmpCmds:
            with Pool(5) as p:
                clmpReturns = p.map(clmpRunner, clmpCmds)
                if clmpReturns.count((0, 0)) == len(clmpReturns):
                    log.info("Clumping done for {}.".format(project))
                else:
                    log.critical(
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
            log.info("No clump run for {}".format(project))
    else:
        log.info("No clump for MiSeq.")


def fqScreenRunner(cmd):
    cmds = cmd.split(" ")
    fqScreenRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = fqScreenRun.wait()
    return exitcode


def fastqscreen(project, laneFolder, sampleIDs, config):
    log.info("Fastq_screen for {}".format(project))
    screenRunnerCmds = []
    for ID in sampleIDs:
        IDfolder = os.path.join(
            laneFolder,
            "FASTQC_Project_" + project,
            "Sample_" + ID
        )
        if len(glob.glob(
            os.path.join(
                IDfolder,
                '*screen.txt'
            )
        )) == 0:
            sampleFolder = os.path.join(
                laneFolder,
                "Project_" + project,
                "Sample_" + ID
            )
            fqFile = screenFqFetcher(sampleFolder)
            screenRunnerCmds.append(
                config['software']['fastq_screen'] + " " +
                '-conf' + " " +
                os.path.join(
                    os.path.expanduser("~"),
                    'fastq_screen.conf'
                ) + " " +
                '--outdir' + " " +
                IDfolder + " " +
                '--subset' + " " +
                '1000000' + " " +
                '--quiet' + " " +
                '--threads' + " " +
                '4' + " " +
                fqFile
            )
    if screenRunnerCmds:
        with Pool(10) as p:
            screenReturns = p.map(fqScreenRunner, screenRunnerCmds)
            if screenReturns.count(0) == len(screenReturns):
                log.info("fastqScreen ran {}.".format(project))
            else:
                log.critical(
                    "fastqScreen failed for {}. exiting.".format(project)
                )
                mailHome(
                    laneFolder,
                    "fastqscreen runs failed. Investigate.",
                    config,
                    toCore=True
                )
                sys.exit(1)
    else:
        log.info(
            "fastqscreen files already present. Skipping {}".format(project)
        )


def multiqc(project, laneFolder, config, flowcell, sampleSheet):
    log.info("multiqc for {}".format(project))
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
        config['software']['multiqc'],
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
        log.info('multiqc ran for {}'.format(project))
        os.remove(confOut)
        os.remove(dataOut)
        os.remove(seqrepOut)
        os.remove(indexrepOut)
    else:
        log.critical("multiqc failed for {}".format(project))
        mailHome(
            laneFolder,
            "multiQC runs failed. Investigate.",
            config,
            toCore=True
        )
        sys.exit(1)
    return exitcode


def postmux(flowcell, sampleSheet, config):
    log.warning("Postmux module")
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
            log.info("Renaming {}".format(outLane))
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
            Path(
                os.path.join(laneFolder, 'renamed.done')
            ).touch()
        # postMux module
        if not os.path.exists(postmuxFlag):
            log.info("FastQC pool {}".format(outLane))
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
                # fastq_screen
                fastqscreen(
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
            log.info("Moving optical dup txt into FASTQC folder")
        moveOptDup(laneFolder)
        Path(
                os.path.join(laneFolder, 'postmux.done')
        ).touch()
    return (0)
