from dissectBCL.fakeNews import mailHome
from dissectBCL.misc import krakenfqs
from dissectBCL.misc import multiQC_yaml
import hashlib
import logging
import os
from multiprocessing import Pool
from pandas import isna
import ruamel.yaml
import re
import shutil
from subprocess import Popen, DEVNULL
import sys
from pathlib import Path


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
    oldName = fqFile.name
    # 24L002006_S63_L001_R2_001.fastq.gz -> 24L002006
    sampleID = oldName.split('_')[0]
    # 24L002006 -> sample_name.txt
    sampleName = matchIDtoName(sampleID, ssdf)
    sampleIDPath = projectFolder / f"Sample_{sampleID}"
    sampleIDPath.mkdir(exist_ok = True)

    # Create new name
    if laneSplitStatus:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+L[0-9][0-9][0-9]_+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(regstr, r'_\1', newName)
        newName.replace(sampleID, sampleName)
    else:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(regstr, r'_\1', newName)
    logging.debug(f"Postmux - rename - suggesting {oldName} into {newName}")
    return sampleIDPath / newName


def renameProject(projectFolder, ssdf, laneSplitStatus):
    """
    rename and move files under sample_ID folders.
    rename project folder from e.g.
    1906_Hein_B03_Hein -> Project_1906_Hein_B03_Hein
    """

    logging.info(f"Postmux - Renaming {projectFolder}")
    for fq in projectFolder.glob('*fastq.gz'):
        newName = renamefq(fq, projectFolder, ssdf, laneSplitStatus)
        shutil.move(fq, newName)

    # Finally rename the project folder.
    # With Aviti data, the data lives under 'Samples' directory. We don't want to retain this.
    # Remove 'Samples' from the path parts
    parts = [p for p in projectFolder.parts if p != "Samples"]
    projectFolder_clean = Path(*parts)
    print(projectFolder)
    print(projectFolder_clean.with_stem("Project_" + projectFolder.stem))
    shutil.move(
        projectFolder,
        projectFolder_clean.with_stem("Project_" + projectFolder.stem)
    )


def validateFqEnds(pdir, flowcell):
    """
    recursively looks for fastq.gz files,
    validates the ending (e.g. R1, R2, I1, I2)
    ignores 'Undetermined'

    """
    malformat = []
    for f in pdir.rglob('*fastq.gz'):
        if 'Undetermined' not in f:
            e = f.name.split('.')[0]
            if e[-2:] not in [
                'R1', 'R2', 'I1', 'I2'
            ]:
                malformat.append(e)
    if not malformat:
        logging.info(f"Postmux - all fastq files in {pdir} have proper ending.")
    else:
        _msg = f"Improper fastq file format: {malformat}"
        logging.critical(_msg)
        mailHome(
            flowcell.name,
            _msg,
            flowcell.config
        )
        sys.exit(1)


def fqcRunner(cmd):
    cmds = cmd.split(" ")
    qcRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = qcRun.wait()
    return (exitcode)


def qcs(project, laneFolder, sampleIDs, config):
    # make fastqc folder.
    fqcFolder = laneFolder / f"FASTQC_Project_{project}"
    fqcFolder.mkdir(exist_ok=True)
    fastqcCmds = []
    # Decide threading setup - aim to have 2 threads per fastqc instance.
    num_pool_runners = max(1, int(config['misc']['threads']) // 2)
    for ID in sampleIDs:
        # Colliding samples are omitted, and don't have a folder.
        fqFolder = laneFolder / f"Project_{project}" / f"Sample_{ID}"
        if not fqFolder.exists():
            continue
        IDFolder = fqcFolder / f"Sample_{ID}"
        IDFolder.mkdir(exist_ok=True)
        fqFiles = [str(i) for i in fqFolder.glob("*fastq.gz")]
        # Don't do double work.
        if len(list(IDFolder.glob("*zip"))) == 0:
            fastqcCmds.append(
                " ".join([
                    'fastqc',
                    '-a',
                    config['software']['fastqc_adapters'],
                    '-q',
                    '-t',
                    "2",
                    '-o',
                    IDFolder._str
                ] + fqFiles)
            )
    if fastqcCmds:
        logging.info(f"Postmux - FastQC - command example: {project} - {fastqcCmds[0]}")
        with Pool(num_pool_runners) as p:
            fqcReturns = p.map(fqcRunner, fastqcCmds)
            if fqcReturns.count(0) == len(fqcReturns):
                logging.info(f"Postmux - FastQC done for {project}.")
            else:
                logging.critical(f"Postmux - FastQC crashed for {project}. exiting.")
                mailHome(
                    laneFolder,
                    f"FastQC runs failed for project {project}.",
                    config,
                    toCore=True
                )
                sys.exit(1)
    else:
        logging.info(f"Postmux - Seems all FastQCs already done for {project}")


def clmpRunner(cmd):
    cmds = cmd.split(" ")
    effthreads = cmds.pop(-1)
    baseName = cmds.pop(-1)
    PE = str(cmds.pop(-1))
    samplePath = cmds.pop(-1)
    os.chdir(samplePath)
    logging.info(f"Clumpify - {baseName}")
    clumpRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = clumpRun.wait()
    logging.info(f"Clumpify - {baseName} - splitfq")
    splitCmd = ['splitFastq', 'tmp.fq.gz', PE, baseName, effthreads ]
    splitFq = Popen(splitCmd, stdout=DEVNULL, stderr=DEVNULL)
    exitcode_split = splitFq.wait()
    os.remove('tmp.fq.gz')
    return (
        (exitcode, exitcode_split)
    )


def clumper(project, laneFolder, sampleIDs, config, PE, sequencer):
    # Decide threading setup - aim to have 2 threads per fastqc instance.
    configthreads = int(config['misc']['threads'])
    num_pool_runners = max(1, configthreads // 10)
    effthreads = (
        10 if configthreads >= 10
        else configthreads
    )
    clmpOpts = {
        'general': [
            'out=tmp.fq.gz',
            'dupesubs=0',
            'qin=33',
            'markduplicates=t',
            'optical=t',
            '-Xmx400G',
            f'threads={effthreads}',
            'tmpdir={}'.format(config['Dirs']['tempDir'])
        ],
        'NextSeq': [
            'spany=t',
            'adjacent=t',
            'dupedist=40'
        ],
        'NovaSeq': ['dupedist=12000'],
        'aviti': ['dupedist=12000'] # Take same for Aviti as for NovaSeq ?
    }
    clmpCmds = []
    if sequencer != 'MiSeq':
        for ID in sampleIDs:
            sampleDir = laneFolder / f"Project_{project}" / f"Sample_{ID}"
            if sampleDir.exists() and len(list(sampleDir.glob("*optical_duplicates*"))) == 0:
                fqFiles = list(sampleDir.glob("*fastq.gz"))
                if len(fqFiles) < 3:
                    if PE and len(fqFiles) == 2:
                        for i in fqFiles:
                            if '_R1.fastq.gz' in str(i):
                                in1 = "in=" + str(i)
                                baseName = i.name.replace('_R1.fastq.gz', '')
                            elif '_R2.fastq.gz' in str(i):
                                in2 = "in2=" + str(i)
                        clmpCmds.append(
                            'clumpify.sh' + " " +
                            in1 + " " +
                            in2 + " " +
                            " ".join(clmpOpts['general']) + " " +
                            " ".join(clmpOpts[sequencer]) + " " +
                            str(sampleDir) + " " +
                            "1" + " " +
                            baseName + " " + f"{effthreads}"
                        )
                    elif not PE and len(fqFiles) == 1:
                        if '_R1.fastq.gz' in str(fqFiles[0]):
                            in1 = "in=" + str(fqFiles[0])
                            baseName = fqFiles[0].name.replace('_R1.fastq.gz', '')
                            clmpCmds.append(
                                'clumpify.sh' + " " +
                                in1 + " " +
                                " ".join(clmpOpts['general']) + " " +
                                " ".join(clmpOpts[sequencer]) + " " +
                                sampleDir + " " +
                                "0" + " " +
                                baseName + " " + f"{effthreads}"
                            )
                        else:
                            logging.info(f"Not clumping {ID}")
        if clmpCmds:
            logging.info(f"Postmux - Clump - command example: {project} - {clmpCmds[0]}")
            with Pool(num_pool_runners) as p:
                clmpReturns = p.map(clmpRunner, clmpCmds)
                if clmpReturns.count((0, 0)) == len(clmpReturns):
                    logging.info(f"Postmux - Clumping done for {project}.")
                else:

                    logging.critical(f"Postmux - Clumping failed for {project}. Exiting.")
                    mailHome(
                        laneFolder,
                        f"Clump runs failed for {project}.",
                        config,
                        toCore=True
                    )
                    sys.exit(1)
        else:
            logging.info(f"Postmux - Clump - No clump run for {project}")
    else:
        logging.info("Postmux - Clump - no clumping for MiSeq.")


def krakRunner(cmd):
    cmds = cmd.split(" ")
    krakRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = krakRun.wait()
    return exitcode


def kraken(project, laneFolder, sampleIDs, config):
    configthreads = int(config['misc']['threads'])
    num_pool_runners = max(1, configthreads // 5)
    effthreads = (
        5 if configthreads >= 5
        else configthreads
    )
    krakenCmds = []
    for ID in sampleIDs:
        IDfolder = laneFolder / f"FASTQC_Project_{project}" / f"Sample_{ID}"
        if IDfolder.exists() and len(list(IDfolder.glob("*.rep"))) == 0:
            sampleFolder = laneFolder / f"Project_{project}" / f"Sample_{ID}"
            reportname, fqs = krakenfqs(sampleFolder)
            krakenCmds.append(
                ' '.join([
                    'kraken2',
                    '--db',
                    config['software']['kraken2db'],
                    '--out',
                    '-',
                    '--threads',
                    f'{effthreads}',
                    '--report',
                    reportname
                ] + fqs)
            )
    if krakenCmds:
        logging.info(f"Postmux - Kraken - command example: {project} - {krakenCmds[0]}")
        with Pool(num_pool_runners) as p:
            screenReturns = p.map(krakRunner, krakenCmds)
            if screenReturns.count(0) == len(screenReturns):
                logging.info(f"Postmux - Kraken done for {project}.")
            else:
                logging.critical(f"Postmux - Kraken failed for {project}. Exiting")
                mailHome(
                    laneFolder,
                    f"Kraken runs failed for {project}.",
                    config,
                    toCore=True
                )
                sys.exit(1)
    else:
        logging.info(f"Postmux - Kraken - No kraken run for {project}")


def md5Runner(fqfile):
    return (fqfile.name, hashlib.md5(open(fqfile, 'rb').read()).hexdigest())


def moveOptDup(laneFolder):
    for txt in laneFolder.glob('*/*/*duplicate.txt'):
        # Field -3 == project folder
        # escape those already in a fastqc folder (reruns)
        if 'FASTQC' not in str(txt):
            pathLis = str(txt).split('/')
            pathLis[-3] = 'FASTQC_' + pathLis[-3]
            ofile = "/".join(pathLis)
            os.rename(txt, ofile)


def md5_multiqc(project, laneFolder, flowcell):
    QCFolder = laneFolder / f"FASTQC_Project_{project}"
    projectFolder = laneFolder / f"Project_{project}"

    # md5sums
    logging.info(f"Postmux - md5sums - {project}")
    md5out = projectFolder / 'md5sums.txt'

    if not md5out.exists():
        with Pool(20) as p:
            _m5sums = p.map(md5Runner, list(projectFolder.glob("*/*fastq.gz")))
        with open(md5out, 'w') as f:
            for _m5sum in sorted(_m5sums, key=lambda x:x[0]):
                f.write(f"{_m5sum[0]}\t{_m5sum[1]}\n")

    # Always overwrite the multiQC reports. RunTimes are marginal anyway.
    mqcConf, mqcData, seqrepData, indexreportData = multiQC_yaml(
        flowcell,
        project,
        laneFolder
    )

    yaml = ruamel.yaml.YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    confOut = projectFolder / 'mqc.yaml'
    dataOut = QCFolder / 'parkour_mqc.tsv'
    seqrepOut = QCFolder / 'Sequencing_Report_mqc.tsv'
    indexrepOut = QCFolder / 'Index_Info_mqc.tsv'
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
        logging.info(f"Postmux - multiqc done for {project}")
        os.remove(confOut)
        os.remove(dataOut)
        os.remove(seqrepOut)
        os.remove(indexrepOut)
    else:
        logging.critical(f"Postmux - multiqc failed for {project}")
        mailHome(
            laneFolder,
            f"multiQC runs failed for {project}.",
            flowcell.config,
            toCore=True
        )
        sys.exit(1)