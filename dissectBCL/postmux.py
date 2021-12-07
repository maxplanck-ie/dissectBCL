import os
import glob
from dissectBCL.fakeNews import log
from pathlib import Path
import re
import shutil
import sys
from multiprocessing import Pool
from subprocess import Popen


def matchIDtoName(ID, ssdf):
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
        else:
            return name[0]
    else:
        return name[0]


def renamefq(fqFile, projectFolder, ssdf):
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
    newName = oldName.replace(sampleID, sampleName)
    newName = re.sub(
        r"_S[0-9]?[0-9]_+L[0-9][0-9][0-9]_+([IR][123])+_[0-9][0-9][0-9]",
        r'_\1',
        newName
    )
    newName.replace(sampleID, sampleName)
    return os.path.join(sampleIDPath, newName)


def renameProject(projectFolder, ssdf):
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
        newName = renamefq(fq, projectFolder, ssdf)
        shutil.move(fq, newName)
    # Finally rename the project folder.
    projID = projectFolder.split('/')[-1]
    shutil.move(
        projectFolder,
        projectFolder.replace(projID, "Project_" + projID)
    )


def fqcRunner(cmd):
    cmds = cmd.split(" ")
    qcRun = Popen(cmds)
    exitcode = qcRun.wait()
    return(exitcode)


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
                log.critical("FastQC runs failed for {}. exiting.".format(project))
                sys.exit(1)
    else:
        log.info("FastQCs already done for {}".format(project))


def clmpRunner(cmd):
    cmds = cmd.split(" ")
    samplePath = cmds.pop(-1)
    os.chdir(samplePath)
    clumpRun = Popen(cmds)
    exitcode = clumpRun.wait()
    return(exitcode)

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
    for ID in sampleIDs:
        sampleDir = os.path.join(
            laneFolder,
            "Project_" + project,
            "Sample_" + ID
        )
        fqFiles = glob.glob(
            os.path.join(
                sampleDir,
                "*fastq.gz"
            )
        )
        if len(fqFiles) < 3:
            if PE and len(fqFiles) == 2:
                for i in fqFiles:
                    if 'R1' in i:
                        in1 = "in=" + i
                    elif 'R2' in i:
                        in2 = "in2=" + i
                clmpCmds.append(
                    config['software']['clumpify'] + " " +\
                    in1 + " " +\
                    in2 + " " +\
                    " ".join(clmpOpts['general']) + " " +\
                    " ".join(clmpOpts[sequencer]) + " " +\
                    sampleDir
                )
            elif not PE and len(fqFiles) == 1:
                if 'R1' in fqFiles[0]:
                    in1 = "in=" + fqFiles[0]
                    clmpCmds.append(
                        config['software']['clumpify'] + " " +\
                        in1 + " " +\
                        " ".join(clpmOpts['general']) + " " +\
                        " ".join(clmpOpts[sequencer]) + " " +\
                        sampleDir
                    )
                else:
                    log.info("Not clumping {}".format(ID))
    if clmpCmds:
        with Pool(5) as p:
            clmpReturns = p.map(clmpRunner, clmpCmds)
            if clmpReturns.count(0) == len(clmpReturns):
                log.info("Clumping done for {}.".format(project))
            else:
                log.critical("Clumping failed for {}. exiting.".format(project))
                sys.exit(1)
    else:
        log.info("No clump run for {}".format(project))


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
                    df
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

