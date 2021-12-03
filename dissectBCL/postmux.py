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
        print(fqFiles)
        fastqcCmds.append(
            " ".join([
                config['software']['fastqc'],
                '-t',
                str(len(fqFiles)),
                '-o',
                IDfolder
            ] + fqFiles)
        )
    with Pool(20) as p:
        fqcReturns = p.map(fqcRunner, fastqcCmds)
        if fqcReturns.count(1) == len(fqcReturns):
            log.info("FastQC run for {}.".format(project))
        else:
            log.critical("FastQC runs failed. exiting.")
            sys.exit(1)


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
            # FastQC per project.
            for project in projects:
                qcs(
                    project,
                    laneFolder,
                    set(
                        df[df['Sample_Project'] == project]['Sample_ID']
                    ),
                    config
                )
