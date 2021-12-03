import os
import glob
from dissectBCL.fakeNews import log
import re
import shutil
import sys


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

def qcs(project)
def postmux(flowcell, sampleSheet, config):
    log.warning("Postmux module")
    for outLane in sampleSheet.ssDic:
        log.info("Processing {}".format(outLane))
        df = sampleSheet.ssDic[outLane]['sampleSheet']
        for project in df['Sample_Project'].unique():
            renameProject(
                os.path.join(
                    flowcell.outBaseDir,
                    outLane,
                    project
                ),
                sampleSheet.ssDic[outLane]['sampleSheet']
            )