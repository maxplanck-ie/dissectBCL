import glob
import os
import sys
from pathlib import Path
from subprocess import check_output

import requests
from rich import print

from dissectBCL.misc import projectPI


def fetchLatestSeqDir(pref, PI, postfix):
    globStr = os.path.join(pref, PI, postfix + "*")
    if len(glob.glob(globStr)) == 1:
        return glob.glob(globStr)[0]
    else:
        maxFolder = 0
        seqInt = 0
        for seqDir in glob.glob(os.path.join(pref, PI, postfix + "*")):
            seqDirStrip = seqDir.split("/")[-1].replace("sequencing_data", "")
            if seqDirStrip:
                seqInt = int(seqDirStrip)
            if seqInt > maxFolder:
                maxFolder = seqInt
        return os.path.join(pref, PI, postfix + str(maxFolder))


def fetchFolders(
    flowcellPath, piList, prefix, postfix, fexBool, parkourVars, deliverTo=None
):
    parkourURL, parkourAuth, parkourCert, fromAddress = parkourVars
    # piList is the comma-joined PI string from config[Internals][PIs]; split it
    # so membership is an exact per-name match, not a substring test.
    institute_PIs = piList.split(",")
    # deliverTo maps a PI name to its IT filesystem dir token when it differs
    # from the name (parkour2 #317); membership still keys on the PI name.
    deliverTo = deliverTo or {}
    flowcellPath = os.path.abspath(flowcellPath)
    FID = flowcellPath.split("/")[-1]
    projDic = {}
    try:
        int(FID[:6])
        print("[green]Valid flowcell folder.[/green]")
    except ValueError:
        sys.exit("First 6 digits of flowcellpath don't convert to an int. Exiting.")
    for projF in glob.glob(os.path.join(flowcellPath, "Project_*")):
        proj = projF.split("/")[-1]
        PI = projectPI(proj)
        # A project is internal either when Parkour currently lists PI, or
        # when PI is a known deliver_to key - the latter covers a PI whose
        # Parkour name changed after some of their projects were already
        # created, so older projects still carry the PI's previous name.
        if PI in institute_PIs or PI in deliverTo:
            deliverDir = deliverTo.get(PI, PI)
            seqFolder = fetchLatestSeqDir(prefix, deliverDir, postfix)
            if os.path.exists(os.path.join(seqFolder, FID)):
                projDic[proj] = [
                    deliverDir + "grp",
                    [
                        os.path.join(seqFolder, FID),
                        os.path.join(seqFolder, FID, proj),
                        os.path.join(seqFolder, FID, "FASTQC_" + proj),
                        os.path.join(
                            seqFolder,
                            FID,
                            "Analysis_" + proj.replace("Project_", ""),
                        ),
                    ],
                ]
            else:
                print(
                    f"[red]{os.path.join(seqFolder, FID)} not found! Double check[/red]"
                )
        else:
            print(f"Assuming {proj} is fex'ed.")
            fexList = (
                check_output(["fexsend", "-l", fromAddress])
                .decode("utf-8")
                .replace("\n", " ")
                .split(" ")
            )

            tarBall = FID + "_" + proj + "_ro_crate.zip"
            if tarBall in fexList:
                if fexBool:
                    d = {"data": tarBall, "metadata": None}
                    print(
                        f"{tarBall} found in fexlist. Added filepaths to Parkour2: ",
                        requests.post(
                            parkourURL
                            + "/api/requests/"
                            + proj.split("_")[1]
                            + "/put_filepaths/",
                            auth=parkourAuth,
                            data=d,
                            verify=parkourCert,
                        ),
                    )
            else:
                print(f"{tarBall} not found in fex. Parkour not updated.")
    return projDic


def release_folder(grp, lis):
    flowcellF = lis[0]
    projectF = lis[1]
    fastqcF = lis[2]
    analysisF = lis[3]
    # flowcellF
    gotgrp = Path(flowcellF).group()
    if grp != gotgrp:
        print(f"[bold red]wrong grp for {grp}! change manually![/bold red]!")
    os.chmod(flowcellF, 0o750)
    os.chmod(projectF, 0o750)
    os.chmod(fastqcF, 0o750)
    succes_project = release_rights(projectF, grp)
    succes_fqc = release_rights(fastqcF, grp)
    if os.path.exists(analysisF):
        os.chmod(analysisF, 0o750)
        succes_analysis = release_rights(analysisF, grp)
        return [succes_project, succes_fqc, succes_analysis]
    return [succes_project, succes_fqc]


def release_rights(F, grp):
    changed = 0
    failed = 0
    failedfiles = []
    faileddirs = []
    grouperror = False
    groupfiles = []
    for r, dirs, files in os.walk(F):
        for d in dirs:
            try:
                os.chmod(os.path.join(r, d), 0o750)
                changed += 1
            except PermissionError:
                print(f"Permission error for {d}")
                failed += 1
                faileddirs.append(d)
        for f in files:
            fil = os.path.join(r, f)
            if not os.path.islink(fil):
                if grp != Path(fil).group():
                    grouperror = True
                    groupfiles.append(fil)
                try:
                    os.chmod(fil, 0o750)
                    changed += 1
                except PermissionError:
                    print(f"Permission error for {f}")
                    failed += 1
                    failedfiles.append(f)
    successRate = changed / (changed + failed)
    if grouperror:
        print(f"[bold red]wrong grp (for some) {F}! change it![/bold red]!")
        print(f"files that have wrong group: {groupfiles}")
    return successRate


def checkBRBDone(flowcellPath):
    """
    BigRedButton (BRB) touches an 'analysis.done' flag in the flowcell folder
    once it has finished writing to the Analysis_* folders. If that flag is
    missing, BRB may still be running (or hasn't picked up the flowcell yet)
    and releasing now can chmod/chown files out from under it, leaving the
    periphery copy with inconsistent permissions.
    """
    flowcellPath = os.path.abspath(flowcellPath)
    doneFlag = os.path.join(flowcellPath, "analysis.done")
    if not os.path.exists(doneFlag):
        print(
            "[bold red]Warning:[/bold red] no 'analysis.done' flag found in "
            f"{flowcellPath}. BigRedButton may still be running for this "
            "flowcell - releasing now can leave wrong permissions on files "
            "BRB writes afterwards. Double check before proceeding."
        )


def rel(
    flowcellPath,
    piList,
    prefix,
    postfix,
    parkourURL,
    parkourAuth,
    parkourCert,
    fexBool,
    fromAddress,
    deliverTo=None,
):
    # PI membership/lookup for the put_filepaths step below shares this map
    # with fetchFolders (parkour2 #317/#294) - keep it in lockstep, don't
    # hardcode individual PI overrides here.
    deliverTo = deliverTo or {}
    checkBRBDone(flowcellPath)
    projDic = fetchFolders(
        flowcellPath,
        piList,
        prefix,
        postfix,
        fexBool,
        (parkourURL, parkourAuth, parkourCert, fromAddress),
        deliverTo,
    )
    print("Print number of changed/(changed+unchanged)!")
    for proj in projDic:
        """
        every projDic[proj] is a nested list of:
        [grp, [flowcell, project, fastqc]]
        """
        successes = release_folder(projDic[proj][0], projDic[proj][1])
        if len(successes) == 2:
            print(
                f"[green]Project[/green] {proj},{successes[0]} proj,{successes[1]} fqc"
            )
        else:
            print(
                f"[green]Project[/green] {proj},{successes[0]} proj,{successes[1]} fqc,{successes[2]} analysis"
            )
        projectPath = projDic[proj][1][1].split("/")[-1]
        PI = projectPath.split("_")[-1].lower()
        # institute_PIs check mirrors fetchFolders(): PI is internal either
        # by matching Parkour's current PI list, or by being a known
        # deliver_to key (a PI whose Parkour name changed after some of
        # their projects were already created). Split piList for an exact
        # per-name match - `PI in piList` was comparing against the raw
        # comma-joined string, a substring match that could false-positive.
        institute_PIs = piList.split(",")
        isInternal = PI in institute_PIs or PI in deliverTo
        PI = deliverTo.get(PI, PI)
        d = None
        if isInternal:
            d = {
                "data": projDic[proj][1][1],
                "metadata": projDic[proj][1][1] + "/multiqc_report.html",
            }
            print(
                "Adding filepaths to Parkour2:",
                requests.post(
                    parkourURL
                    + "/api/requests/"
                    + proj.split("_")[1]
                    + "/put_filepaths/",
                    auth=parkourAuth,
                    data=d,
                    verify=parkourCert,
                ),
            )  # print the returned answer from the API
        else:
            print(f"{PI} not in piList for {proj}.")
