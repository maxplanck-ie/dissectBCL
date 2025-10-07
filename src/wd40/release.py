import os
import requests
from subprocess import check_output
import sys
import glob
from pathlib import Path
from rich import print

def fetchLatestSeqDir(pref, PI, postfix):
    globStr = os.path.join(pref, PI, postfix + "*")
    if len(glob.glob(globStr)) == 1:
        return glob.glob(globStr)[0]
    else:
        maxFolder = 0
        seqInt = 0
        for seqDir in glob.glob(os.path.join(pref, PI, postfix + "*")):
            seqDirStrip = seqDir.split('/')[-1].replace('sequencing_data', '')
            if seqDirStrip:
                seqInt = int(seqDirStrip)
            if seqInt > maxFolder:
                maxFolder = seqInt
        return os.path.join(pref, PI, postfix + str(maxFolder))


def fetchFolders(flowcellPath, piList, prefix, postfix, fexBool, parkourVars):
    parkourURL, parkourAuth, parkourCert, fromAddress = parkourVars
    institute_PIs = piList
    flowcellPath = os.path.abspath(flowcellPath)
    FID = flowcellPath.split("/")[-1]
    projDic = {}
    try:
        int(FID[:6])
        print("[green]Valid flowcell folder.[/green]")
    except ValueError:
        sys.exit(
            "First 6 digits of flowcellpath don't convert to an int. Exiting."
        )
    for projF in glob.glob(os.path.join(flowcellPath, "Project_*")):
        proj = projF.split("/")[-1]
        PI = proj.split("_")[-1].lower()
        if PI == "cabezas-wallscheid":
            PI = "cabezas"
        if PI in institute_PIs:
            seqFolder = fetchLatestSeqDir(prefix, PI, postfix)
            if os.path.exists(os.path.join(seqFolder, FID)):
                projDic[proj] = [
                    PI + "grp",
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
                    "[red]{} not found! Double check[/red]".format(
                        os.path.join(seqFolder, FID)
                    )
                )
        else:
            print(f"Assuming {proj} is fex'ed.")
            fexList = (
                check_output(["fexsend", "-l", fromAddress])
                .decode("utf-8")
                .replace("\n", " ")
                .split(" ")
            )

            tarBall = FID + "_" + proj + ".tar"
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
        print(
            "[bold red]wrong grp for {}! change manually![/bold red]!".format(
                grp
            )
        )
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
                print("Permission error for {}".format(d))
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
                    print("Permission error for {}".format(f))
                    failed += 1
                    failedfiles.append(f)
    successRate = changed / (changed + failed)
    if grouperror:
        print(f"[bold red]wrong grp (for some) {F}! change it![/bold red]!")
        print(f"files that have wrong group: {groupfiles}")
    return successRate


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
):
    projDic = fetchFolders(
        flowcellPath, piList, prefix, postfix, fexBool,
        (parkourURL, parkourAuth, parkourCert, fromAddress)
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
                "[green]Project[/green] {},{} proj,{} fqc".format(
                    proj, successes[0], successes[1]
                )
            )
        else:
            print(
                "[green]Project[/green] {},{} proj,{} fqc,{} analysis".format(
                    proj, successes[0], successes[1], successes[2]
                )
            )
        projectPath = projDic[proj][1][1].split("/")[-1]
        PI = (
            projectPath.split("_")[-1]
            .lower()
            .replace("cabezas-wallscheid", "cabezas")
        )
        d = None
        if PI in piList:
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
