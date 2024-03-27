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
        for seqDir in glob.glob(os.path.join(pref, PI, postfix + "*")):
            try:
                seqInt = int(seqDir[-1])
            except ValueError:
                seqInt = 0
                continue
            if seqInt > maxFolder:
                maxFolder = seqInt
        return os.path.join(pref, PI, postfix + str(maxFolder))


def fetchFolders(flowcellPath, piList, prefix, postfix):
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
            print("[bold cyan]Ignoring {}.[/bold cyan]".format(proj))
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
        succes_analysis = release_rights(analysisF, grp)
        return [succes_project, succes_fqc, succes_analysis]
    return [succes_project, succes_fqc]


def release_rights(F, grp):
    changed = 0
    failed = 0
    grouperror = False
    for r, dirs, files in os.walk(F):
        for d in dirs:
            try:
                os.chmod(os.path.join(r, d), 0o750)
                changed += 1
            except PermissionError:
                print("Permission error for {}".format(d))
                failed += 1
        for f in files:
            fil = os.path.join(r, f)
            if not os.path.islink(fil):
                if grp != Path(fil).group():
                    grouperror = True
                try:
                    os.chmod(fil, 0o750)
                    changed += 1
                except PermissionError:
                    print("Permission error for {}".format(f))
                    failed += 1
    successRate = changed / (changed + failed)
    if grouperror:
        print(
            "[bold red]wrong grp (for some) {}! change it![/bold red]!".format(
                F
            )
        )
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
    projDic = fetchFolders(flowcellPath, piList, prefix, postfix)
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
        elif fexBool:
            fexList = (
                check_output(["fexsend", "-l", fromAddress])
                .decode("utf-8")
                .replace("\n", " ")
                .split(" ")
            )
            tar_lane, tar_proj = projDic[proj][1][1].split("/")[-2:]
            tarBall = tar_lane + "_" + tar_proj + ".tar"
            if tarBall in fexList:
                d = {"data": tarBall, "metadata": None}
            else:
                print("fexLink: ", tarBall, " not found!")
        if d:
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
