import os
import sys
import glob
from pathlib import Path
from rich import print


def fetchLatestSeqDir(pref, PI, postfix):
    globStr = os.path.join(
        pref,
        PI,
        postfix + '*'
    )
    if len(glob.glob(globStr)) == 1:
        return glob.glob(globStr)[0]
    else:
        maxFolder = 0
        for seqDir in glob.glob(
            os.path.join(
                pref,
                PI,
                postfix + '*'
            )
        ):
            try:
                seqInt = int(seqDir[-1])
            except ValueError:
                seqInt = 0
                continue
            if seqInt > maxFolder:
                maxFolder = seqInt
        return (os.path.join(
            pref,
            PI,
            postfix + str(maxFolder)
        ))


def fetchFolders(flowcellPath, piList, prefix, postfix):
    institute_PIs = piList
    flowcellPath = os.path.abspath(flowcellPath)
    FID = flowcellPath.split('/')[-1]
    projDic = {}
    try:
        int(FID[:6])
        print("[green]Valid flowcell folder.[/green]")
    except ValueError:
        sys.exit(
            "First 6 digits of flowcellpath don't convert to an int. Exiting."
        )
    for projF in glob.glob(
        os.path.join(
            flowcellPath,
            'Project_*'
        )
    ):
        proj = projF.split('/')[-1]
        PI = proj.split("_")[-1].lower()
        if PI == 'cabezas-wallscheid':
            PI = 'cabezas'
        if PI in institute_PIs:
            seqFolder = fetchLatestSeqDir(prefix, PI, postfix)
            if os.path.exists(
                os.path.join(seqFolder, FID)
            ):
                projDic[proj] = [
                    PI + 'grp',
                    [
                        os.path.join(seqFolder, FID),
                        os.path.join(seqFolder, FID, proj),
                        os.path.join(seqFolder, FID, 'FASTQC_' + proj),
                        os.path.join(
                            seqFolder,
                            FID,
                            'Analysis_' + proj.replace('Project_', '')
                        )
                    ]
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
        print("[bold red]wrong grp for {}! change this manually.[/bold red]!")
    os.chmod(flowcellF, 0o750)
    os.chmod(projectF, 0o750)
    os.chmod(fastqcF, 0o750)
    succes_project = release_rights(projectF)
    succes_fqc = release_rights(fastqcF)
    if os.path.exists(analysisF):
        succes_analysis = release_rights(analysisF)
<<<<<<< HEAD
        return (
            [succes_project, succes_fqc, succes_analysis]
        )
    return (
=======
        return(
            [succes_project, succes_fqc, succes_analysis]
        )
    return(
>>>>>>> f49305dbbf3abf78ca0cca5aa69ea9ce9cab8afb
        [succes_project, succes_fqc]
    )


def release_rights(F):
    changed = 0
    failed = 0
    for r, dirs, files in os.walk(F):
        for d in dirs:
            try:
                os.chmod(
                    os.path.join(r, d),
                    0o750
                )
                changed += 1
            except PermissionError:
                print("Permission error for {}".format(d))
                failed += 1
        for f in files:
            fil = os.path.join(r, f)
            if not os.path.islink(fil):
                try:
                    os.chmod(
                        fil,
                        0o750
                    )
                    changed += 1
                except PermissionError:
                    print("Permission error for {}".format(f))
                    failed += 1
    successRate = changed / (changed + failed)
    return (successRate)


def rel(flowcellPath, piList, prefix, postfix):
    projDic = fetchFolders(
        flowcellPath,
        piList,
        prefix,
        postfix
    )
    print("Print number of changed/(changed+unchanged)!")
    for proj in projDic:
        '''
        every projDic[proj] is a nested list of:
        [grp, [flowcell, project, fastqc]]
        '''
        successes = release_folder(projDic[proj][0], projDic[proj][1])
        if len(successes) == 2:
            print(
                "[green]Project[/green] {},{} proj,{} fqc".format(
                    proj,
                    successes[0],
                    successes[1]
                )
            )
        else:
            print(
                "[green]Project[/green] {},{} proj,{} fqc,{} analysis".format(
                    proj,
                    successes[0],
                    successes[1],
                    successes[2]
                )
            )
