import hashlib
import logging
import os
import re
import shutil
import sys
from multiprocessing import Pool
from pathlib import Path
from subprocess import DEVNULL, Popen

import ruamel.yaml
from pandas import isna

from dissectBCL import screening
from dissectBCL.fakeNews import mailHome
from dissectBCL.misc import krakenfqs, multiQC_yaml


def matchIDtoName(ID, ssdf):
    if ID not in set(ssdf["Sample_ID"]):
        # can happen if filename is not legit
        # e.g. if demuxSheet did not match SampleSheet
        logging.critical(f"ID {ID} is not defined in SampleSheet.")
        sys.exit(1)

    name = ssdf[ssdf["Sample_ID"] == ID]["Sample_Name"].values

    if len(name) > 1:
        # It can happen one sample sits in 2 lanes.
        if len(set(name)) > 1:
            logging.critical(f"SampleID {ID} has multiple names {name}, exiting.")
            sys.exit(1)

        # can happen if ID is not listed in SampleSheet --> no Sample_Name
        elif isna(name[0]):
            logging.critical(f"Sample_Name is not defined for ID {ID} .")
            sys.exit(1)
        else:
            return name[0]
    else:
        return name[0]


def renamefq(fqFile, projectFolder, ssdf, laneSplitStatus):
    oldName = fqFile.name
    # 24L002006_S63_L001_R2_001.fastq.gz -> 24L002006
    sampleID = oldName.split("_")[0]
    # 24L002006 -> sample_name.txt
    sampleName = matchIDtoName(sampleID, ssdf)
    sampleIDPath = projectFolder / f"Sample_{sampleID}"
    sampleIDPath.mkdir(exist_ok=True)

    try:
        # check if any _stats.json exist, if yes, move to the Sample_XXXX directory
        if any(projectFolder.glob("*_stats.json")):
            sampleID_json = projectFolder / f"{sampleID}_stats.json"
            shutil.move(sampleID_json, sampleIDPath)
            print(f"File moved. Sample path: {sampleIDPath}")
        else:
            print("No .json file found in projectFolder. Skipping move.")
    except Exception as e:
        print(f"Unexpected error: {e}")

    # Create new name
    if laneSplitStatus:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+L[0-9][0-9][0-9]_+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(regstr, r"_\1", newName)
        newName.replace(sampleID, sampleName)
    else:
        newName = oldName.replace(sampleID, sampleName)
        regstr = r"_S[0-9]?[0-9]?[0-9]?[0-9]_"
        regstr += r"+([IR][123])+_[0-9][0-9][0-9]"
        newName = re.sub(regstr, r"_\1", newName)
    logging.debug(f"Postmux - rename - suggesting {oldName} into {newName}")
    return sampleIDPath / newName


def renameProject(projectFolder, ssdf, laneSplitStatus):
    """
    rename and move files under sample_ID folders.
    rename project folder from e.g.
    1906_Hein_B03_Hein -> Project_1906_Hein_B03_Hein
    """

    logging.info(f"Postmux - Renaming {projectFolder}")
    for fq in projectFolder.glob("*fastq.gz"):
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
        projectFolder, projectFolder_clean.with_stem("Project_" + projectFolder.stem)
    )


def validateFqEnds(pdir, flowcell):
    """
    recursively looks for fastq.gz files,
    validates the ending (e.g. R1, R2, I1, I2)
    ignores 'Undetermined'

    """
    malformat = []
    for f in pdir.rglob("*fastq.gz"):
        if "Undetermined" not in f:
            e = f.name.split(".")[0]
            if e[-2:] not in ["R1", "R2", "I1", "I2"]:
                malformat.append(e)
    if not malformat:
        logging.info(f"Postmux - all fastq files in {pdir} have proper ending.")
    else:
        _msg = f"Improper fastq file format: {malformat}"
        logging.critical(_msg)
        mailHome(flowcell.name, _msg, flowcell.config)
        sys.exit(1)


def fqcRunner(cmd):
    cmds = cmd.split(" ")
    qcRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = qcRun.wait()
    return exitcode


def qcs(project, laneFolder, sampleIDs, config):
    # make fastqc folder.
    fqcFolder = laneFolder / f"FASTQC_Project_{project}"
    fqcFolder.mkdir(exist_ok=True)
    fastqcCmds = []
    # Decide threading setup - aim to have 2 threads per fastqc instance.
    num_pool_runners = max(1, int(config["misc"]["threads"]) // 2)
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
                " ".join(
                    [
                        "fastqc",
                        "-a",
                        config["software"]["fastqc_adapters"],
                        "-q",
                        "-t",
                        "2",
                        "-o",
                        IDFolder._str,
                    ]
                    + fqFiles
                )
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
                    toCore=True,
                )
                sys.exit(1)
    else:
        logging.info(f"Postmux - Seems all FastQCs already done for {project}")


def clmpRunner(cmd):
    cmds = cmd.split(" ")
    splitFastqBin = cmds.pop(-1)
    effthreads = cmds.pop(-1)
    baseName = cmds.pop(-1)
    PE = str(cmds.pop(-1))
    samplePath = cmds.pop(-1)
    os.chdir(samplePath)
    logging.info(f"Clumpify - {baseName}")
    clumpRun = Popen(cmds, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = clumpRun.wait()
    if exitcode != 0 or not os.path.exists("tmp.fq.gz"):
        logging.critical(
            f"Clumpify - {baseName} - clumpify failed (exit {exitcode}), "
            "no tmp.fq.gz produced."
        )
        return (exitcode if exitcode != 0 else 1, 1)
    logging.info(f"Clumpify - {baseName} - splitfq")
    splitCmd = [splitFastqBin]
    if PE == "0":
        splitCmd.append("--SE")
    splitCmd += ["--pigzThreads", str(effthreads), "tmp.fq.gz", baseName]
    splitFq = Popen(splitCmd, stdout=DEVNULL, stderr=DEVNULL)
    exitcode_split = splitFq.wait()
    if os.path.exists("tmp.fq.gz"):
        os.remove("tmp.fq.gz")
    return (exitcode, exitcode_split)


def clumper(project, laneFolder, sampleIDs, config, PE, sequencer):
    # Decide threading setup - aim to have 2 threads per fastqc instance.
    configthreads = int(config["misc"]["threads"])
    num_pool_runners = max(1, configthreads // 10)
    effthreads = 10 if configthreads >= 10 else configthreads
    clmpOpts = {
        "general": [
            "out=tmp.fq.gz",
            "dupesubs=0",
            "qin=33",
            "markduplicates=t",
            "optical=t",
            "-Xmx650G",
            f"threads={effthreads}",
            "tmpdir={}".format(config["Dirs"]["tempDir"]),
        ],
        "NextSeq": ["spany=t", "adjacent=t", "dupedist=40"],
        "NovaSeq": ["dupedist=12000"],
    }
    clmpOpts["aviti"] = clmpOpts["NextSeq"].copy()

    clmpCmds = []
    if sequencer != "MiSeq":
        for ID in sampleIDs:
            sampleDir = laneFolder / f"Project_{project}" / f"Sample_{ID}"
            if (
                sampleDir.exists()
                and len(list(sampleDir.glob("*optical_duplicates*"))) == 0
            ):
                fqFiles = list(sampleDir.glob("*fastq.gz"))
                if len(fqFiles) < 3:
                    if PE and len(fqFiles) == 2:
                        for i in fqFiles:
                            if "_R1.fastq.gz" in str(i):
                                in1 = "in=" + str(i)
                                baseName = i.name.replace("_R1.fastq.gz", "")
                            elif "_R2.fastq.gz" in str(i):
                                in2 = "in2=" + str(i)
                        clmpCmds.append(
                            "clumpify.sh"
                            + " "
                            + in1
                            + " "
                            + in2
                            + " "
                            + " ".join(clmpOpts["general"])
                            + " "
                            + " ".join(clmpOpts[sequencer])
                            + " "
                            + str(sampleDir)
                            + " "
                            + "1"
                            + " "
                            + baseName
                            + " "
                            + f"{effthreads}"
                            + " "
                            + config["software"]["splitFastq"]
                        )
                    elif not PE and len(fqFiles) == 1:
                        if "_R1.fastq.gz" in str(fqFiles[0]):
                            in1 = "in=" + str(fqFiles[0])
                            baseName = fqFiles[0].name.replace("_R1.fastq.gz", "")
                            clmpCmds.append(
                                "clumpify.sh"
                                + " "
                                + in1
                                + " "
                                + " ".join(clmpOpts["general"])
                                + " "
                                + " ".join(clmpOpts[sequencer])
                                + " "
                                + str(sampleDir)
                                + " "
                                + "0"
                                + " "
                                + baseName
                                + " "
                                + f"{effthreads}"
                                + " "
                                + config["software"]["splitFastq"]
                            )
                        else:
                            logging.info(f"Not clumping {ID}")
        if clmpCmds:
            logging.info(
                f"Postmux - Clump - command example: {project} - {clmpCmds[0]}"
            )
            with Pool(num_pool_runners) as p:
                clmpReturns = p.map(clmpRunner, clmpCmds)
                if clmpReturns.count((0, 0)) == len(clmpReturns):
                    logging.info(f"Postmux - Clumping done for {project}.")
                else:
                    logging.critical(
                        f"Postmux - Clumping failed for {project}. Exiting."
                    )
                    mailHome(
                        laneFolder,
                        f"Clump runs failed for {project}.",
                        config,
                        toCore=True,
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


def kraken(project, laneFolder, sampleIDs, ssdf, config):
    configthreads = int(config["misc"]["threads"])
    num_pool_runners = max(1, configthreads // 5)
    effthreads = 5 if configthreads >= 5 else configthreads
    krakenCmds = []
    for ID in sampleIDs:
        IDfolder = laneFolder / f"FASTQC_Project_{project}" / f"Sample_{ID}"
        if IDfolder.exists() and len(list(IDfolder.glob("*.rep"))) == 0:
            sampleFolder = laneFolder / f"Project_{project}" / f"Sample_{ID}"
            reportname, fqs = krakenfqs(sampleFolder)
            krakenCmds.append(
                " ".join(
                    [
                        "kraken2",
                        "--db",
                        config["software"]["kraken2db"],
                        "--out",
                        "-",
                        "--threads",
                        f"{effthreads}",
                        "--report",
                        reportname,
                    ]
                    + fqs
                )
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
                    toCore=True,
                )
                sys.exit(1)
    else:
        logging.info(f"Postmux - Kraken - No kraken run for {project}")

    # PlusPF escalation: re-screen any sample whose unclassified fraction
    # exceeds its Library_Type's threshold against the broader PlusPF db.
    # Deployed configs that predate this feature won't have [screening] --
    # degrade to a no-op rather than crash the flowcell.
    if not config.has_section("screening"):
        return
    plusPFdb = config["screening"].get("plusPFdb", fallback="")
    if not plusPFdb or not Path(plusPFdb).exists():
        logging.info("Postmux - PlusPF escalation skipped: plusPFdb not configured.")
        return
    if "Library_Type" not in ssdf.columns:
        return
    escalateIDs = []
    for ID in sampleIDs:
        IDfolder = laneFolder / f"FASTQC_Project_{project}" / f"Sample_{ID}"
        sampleFolder = laneFolder / f"Project_{project}" / f"Sample_{ID}"
        if not sampleFolder.exists() or not IDfolder.exists():
            continue
        if list(IDfolder.glob("*.plusPF.krakenreport")):
            continue  # already escalated in a prior run
        try:
            fqInfo = krakenfqs(sampleFolder)
        except IndexError:
            # krakenfqs() indexes into an empty fastq list when a sample
            # folder has zero matching fastq files -- treat the same as
            # its "no usable fastqs" None return, below.
            fqInfo = None
        if not fqInfo:
            continue
        reportname, _ = fqInfo
        if not Path(reportname).exists():
            continue  # kraken2 hasn't produced a report for this sample yet
        libraryTypes = ssdf[ssdf["Sample_ID"] == ID]["Library_Type"].values
        libraryType = libraryTypes[0] if len(libraryTypes) else None
        if screening.needsEscalation(Path(reportname), libraryType, config):
            escalateIDs.append(ID)
    if escalateIDs:
        logging.info(
            f"Postmux - Kraken - PlusPF escalation flagged for {project}: {escalateIDs}"
        )
        runPlusPF(project, laneFolder, escalateIDs, config)


def runPlusPF(project, laneFolder, sampleIDs, config):
    """
    Re-screens sampleIDs (already flagged by screening.needsEscalation)
    against the broader PlusPF kraken2 database, writing '<sample>.plusPF.krakenreport'
    next to the routine '<sample>.rep'. Unlike kraken(), a failed run here
    does not abort the flowcell -- PlusPF is a supplementary check on
    already-demuxed, already-shippable data. Any report left behind by a
    failed run is removed, so a later run doesn't mistake a partial report
    for a completed escalation (kraken2 can write a --report file before
    later failing).
    """
    configthreads = int(config["misc"]["threads"])
    # Unlike kraken()'s small contaminomedb, the PlusPF index is ~75-80GB
    # and kraken2 loads the whole hash into a private per-process heap
    # without --memory-mapping. Escalation only flags ~1-2 samples a week
    # in practice, so there's no throughput pressure to run many of these
    # concurrently -- half of kraken()'s configthreads // 5 pooling
    # constant keeps a wide margin against several ~80GB indices loading
    # into RAM at once, without forcing every escalation onto one worker.
    # --memory-mapping lets repeat/concurrent runs share the index via
    # page cache instead of each reloading it from scratch.
    num_pool_runners = max(1, configthreads // 10)
    effthreads = 5 if configthreads >= 5 else configthreads
    krakenCmds = []
    reportPaths = []
    for ID in sampleIDs:
        sampleFolder = laneFolder / f"Project_{project}" / f"Sample_{ID}"
        reportname, fqs = krakenfqs(sampleFolder)
        # reportname always ends in ".rep" (see krakenfqs) -- slice off
        # just that suffix rather than a global .replace(), which could
        # also rewrite an unrelated ".rep" earlier in the path.
        plusReportname = reportname[: -len(".rep")] + ".plusPF.krakenreport"
        reportPaths.append(plusReportname)
        krakenCmds.append(
            " ".join(
                [
                    "kraken2",
                    "--db",
                    config["screening"]["plusPFdb"],
                    "--out",
                    "-",
                    "--threads",
                    f"{effthreads}",
                    "--memory-mapping",
                    "--report",
                    plusReportname,
                ]
                + fqs
            )
        )
    if krakenCmds:
        logging.info(
            f"Postmux - PlusPF escalation - command example: {project} - {krakenCmds[0]}"
        )
        with Pool(num_pool_runners) as p:
            screenReturns = p.map(krakRunner, krakenCmds)
        if screenReturns.count(0) == len(screenReturns):
            logging.info(f"Postmux - PlusPF escalation done for {project}.")
        else:
            logging.critical(f"Postmux - PlusPF escalation failed for {project}.")
            for returncode, reportPath in zip(screenReturns, reportPaths, strict=True):
                if returncode != 0:
                    Path(reportPath).unlink(missing_ok=True)
            mailHome(
                laneFolder,
                f"PlusPF escalation runs failed for {project}.",
                config,
                toCore=True,
            )
    else:
        logging.info(f"Postmux - PlusPF escalation - no samples flagged for {project}")


def md5Runner(fqfile):
    md5 = hashlib.md5()
    with open(fqfile, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5.update(chunk)
    return (fqfile.name, md5.hexdigest())


def moveOptDup(laneFolder):
    for txt in laneFolder.glob("*/*/*.metrics"):
        # Field -3 == project folder
        # escape those already in a fastqc folder (reruns)
        if "FASTQC" not in str(txt):
            pathLis = str(txt).split("/")
            pathLis[-3] = "FASTQC_" + pathLis[-3]
            ofile = "/".join(pathLis)
            os.rename(txt, ofile)


def md5_multiqc(project, laneFolder, flowcell):
    QCFolder = laneFolder / f"FASTQC_Project_{project}"
    projectFolder = laneFolder / f"Project_{project}"

    # md5sums
    logging.info(f"Postmux - md5sums - {project}")
    md5out = projectFolder / "md5sums.txt"

    if not md5out.exists():
        with Pool(20) as p:
            _m5sums = p.map(md5Runner, list(projectFolder.glob("*/*fastq.gz")))
        with open(md5out, "w") as f:
            for _m5sum in sorted(_m5sums, key=lambda x: x[0]):
                f.write(f"{_m5sum[0]}\t{_m5sum[1]}\n")

    # Always overwrite the multiQC reports. RunTimes are marginal anyway.
    mqcConf, mqcData, seqrepData, indexreportData, plusPFData = multiQC_yaml(
        flowcell, project, laneFolder
    )

    yaml = ruamel.yaml.YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    confOut = projectFolder / "mqc.yaml"
    dataOut = QCFolder / "parkour_mqc.tsv"
    seqrepOut = QCFolder / "Sequencing_Report_mqc.tsv"
    indexrepOut = QCFolder / "Index_Info_mqc.tsv"
    plusPFOut = QCFolder / "PlusPF_Escalation_mqc.tsv"
    with open(confOut, "w") as f:
        yaml.dump(mqcConf, f)
    with open(seqrepOut, "w") as f:
        f.write(seqrepData)
    with open(dataOut, "w") as f:
        f.write(mqcData)
    with open(indexrepOut, "w") as f:
        f.write(indexreportData)
    # Only write (and later remove) this one when samples were actually
    # escalated -- an always-present, always-empty section is just noise.
    if plusPFData:
        with open(plusPFOut, "w") as f:
            f.write(plusPFData)
    multiqcCmd = [
        "multiqc",
        "--quiet",
        "--no-data-dir",
        "-f",
        "-o",
        projectFolder,
        "-c",
        confOut,
        QCFolder,
    ]
    multiqcRun = Popen(multiqcCmd, stdout=DEVNULL, stderr=DEVNULL)
    exitcode = multiqcRun.wait()
    if exitcode == 0:
        logging.info(f"Postmux - multiqc done for {project}")
        os.remove(confOut)
        os.remove(dataOut)
        os.remove(seqrepOut)
        os.remove(indexrepOut)
        if plusPFData:
            os.remove(plusPFOut)
    else:
        logging.critical(f"Postmux - multiqc failed for {project}")
        mailHome(
            laneFolder,
            f"multiQC runs failed for {project}.",
            flowcell.config,
            toCore=True,
        )
        sys.exit(1)
