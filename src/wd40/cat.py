import os
import subprocess
from rich import print
import glob
import sys
import shutil


def getSamples(p1, p2):
    samples_1 = glob.glob(
        os.path.join(
            p1,
            "Sample_*"
        )
    )
    samples_2 = glob.glob(
        os.path.join(
            p2,
            "Sample_*"
        )
    )
    samples_1 = sorted(os.path.basename(i) for i in samples_1)
    samples_2 = sorted(os.path.basename(i) for i in samples_2)

    if set(samples_1) != set(samples_2):
        print("Error")
    else:
        print(":thumbs_up: Moving on with {}".format(samples_1))
        return(samples_1)


def catRun(project, p1, p2, outDir):
    samples = getSamples(
        p1, p2
    )
    if not os.path.exists(os.path.join(outDir, project)):
        os.mkdir(os.path.join(outDir, project))
    for sample in samples:
        oDir = os.path.join(outDir, project, sample)
        if os.path.exists(oDir):
            shutil.rmtree(oDir)
        os.mkdir(oDir)
        for fqFile in glob.glob(
            os.path.join(
                p1,
                sample,
                '*fastq.gz'
            )
        ):
            fqFile = os.path.basename(fqFile)
            fq1 = os.path.join(
                p1,
                sample,
                fqFile
            )
            fq2 = os.path.join(
                p2,
                sample,
                fqFile
            )
            if not os.path.exists(fq1):
                sys.exit("{} not found.".format(fq1))
            if not os.path.exists(fq2):
                sys.exit("{} not found.".format(fq2))
            fqO = os.path.join(
                oDir,
                fqFile
            )
            print("writing to {}".format(fqO))
            with open(fqO, 'w') as f:
                subprocess.run(
                    ['cat', fq1, fq2],
                    stdout=f
                )
