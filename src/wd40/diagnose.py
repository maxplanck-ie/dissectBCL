import json
from zipfile import ZipFile
import io
# import pandas as pd
import os
import numpy as np
from Bio.Seq import Seq
from rich import print, inspect
import glob
# import sys


def diagnose(fPath, solDir):
    print("[bold cyan]flowcell diagnosis[/bold cyan]")
    diagcls = barDiag(fPath)
    inspect(diagcls)


class barDiag:
    '''
    Tries to diagnose barcode issues.
    '''
    # @staticmethod
    # def fileNotExist(fileList):
    #     for f in fileList:
    #         if not os.path.exists(f):
    #             sys.exit(1)

    # @staticmethod
    # def grabreadCounts(ssdf):
    #     readCount = []
    #     for index, row in ssdf.iterrows():
    #         projID = "Project_" + row['Sample_Project']
    #         sampleID = "Sample_" + row['Sample_ID']
    #         sampleName = row['Sample_Name']
    #         zipStr = os.path.join("FASTQC_" + projID,
    #                             sampleID,
    #                             sampleName + "_R1_fastqc.zip")
    #         if os.path.exists(zipStr):
    #             readCount.append(barDiag.grabZipCount(zipStr))
    #         else:
    #             print("Warning, {} not found.".format(zipStr))
    #             readCount.append(0)
    #         ssdf['readCount'] = readCount
    #         return ssdf

    @staticmethod
    def grabZipCount(inputzip):
        dataStr = inputzip.split('/')[-1].split('.')[0] + "/fastqc_data.txt"
        with ZipFile(inputzip) as z:
            with z.open(dataStr) as f:
                qcData = io.TextIOWrapper(f, newline='\n')
                for line in qcData:
                    if line.startswith('Total'):
                        return int(line.strip().split()[2])

    @staticmethod
    def IDtoName(id, project, flowcellPath):
        fqcDir = os.path.join(
            flowcellPath,
            'FASTQC_' + project
        )
        idDir = os.path.join(
            fqcDir,
            'Sample_' + id
        )
        if os.path.exists(idDir):
            glob.glob("*_screen.txt")[0].replace()

    @staticmethod
    def parseSS(ssPath, convert=False):
        if convert:
            sampleDic = {}
            with open(ssPath) as f:
                for line in f:
                    if line.startswith('22L'):
                        lis = line.strip().split(',')
                        lis[2] = 'Project_' + lis[2]
                        if lis[0] not in sampleDic:
                            sampleDic[lis[0]] = lis[1:]
            print(sampleDic)
            return sampleDic
        return None

    def __init__(self, flowcellPath):
        _ssNames = ['SampleSheet.csv', 'demuxSheet.csv']
        for f in _ssNames:
            fPath = os.path.join(flowcellPath, f)
            if os.path.exists(fPath):
                if f == 'demuxSheet.csv':
                    self.ssPath = fPath
                    self.ssdf = barDiag.parseSS(f, convert=True)


def grabZipCount(inputzip):
    dataStr = inputzip.split('/')[-1].split('.')[0] + "/fastqc_data.txt"
    with ZipFile(inputzip) as z:
        with z.open(dataStr) as f:
            qcData = io.TextIOWrapper(f, newline='\n')
            for line in qcData:
                if line.startswith('Total'):
                    return int(line.strip().split()[2])


def revC(string):
    return str(Seq(string).reverse_complement())


def parseUnd(statFile, depth):
    with open(statFile) as f:
        stats = json.load(f)
        candidates = {}
        for lane in stats['UnknownBarcodes']:
            for comb in lane['Barcodes']:
                if np.round(
                    depth/lane['Barcodes'][comb]) == 1 or \
                   depth < lane['Barcodes'][comb]:
                    if comb not in candidates:
                        candidates[comb] = lane['Barcodes'][comb]
        return candidates


def crapMatcher(ssdf, pairedStatus, candidates, depth):
    # Fetch the combinations in candidates.
    # if pairedStatus == True:
    candNes = []
    for indexPair in candidates:
        candNes.append([
            indexPair.split('+')[0],
            indexPair.split('+')[1]
        ])
    samplesDic = {}
    for index, row in ssdf[ssdf['readCount'] < depth].iterrows():
        samplesDic[row['Sample_Name']] = [
            row['index'],
            row['index2']
        ]
    updateDic = {}
    for failure in samplesDic:
        for candidate in candNes:
            if samplesDic[
                failure
            ][0] in candidate and revC(
                samplesDic[failure][1]
            ) in candidate:
                updateDic[failure] = candidate
    updateDF = ssdf
    for update in updateDic:
        updateDF.loc[
            updateDF['Sample_Name'] == update, 'index'
        ] = updateDic[update][0]
        updateDF.loc[
            updateDF['Sample_Name'] == update, 'index2'
        ] = updateDic[update][1]
    del updateDF['readCount']
    return updateDF, len(updateDic)
