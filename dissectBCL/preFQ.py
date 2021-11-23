from dissectBCL.fakeNews import log

def cycleOverwrite(readLens):
    log.info("parsing read lengths.")
    cycleList = []
    for read in readLens:
        if read[1] == 'Read':
            cycleList.append("Y"+str(read[0]))
        if read[1] == 'Index' and read[0] < 16:
            cycleList.append("I "+str(read[0]))
        elif read[1] == 'Index' and read[0] >= 16:
            cycleList.append("Y"+str(read[0]))
    return cycleList
        
def detMask(cycleList, sampleSheetDF):
    log.info("determine masking")
    maskStr = cycleList[0]
    minP7 = sampleSheetDF['index'].str.len().min()
    P7seqLen = int(cycleList[1].split(' ')[1])
    if P7seqLen > minP7:
        maskStr += ";I" + str(minP7) + "N" + str(P7seqLen - minP7)
    if 'I' in cycleList[2]:
        minP5 = sampleSheetDF['index2'].str.len().min()
        P5seqLen = int(cycleList[2].split(' ')[1])
        if P5seqLen > minP5:
            maskStr += ";I" + str(minP5) + "N" + str(P5seqLen - minP5)
    else:
        maskStr += ';' + cycleList[2]
    maskStr += ';' + cycleList[3]
    return maskStr
    



def prepConvert(flowcell, sampleSheet):
    log.warning("PreFQ module")
    cycleList = cycleOverwrite(flowcell.readLens)
    log.info("determine masking")
    for ss in sampleSheet.ssDic:
        sampleSheet.ssDic[ss]['maskStr'] = detMask(cycleList, sampleSheet.ssDic[ss]['sampleSheet'])
    return sampleSheet
