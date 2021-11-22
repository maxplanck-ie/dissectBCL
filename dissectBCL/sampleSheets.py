import pandas as pd

def ixMask(ss, indexLen, index2Len=None, readLen1, readLen2):
    "Returns the masking parameter based on Index Length and actual sequence as listed in the sampleSheet."
    shortestIx = min( ss['index'].str.len() )
    if shortestIx < indexLen:
        ix1 = "I" + str(shortestIx) + "N" + str(indexLen-shortestIx)
        if not index2Len:
            return "OverrideCycles,Y" + str(readLen1)
    shortestIx2 = min(ss['index2'].str.len())
    if shortestIx2 < indexLen:
        ix2 = "I" + str(shortestIx) + "N" + str(indexLen-shortestIx)
        return ix1 + ";"
