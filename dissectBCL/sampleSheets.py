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


# Parse sampleSheet
def parseSS(ssPath):
    """
    We read the sampleSheet csv, and remove the stuff above the header.
    """
    ssdf = pd.read_csv(ssPath, sep=',')
    # There is a bunch of header 'junk' that we don't want.
    # subset the df from [data] onwards.
    startIx = ssdf[ssdf.iloc[:, 0] == '[Data]'].index.values[0] + 1
    # only take those with actual sample information.
    ssdf = ssdf.iloc[startIx:, :]
    ssdf.columns = ssdf.iloc[0]
    ssdf = ssdf.drop(ssdf.index[0])
    # Reset index
    ssdf.reset_index(inplace=True)
    ssdf.index.name = None
    # Remove 'level0' column
    ssdf.drop('level_0', axis=1, inplace=True)
    ssdf = ssdf.dropna(axis=1, how='all')
    return ssdf

def decideSplit(ss, lanes):
    """
    Do we need to split per lane ? 
    We like to split per lane because, less barcodes = bigger chance for allowing more mismatches.
    We can't split per lane if:
      - 1 sample is loaded on multiple lanes
    or
      - 1 project is loaded on multiple lanes
    or
      - there are more then 1 lanes, but only 1 is specified in sampleSheet
    """
    laneSplitStatus = True
    # Do we need lane splitting or not ?
    # If there is at least one sample in more then 1 lane, we cannot split:
    if sum(ss['Sample_Name'].value_counts() > 1) > 0:
        laneSplitStatus = False
    # If one project is split over multiple lanes, we also don't split:
    projects = list(ss['Sample_Project'].unique())
    for project in projects:
        if len(
            list(
                ss[ss['Sample_Project'] == project]['Lane'].unique()
                )) > 1:
            laneSplitStatus = False
    if len(list(ss['Lane'].unique())) < lanes:
        laneSplitStatus = False
    return laneSplitStatus

def singleIndex(ss):
    if 'index2' in list(ss.columns):
        return False
    else:
        return True