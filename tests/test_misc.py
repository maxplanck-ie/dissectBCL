import pytest
import pandas as pd
from dissectBCL.misc import joinLis
from dissectBCL.misc import hamming
from dissectBCL.misc import lenMask
from dissectBCL.misc import P5Seriesret

# Fake DFs
adf = pd.DataFrame(
    {
        'index': [1,2,3,4,5],
        'index2': [1,2,3,4,5]
    }
)

bdf = pd.DataFrame(
    {
        'index':[1,2,3,4,5]
    }
)


def test_joinLis():
    assert joinLis([1, 2, 'A']) == "12A"
    assert joinLis([1, 2, 3]) == "123"
    assert joinLis(['a', 2, 'b'], joinStr = ',') == 'a,2,b'

def test_hamming():
    assert hamming(float(1), float(2)) == 0
    assert hamming('aabb', 'aaba') == 1
    assert hamming('aaaa', 'bbbb') == 4

def test_lenMask():
    assert lenMask(8,4) == "I4N4"
    assert lenMask(10,10) == "I10"

def test_P5Seriesret():
    serRet = P5Seriesret(adf)
    for i,j in zip(serRet, adf['index2']):
        assert i == j
    serRet = P5Seriesret(bdf)
    assert serRet.empty


# tests that need files/folders.
# getConf
# getNewFlowCell
# parseRunInfo
# screenFqFetcher
# moveoptdup
# fetchLatestSeqDir


