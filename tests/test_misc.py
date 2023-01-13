import pandas as pd
import os
from dissectBCL.misc import joinLis
from dissectBCL.misc import hamming
from dissectBCL.misc import lenMask
from dissectBCL.misc import P5Seriesret
from dissectBCL.misc import moveOptDup
from dissectBCL.misc import retBCstr
from dissectBCL.misc import retIxtype
from dissectBCL.misc import retMean_perc_Q
from dissectBCL.misc import formatSeqRecipe
from dissectBCL.misc import formatMisMatches
from dissectBCL.misc import umlautDestroyer
from dissectBCL.misc import parseRunInfo

class Test_misc_data():
    def test_joinLis(self):
        assert joinLis([1, 2, 'A']) == "12A"
        assert joinLis([1, 2, 3]) == "123"
        assert joinLis(['a', 2, 'b'], joinStr=',') == 'a,2,b'

    def test_hamming(self):
        assert hamming(float(1), float(2)) == 0
        assert hamming('aabb', 'aaba') == 1
        assert hamming('aaaa', 'bbbb') == 4

    def test_lenMask(self):
        assert lenMask(8, 4) == "I4N4"
        assert lenMask(10, 10) == "I10"

    def test_P5Seriesret(self):
        adf = pd.DataFrame(
            {
                'index': [1, 2, 3, 4, 5],
                'index2': [1, 2, 3, 4, 5]
            }
        )

        bdf = pd.DataFrame(
            {
                'index': [1, 2, 3, 4, 5]
            }
        )
        for i, j in zip(P5Seriesret(adf), adf['index2']):
            assert i == j
        assert P5Seriesret(bdf).empty

    def test_retBCstr(self):
        _a = pd.Series(
            data=[1, 1],
            index=['index', 'index2']
        )
        _b = pd.Series(
            data=[1, 'bar'],
            index=['index', 'foo']
        )
        _c = pd.Series(
            data=1,
            index=['foo']
        )
        assert retBCstr(_a) == '1\t1'
        assert retBCstr(_b) == '1'
        assert retBCstr(_c) == 'nan'

    def test_retIxtype(self):
        _a = pd.Series(
            data=['I7type', 'I5type'],
            index=['I7_Index_ID', 'I5_Index_ID']
        )
        _b = pd.Series(
            data=['I7type'],
            index=['I7_Index_ID']
        )
        _c = pd.Series(
            data=['foo'],
            index=['bar']
        )
        assert retIxtype(_a) == 'I7type\tI5type'
        assert retIxtype(_b) == 'I7type'
        assert retIxtype(_c) == 'NA'

    def test_retMean_perc_Q(self):
        _a = pd.Series(
            data=['1:26.1,I1:22.8,2:23.0'],
            index=['meanQ']
        )
        _b = pd.Series(
            data=['1:48.2'],
            index=['meanQ']
        )
        assert retMean_perc_Q(_a) == '26.0\t23.0\t23.0'
        assert retMean_perc_Q(_b) == '48.0'
        assert retMean_perc_Q(_a, returnHeader=True) == (
            'R1_meanQ\tI1_meanQ\tR2_meanQ',
            '26.0\t23.0\t23.0'
        )

    def test_formatSeqRecipe(self):
        _a = {
            'Read1': ['Y', 100],
            'Read2': ['Y', 100]
        }
        _b = {
            'Read1': ['Y', 51],
            'Index1': ['I', 10]
        }
        assert formatSeqRecipe(_a) == "Read1:100; Read2:100"
        assert formatSeqRecipe(_b) == "Read1:51; Index1:10"

    def test_formatMisMatches(self):
        _a = {
            'BarcodeMismatchesIndex1': 2,
            'BarcodeMismatchesIndex2': 1
        }
        _b = {
            'BarcodeMismatchesIndex1': 1
        }
        assert formatMisMatches(
            _a
        ) == "BarcodeMismatchesIndex1:2, BarcodeMismatchesIndex2:1"
        assert formatMisMatches(_b) == "BarcodeMismatchesIndex1:1"

    def test_umlautDestroyer(self):
        _a = "ö"
        _b = "ä"
        _c = "ß"
        assert umlautDestroyer(_a) == "o"
        assert umlautDestroyer(_b) == "a"
        assert umlautDestroyer(_c) == "ss"


class Test_misc_Files():
    def RTF(self, testFile):
        return os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)
            ),
            'test_misc',
            testFile
        )

    def test_moveOptDup(self):
        moveOptDup(self.RTF('myLane'))
        _file_in = "myLane/x1/x2/x3_duplicate.txt"
        _file_out = "myLane/FASTQC_x1/x2/x3_duplicate.txt"
        assert os.path.exists(self.RTF(_file_out))
        os.rename(self.RTF(_file_out), self.RTF(_file_in))
        assert os.path.exists(self.RTF(_file_in))

    def test_parseRunInfo(self):
        _runInfo = parseRunInfo(
            self.RTF("RunInfo.xml")
        )
        _readDic = {
            'Read1': ['150', 'Read'],
            'Read2': ['8', 'Index'],
            'Read3': ['8', 'Index'],
            'Read4': ['150', 'Read']
        }
        assert _runInfo['instrument'] == 'NB000000'
        assert _runInfo['readDic'] == _readDic
        assert _runInfo['lanes'] == 4
        assert _runInfo['flowcellID'] == 'HHHHHHHHH'

