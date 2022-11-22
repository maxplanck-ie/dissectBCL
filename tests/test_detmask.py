import pandas as pd
import os
from dissectBCL.demux import detMask

class Test_detmask_Files():
    def readss(self, ss):
        sspath = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'test_detmask',
            ss
        )
        return (
            pd.read_csv(
                sspath,
                sep='\t',
                keep_default_na=False
            )
        )

    def test_noindex(self):
        sR = {
            'Read1': ['Y', 37]
        }
        ss = self.readss('noindex.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y37'
        assert dualIx == False
        assert PE == False
        assert convOpts == []
        assert minP5 == None
        assert minP7 == None

    def test_scATAC(self):
        sR = {
            'Read1': ['Y', 101],
            'Index1': ['I', 10],
            'Index2': ['I', 24],
            'Read2': ['Y', 101]
        }
        ss = self.readss('scatac.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y101;I8N2;U24;Y101' # P5 = 24bp UMI.
        assert dualIx == False
        assert PE == True
        assert convOpts == ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
    
    def test_nugen(self):
        sR = {
            'Read1': ['Y', 51],
            'Index1': ['I', 16],
            'Read2': ['Y', 51]
        }
        ss = self.readss('nugen.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y51;I8U8;Y51' # P7 = 8bp index, 8 index UMI
        assert dualIx == False
        assert PE == True
        assert convOpts == ['CreateFastQForIndexReads,1,,', 'TrimUMI,0,,']
    
    def test_dualix1(self):
        sR = {
            'Read1': ['Y', 101],
            'Index1': ['I', 10],
            'Index2': ['I', 10],
            'Read2': ['Y', 101]
        }
        ss = self.readss('dualIx_1.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y101;I8N2;I8N2;Y101'
        assert dualIx == True
        assert PE == True
        assert convOpts == []
        assert minP5 == 8
        assert minP7 == 8
    
    def test_dualix2(self):
        sR = {
            'Read1': ['Y', 150],
            'Index1': ['I', 8],
            'Index2': ['I', 8],
            'Read2': ['Y', 150]
        }
        ss = self.readss('dualIx_2.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y150;I8;I8;Y150'
        assert dualIx == True
        assert PE == True
        assert convOpts == []
        assert minP5 == 8
        assert minP7 == 8
    
    def test_dualix3(self):
        sR = {
            'Read1': ['Y', 101],
            'Index1': ['I', 10],
            'Index2': ['I', 24],
            'Read2': ['Y', 101]
        }
        ss = self.readss('dualIx_3.tsv')
        mask, dualIx, PE, convOpts, minP5, minP7 = detMask(sR, ss, 'out')
        assert mask == 'Y101;I8N2;I8N16;Y101'
        assert dualIx == True
        assert PE == True
        assert convOpts == []
        assert minP5 == 8
        assert minP7 == 8
    


