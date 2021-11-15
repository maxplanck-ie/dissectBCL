import rich
import os
import sys

def bclConverter(flowCellClass):
    if flowCellClass.inferredVars['laneSplitStatus']:
        

    [Header],,,
FileFormatVersion,2,,
,,,
[BCLConvert_Settings],,,
BarcodeMismatchesIndex1,1,,
BarcodeMismatchesIndex2,1,,
OverrideCycles,Y101;I8N2;I8N16;Y101,,

,,,
[BCLConvert_Data],,,,,,
Lane,Sample_ID,index,index2,Sample_Project
1,21L003793,TAAGGCGA,TACTCCTT,2075_Schoenberger_Cabezas-Wallscheid
1,21L003794,CGTACTAG,AGGCTTAG,2075_Schoenberger_Cabezas-Wallscheid
1,21L003795,AGGCAGAA,ATTAGACG,2075_Schoenberger_Cabezas-Wallscheid