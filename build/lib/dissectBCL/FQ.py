import rich
import os
import sys

def dirCreator(flowClass):
    # Do we split lanes ?
    if flowClass.laneSplitStatus:
        tmpLis = []
        for lane in range(flowClass.lanes):
            tmpLis.append( flowClass.name + '_lanes' + str(lane + 1) )
        flowClass.outDir = tmpLis
        rich.print(flowClass.outDir)