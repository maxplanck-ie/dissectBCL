#!/usr/bin/env python3
from dissectBCL import preFQ
from dissectBCL import FQ
from rich import print, inspect
from rich.pretty import pprint
import pandas as pd


def main():
    # Read config
    config = preFQ.getConf()
    # Initiate flowcell class (None if no new flowcell found).
    a = preFQ.getNewFlowCell(config)
    if a:
        #################### PREFQ MODULE ####################
        pprint("PREFQ MODULE")
        # Initiate the inferredVars dictionary (until now it was defined but empty)
        a.inferredVars = preFQ.parseRunInfo(a.inVars['runInfo'])
        # Add the sampleSheet vars to inferredVars dictionary
        sampleSheet = preFQ.parseSS(a.inVars['origSS'])
        pprint(sampleSheet)
        a.inferredVars['laneSplitStatus'] = preFQ.decideSplit(sampleSheet, a.inferredVars['lanes'])
        a.inferredVars['singleIndex'] = preFQ.singleIndex(sampleSheet)
        # Create the directories under the output dir.
        a = preFQ.dirCreator(a)
        # Write sampleSheet in the correct way for bclconvert
        inspect(a)

        #################### FQ MODULE ####################
        pprint("FQ MODULE")
    else:
        pprint("Nothing to do. Moving on.")
        

if __name__ == "__main__":
    main()
