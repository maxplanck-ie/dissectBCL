#!/usr/bin/env python3
from dissectBCL import preFQ, FQ, fakeNews
from rich import print, inspect
from rich.pretty import pprint
import pandas as pd
import argparse


def main():
    # Read config
    config = preFQ.getConf()
    # Initiate flowcell class (None if no new flowcell found).
    a = preFQ.getNewFlowCell(config)
    if a:
        #################### PREFQ MODULE ####################

        # Initiate the inferredVars dictionary (until now it was defined but empty)
        a.inferredVars = preFQ.parseRunInfo(a.inVars['runInfo'])
        # Read the sampleSheet.
        sampleSheet = preFQ.parseSS(a.inVars['origSS'])
        pprint(sampleSheet)
        # Infer if we need to split lanes
        a.inferredVars['laneSplitStatus'] = preFQ.decideSplit(sampleSheet, a.inferredVars['lanes'])
        # Test if we have single index or dual index
        a.inferredVars['singleIndex'] = preFQ.singleIndex(sampleSheet)
        # Create the directories under the output dir.
        a = preFQ.dirCreator(a)
        # Read parkour
        a.pullParkour = fakeNews.pullParkour(a.inferredVars['flowcellID'], config['parkour']['user'], config['parkour']['password'], config['parkour']['pullURL'])
        inspect(a)

        #################### FQ MODULE ####################
        pprint("FQ MODULE")
    else:
        pprint("Nothing to do. Moving on.")
        

if __name__ == "__main__":
    main()
