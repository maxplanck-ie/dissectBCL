#!/usr/bin/env python3
from dissectBCL import preFQ, FQ, fakeNews, misc
from dissectBCL.classes import sampleSheetClass
from rich import print, inspect
from rich.pretty import pprint
import pandas as pd
import argparse
import logging


def main():
    # Read config
    config = misc.getConf()
    # Search flowcell and initiate flowcell Class if found.
    a = preFQ.getNewFlowCell(config)

    if a:
        #################### PREFQ MODULE ####################
        #a.inferredVars = misc.parseRunInfo(a.runInfo)
        # Read the sampleSheet and invoke sampleSheet class.
        #sampleSheet = sampleSheetClass( misc.parseSS(a.origSS), a.lanes )
        # Do we need to split lanes ?
        #sampleSheetClass.decideSplit(sampleSheet)
        # Infer if we need to split lanes
        #a.inferredVars['laneSplitStatus'] = preFQ.decideSplit(sampleSheet, a.inferredVars['lanes'])
        # Test if we have single index or dual index
        #a.inferredVars['singleIndex'] = preFQ.singleIndex(sampleSheet)
        # Create the directories under the output dir.
        #a = preFQ.dirCreator(a)
        # Read parkour
        #a.pullParkour = fakeNews.pullParkour(a.inferredVars['flowcellID'], config['parkour']['user'], config['parkour']['password'], config['parkour']['pullURL'])
        inspect(a)
        #inspect(sampleSheet)

        #################### FQ MODULE ####################
        pprint("FQ MODULE")
    else:
        pprint("Nothing to do. Moving on.")
        

if __name__ == "__main__":
    main()
