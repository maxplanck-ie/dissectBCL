#!/usr/bin/env python3
from dissectBCL import preFQ
from dissectBCL import FQ
from rich import print, inspect


def main():
    config = preFQ.getConf()
    a = preFQ.getNewFlowCell(config)
    if a:
        print("New flowcell found.")
        # Initiate the inferredVars dictionary (until now it was defined but empty)
        a.inferredVars = preFQ.parseRunInfo(a.inVars['runInfo'])
        # Add the sampleSheet vars to inferredVars dictionary
        a.inferredVars = a.inferredVars | preFQ.parseSS(a.inVars['origSS'])
        inspect(a)
        #FQ.dirCreator(a)


if __name__ == "__main__":
    main()
