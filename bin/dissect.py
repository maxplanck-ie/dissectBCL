#!/usr/bin/env python3
from dissectBCL import preFQ
from rich import print, inspect

def main():
    config = preFQ.getConf()
    a = preFQ.getNewFlowCell(config)
    a = preFQ.parseRunInfo( a )
    a = preFQ.parseSS( a )
    inspect( a )

if __name__ == "__main__":
    main()