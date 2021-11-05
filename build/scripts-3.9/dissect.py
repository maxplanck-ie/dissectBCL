#!python
from dissectBCL import preFQ
from dissectBCL import FQ
from rich import print, inspect

def main():
    config = preFQ.getConf()
    a = preFQ.getNewFlowCell(config)
    if a:
        a = preFQ.parseRunInfo( a )
        a = preFQ.parseSS( a )
        inspect( a )
        FQ.dirCreator( a )


if __name__ == "__main__":
    main()