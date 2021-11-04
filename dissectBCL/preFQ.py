import configparser
import os
import sys
import rich

def getConf():
    # Get userDir
    homeDir = os.path.expanduser("~")
    confLoc = os.path.join(homeDir, 'dissectBCL.ini')
    if not os.path.exists(confLoc):
        sys.stderr.write("Ini file not found. Exiting.")
        sys.exit(1)
    else:
        config = configparser.ConfigParser()
        config.read( confLoc )
        return config

def getNewFlowCell(config):
    baseDir = config.baseDir
    outDir = config.outDir
    

config = getConf()
print(config['Dirs']['baseDir'])



