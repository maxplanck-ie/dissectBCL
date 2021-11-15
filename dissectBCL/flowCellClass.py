import os

# Define flowcell class, which will contain all information
class flowCell:
    """This is a flowCell class, which contains:
    - prior information set in the config file and read from the flowcell directory (inVars).
    - inferred variables from the RunInfo.xml and the sampleSheet (inferredVars).
    - variables pulled from parkour (json, pullParkour)
    """

    # fileChecks
    def filesExist(self):
        """
        takes a list of files and returns a tuple (Bool, Dict) with:
        Bool = True (all files in files exist)
        Bool = False (at least 1 file doesn't exist)
        Dict = Dict with bool status of each file.
        """
        checkFiles = [
            self.inVars['bclPath'],
            self.inVars['origSS'],
            self.inVars['runInfo'],
            self.inVars['inBaseDir'],
            self.inVars['outBaseDir']
        ]
        existDic = {}
        existStatus = True
        for f in checkFiles:
            if not os.path.exists(f):
                existDic[f] = False
                existStatus = False
            else:
                existDic[f] = True
        return (existStatus, existDic)


    def __init__(self, name, bclPath, origSS, runInfo, inBaseDir, outBaseDir, logOut, logErr):
        sequencers = {
            'A': 'NovaSeq',
            'N': 'NextSeq',
            'M': 'MiSeq'
        }
        self.inVars = {
            'name': name,
            'sequencer': sequencers[name.split('_')[1][0]],
            'bclPath': bclPath,
            'origSS': origSS,
            'runInfo': runInfo,
            'inBaseDir': inBaseDir,
            'outBaseDir': outBaseDir,
            'logOut': logOut,
            'logErr': logErr
        }
        self.inferredVars = {}
        self.pullParkour = {}