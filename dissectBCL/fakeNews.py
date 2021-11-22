import sys
import requests
import logging

log = logging.getLogger()

def setLog(logFile):
    logging.basicConfig(
        filename = logFile,
        level="DEBUG",
        format="%(levelname)s    %(asctime)s    %(message)s",
        filemode = 'w'
    )
    log = logging.getLogger()

def pullParkour(flowcellID, user, pas, URL):
    """
    Look for the flowcell/lane in parkour for the library type 
    """
    d = {'flowcell_id': flowcellID}
    res = requests.get(URL, auth=(user, pas), params=d)
    if res.status_code == 200:
        return res.json()
    return dict()
