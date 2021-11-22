import sys
import requests
import logging

def pullParkour(flowcellID, user, pas, URL):
    """
    Look for the flowcell/lane in parkour for the library type 
    """
    d = {'flowcell_id': flowcellID}
    res = requests.get(URL, auth=(user, pas), params=d)
    if res.status_code == 200:
        return res.json()
    return dict()

def main(logFile):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(logFile)


