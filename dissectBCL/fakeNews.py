import sys
import requests

def pullParkour(flowcellID, user, pas, URL):
    """
    Look for the flowcell/lane in parkour for the library type 
    """
    d = {'flowcell_id': flowcellID}
    res = requests.get(URL, auth=(user, pas), params=d)
    if res.status_code == 200:
        return res.json()
    return dict()