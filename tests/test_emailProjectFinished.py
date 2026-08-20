import configparser
import json
from unittest.mock import patch

import pytest

from tools.emailProjectFinished import getProjectIDs


def _config(tmp_path, PIs, deliver_to=None):
    config = configparser.ConfigParser()
    config["Internals"] = {
        "PIs": PIs,
        "seqDir": "sequencing_data",
        "deliverTo": json.dumps(deliver_to or {}),
    }
    config["Dirs"] = {"piDir": str(tmp_path)}
    return config


class Test_getProjectIDs:
    @patch("tools.emailProjectFinished.getFlowCell")
    def test_deliver_to_key_resolves_even_when_absent_from_pi_list(
        self, mock_getFlowCell, tmp_path
    ):
        # Project string still carries the PI's previous Parkour name; only
        # a deliverTo override maps it to the current /data directory.
        flowcell = "260818_A00931_0932_BHNHMNDRX7_lanes_2"
        mock_getFlowCell.return_value = flowcell
        seqDir = tmp_path / "cabezas" / "sequencing_data"
        (seqDir / flowcell).mkdir(parents=True)
        config = _config(
            tmp_path, "cabezas,akhtar", {"cabezas-wallscheid": "cabezas"}
        )

        projectID, seqdir = getProjectIDs(
            ["Project_4035_Demollin_Cabezas-Wallscheid"], config
        )

        assert projectID == "4035"
        assert seqdir == "sequencing_data"

    @patch("tools.emailProjectFinished.getFlowCell")
    def test_pi_not_internal_and_no_override_exits(self, mock_getFlowCell, tmp_path):
        mock_getFlowCell.return_value = "260101_M001_0001_AAAA"
        config = _config(tmp_path, "manke")

        with pytest.raises(SystemExit):
            getProjectIDs(["Project_1_jdoe_unknownpi"], config)

    @patch("tools.emailProjectFinished.getFlowCell")
    def test_simple_internal_pi_still_works(self, mock_getFlowCell, tmp_path):
        flowcell = "260101_M001_0001_AAAA"
        mock_getFlowCell.return_value = flowcell
        seqDir = tmp_path / "manke" / "sequencing_data"
        (seqDir / flowcell).mkdir(parents=True)
        config = _config(tmp_path, "manke")

        projectID, seqdir = getProjectIDs(["Project_1_jdoe_manke"], config)

        assert projectID == "1"
        assert seqdir == "sequencing_data"
