import configparser
from unittest.mock import patch

import pytest

from tools.emailProjectFinished import getContactDetails, getFlowCell, getProjectIDs


def _config(pis="goodpi,otherpi"):
    config = configparser.ConfigParser()
    config["Internals"] = {"PIs": pis, "seqDir": "sequencing_data"}
    config["Dirs"] = {"piDir": "/data"}
    config["parkour"] = {
        "URL": "https://parkour.example.com",
        "user": "u",
        "password": "p",
        "cert": "",
    }
    return config


class Test_getFlowCell:
    def test_returns_last_path_component_of_cwd(self):
        with patch("tools.emailProjectFinished.os.getcwd", return_value="/a/b/FC123"):
            assert getFlowCell() == "FC123"


class Test_getContactDetails:
    @patch("tools.emailProjectFinished.requests.get")
    def test_returns_json_on_200(self, mock_get):
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {"email": "x@example.com"}

        result = getContactDetails("1234", _config())

        assert result == {"email": "x@example.com"}

    @patch("tools.emailProjectFinished.requests.get")
    def test_raises_on_non_200(self, mock_get):
        mock_get.return_value.status_code = 500
        mock_get.return_value.json.return_value = {"detail": "boom"}

        with pytest.raises(RuntimeError, match="boom"):
            getContactDetails("1234", _config())


class Test_getProjectIDs:
    @patch("tools.emailProjectFinished.glob.glob")
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    def test_internal_pi_single_project(self, mock_fc, mock_glob):
        mock_glob.return_value = ["/data/goodpi/sequencing_data_2024/FC123"]

        ids, seqdir = getProjectIDs(["Project_1234_jdoe_goodpi"], _config())

        assert ids == "1234"
        assert seqdir == "sequencing_data_2024"

    @patch("tools.emailProjectFinished.glob.glob")
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    def test_internal_pi_multiple_projects_joined_with_and(self, mock_fc, mock_glob):
        mock_glob.return_value = ["/data/goodpi/sequencing_data_2024/FC123"]

        ids, _ = getProjectIDs(
            [
                "Project_1_jdoe_goodpi",
                "Project_2_jdoe_goodpi",
                "Project_3_jdoe_goodpi",
            ],
            _config(),
        )

        assert ids == "1, 2 and 3"

    @patch("tools.emailProjectFinished.glob.glob")
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    def test_external_pi_exits(self, mock_fc, mock_glob):
        with pytest.raises(SystemExit, match="not in the internal PI list"):
            getProjectIDs(["Project_1234_jdoe_externalpi"], _config())
        mock_glob.assert_not_called()

    @patch("tools.emailProjectFinished.glob.glob", return_value=[])
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    def test_no_matching_sequencing_dir_exits(self, mock_fc, mock_glob):
        with pytest.raises(SystemExit, match="No sequencing_data directory"):
            getProjectIDs(["Project_1234_jdoe_goodpi"], _config())
