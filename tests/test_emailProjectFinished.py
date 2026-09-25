import configparser
from unittest.mock import patch

import pytest

from tools.emailProjectFinished import (
    getContactDetails,
    getFlowCell,
    getFexLinks,
    getProjectIDs,
    main,
    parse_force_to,
)


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
    config["communication"] = {
        "fromaddress": "sender@example.com",
        "bioinfoCore": "core@example.com",
        "host": "smtp.example.com",
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


class Test_parse_force_to:
    def test_parses_recipient_and_pi(self):
        assert parse_force_to("mendelevich@ie-freiburg.mpg.de,iovino") == (
            "mendelevich@ie-freiburg.mpg.de",
            "iovino",
        )

    @pytest.mark.parametrize("value", ["mendelevich@example.com", "a@example.com,one,two"])
    def test_rejects_invalid_format(self, value):
        with pytest.raises(ValueError):
            parse_force_to(value)


class Test_getProjectIDs_force_pi:
    @patch("tools.emailProjectFinished.glob.glob")
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    def test_forced_pi_resolves_sequencing_path(self, mock_fc, mock_glob):
        mock_glob.return_value = [
            "/data/iovino/sequencing_data/FC123",
            "/data/iovino/sequencing_data2/FC123",
        ]

        ids, seqdir = getProjectIDs(
            ["Project_4070_user_externalpi"], _config(), forcePI="iovino"
        )

        assert ids == "4070"
        assert seqdir == "sequencing_data2"
        mock_glob.assert_called_once_with(
            "/data/iovino/sequencing_data*/FC123"
        )


class Test_getFexLinks:
    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    @patch("tools.emailProjectFinished.sp.check_output")
    def test_uses_fex_list_and_detail_to_build_link(self, mock_output, mock_fc):
        project = "Project_4070_Hummel_Domschuke"
        archive = f"FC123_{project}_ro_crate.zip"
        mock_output.side_effect = [
            f"#18) 1 MB [60 d] {archive}\n".encode(),
            b"http://fex.example/fop/ABC123/ABC123?LIST\n",
        ]

        links = getFexLinks([project], _config())

        assert links == {
            project: f"http://fex.example/fop/ABC123/{archive}",
        }
        assert mock_output.call_args_list[0].args[0] == [
            "fexsend",
            "-l",
            "sender@example.com",
        ]
        assert mock_output.call_args_list[1].args[0] == [
            "fexsend",
            "-l",
            "18",
            "sender@example.com",
        ]

    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    @patch("tools.emailProjectFinished.sp.check_output")
    def test_accepts_a_direct_link_in_the_list(self, mock_output, mock_fc):
        project = "Project_4070_Hummel_Domschuke"
        archive = f"FC123_{project}_ro_crate.zip"
        link = f"http://fex.example/fop/ABC123/{archive}"
        mock_output.return_value = f"#18) 1 MB [60 d] {link}\n".encode()

        assert getFexLinks([project], _config()) == {project: link}
        assert mock_output.call_count == 1

    @patch("tools.emailProjectFinished.getFlowCell", return_value="FC123")
    @patch("tools.emailProjectFinished.sp.check_output")
    def test_missing_archive_raises(self, mock_output, mock_fc):
        mock_output.return_value = b"#1) 1 MB [1 d] other.zip\n"

        with pytest.raises(RuntimeError, match="Could not find FEX archive"):
            getFexLinks(["Project_4070_Hummel_Domschuke"], _config())


class Test_main_force_to:
    def test_forces_recipient_and_adds_fex_link_to_comments(self):
        project = "Project_4070_Hummel_Domschuke"
        link = "http://fex.example/fop/ABC123/FC123_Project_4070_Hummel_Domschuke_ro_crate.zip"
        config = _config()
        with (
            patch("tools.emailProjectFinished.getConf", return_value=config),
            patch("tools.emailProjectFinished.getFexLinks", return_value={project: link}),
            patch(
                "tools.emailProjectFinished.getProjectIDs",
                return_value=("4070", "sequencing_data2"),
            ) as mock_project_ids,
            patch("tools.emailProjectFinished.getContactDetails") as mock_contact,
            patch("tools.emailProjectFinished.getFlowCell", return_value="FC123"),
            patch(
                "tools.emailProjectFinished.os.path.exists",
                side_effect=lambda path: path == project,
            ),
            patch("tools.emailProjectFinished.smtplib.SMTP") as mock_smtp,
        ):
            main(
                [
                    "--configfile=config.ini",
                    "--force-to=mendelevich@ie-freiburg.mpg.de,iovino",
                    "--fromPerson=Core",
                    "--fromEmail=core@example.com",
                    "--comment=Project was reviewed.",
                    project,
                ]
            )

        mock_contact.assert_not_called()
        mock_project_ids.assert_called_once_with(
            [project], config, forcePI="iovino"
        )
        assert mock_smtp.call_args.args == ("smtp.example.com",)
        message = mock_smtp.return_value.send_message.call_args.args[0]
        assert message["To"] == "mendelevich@ie-freiburg.mpg.de"
        assert message["Bcc"] == "core@example.com"
        payload = message.get_payload()
        assert "available in your group's sequencing_data2 directory" in payload
        assert "available via FEX" not in payload
        assert "Project was reviewed." in payload
        assert link in payload
        assert "===" in payload

    def test_empty_force_to_is_rejected(self):
        with pytest.raises(SystemExit):
            main(["--force-to=", "Project_4070_Hummel_Domschuke"])
