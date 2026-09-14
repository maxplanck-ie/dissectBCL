import configparser
from unittest.mock import patch

from click.testing import CliRunner

from wd40.wd40 import cli


def _config():
    config = configparser.ConfigParser()
    config["Dirs"] = {"piDir": "/pidir"}
    config["Internals"] = {
        "PIs": "goodpi,badpi",
        "seqDir": "sequencing_data",
        "fex": "false",
    }
    config["parkour"] = {
        "URL": "http://parkour",
        "user": "u",
        "password": "p",
        "cert": "",
    }
    config["communication"] = {"fromAddress": "from@x.com"}
    return config


def test_cli_populates_ctx_from_config(tmp_path):
    configfile = tmp_path / "conf.ini"
    configfile.write_text("[dummy]\nkey=val\n")

    with (
        patch("wd40.wd40.getConf", return_value=_config()) as mock_getConf,
        patch("wd40.wd40.release") as mock_release,
    ):
        runner = CliRunner()
        result = runner.invoke(
            cli, ["--configpath", str(configfile), "rel", str(tmp_path)]
        )

    assert result.exit_code == 0, result.output
    mock_getConf.assert_called_once_with(str(configfile), quickload=True)
    mock_release.assert_called_once_with(
        str(tmp_path),
        "goodpi,badpi",
        "/pidir",
        "sequencing_data",
        "http://parkour",
        ("u", "p"),
        "",
        False,
        "from@x.com",
    )


def test_cli_debug_flag_sets_ctx(tmp_path):
    configfile = tmp_path / "conf.ini"
    configfile.write_text("[dummy]\nkey=val\n")

    captured = {}

    def fake_rel(*args):
        captured["args"] = args

    with (
        patch("wd40.wd40.getConf", return_value=_config()),
        patch("wd40.wd40.release", side_effect=fake_rel),
    ):
        runner = CliRunner()
        result = runner.invoke(
            cli, ["--configpath", str(configfile), "--debug", "rel", str(tmp_path)]
        )

    assert result.exit_code == 0, result.output
    assert captured["args"][0] == str(tmp_path)
