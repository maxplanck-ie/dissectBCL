import configparser
from contextlib import suppress
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from dissectBCL.dissect import createFlowcell, dissect, main


def _config(logDir, sequencer):
    config = configparser.ConfigParser()
    config["Dirs"] = {f"flowLogDir_{sequencer}": str(logDir)}
    config["communication"] = {"debug_mode": "false"}
    config["softwareVers"] = {}
    return config


def test_dissect_command_loads_conf_and_calls_main(tmp_path):
    configfile = tmp_path / "conf.ini"
    configfile.write_text("[dummy]\nkey=val\n")

    with (
        patch("dissectBCL.dissect.getConf", return_value="theconfig") as mock_getConf,
        patch("dissectBCL.dissect.main") as mock_main,
        patch("dissectBCL.dissect.getVersion", return_value="1.2.3"),
    ):
        runner = CliRunner()
        result = runner.invoke(
            dissect,
            ["--configfile", str(configfile), "--sequencer", "illumina"],
        )

    assert result.exit_code == 0, result.output
    mock_getConf.assert_called_once_with(str(configfile), sequencer="illumina")
    mock_main.assert_called_once_with("theconfig", None, "illumina", False)


def test_main_processes_one_flowcell_then_stops(tmp_path, monkeypatch):
    logDir = tmp_path / "logs"
    config = _config(logDir, "illumina")

    flowcell = MagicMock()
    calls = {"n": 0}

    def fake_getNewFlowCell(config, flowcellpath, platformFilter):
        calls["n"] += 1
        if calls["n"] == 1:
            return "FC1", "/data/FC1", "illumina"
        raise StopIteration

    with (
        patch("dissectBCL.dissect.getNewFlowCell", side_effect=fake_getNewFlowCell),
        patch("dissectBCL.dissect.flowCellClass", return_value=flowcell) as mock_fc,
        patch("dissectBCL.dissect.getVersion", return_value="1.2.3"),
        suppress(StopIteration),
    ):
        main(config, None, "illumina", False)

    mock_fc.assert_called_once_with(
        name="FC1",
        bclPath="/data/FC1",
        logFile=logDir / "FC1.log",
        config=config,
        sequencer="illumina",
        forceLaneSplit=False,
    )
    flowcell.prepConvert.assert_called_once()
    flowcell.demux.assert_called_once()
    flowcell.demux_aviti.assert_not_called()
    flowcell.postmux.assert_called_once()
    flowcell.fakenews.assert_called_once()
    flowcell.organiseLogs.assert_called_once()


def test_main_aviti_flowcell_nests_logdir_and_uses_demux_aviti(tmp_path):
    logDir = tmp_path / "logs"
    config = _config(logDir, "aviti")

    flowcell = MagicMock()
    calls = {"n": 0}

    def fake_getNewFlowCell(config, flowcellpath, platformFilter):
        calls["n"] += 1
        if calls["n"] == 1:
            return "FC2", "/data/AV251009/FC2", "aviti"
        raise StopIteration

    with (
        patch("dissectBCL.dissect.getNewFlowCell", side_effect=fake_getNewFlowCell),
        patch("dissectBCL.dissect.flowCellClass", return_value=flowcell) as mock_fc,
        patch("dissectBCL.dissect.getVersion", return_value="1.2.3"),
        suppress(StopIteration),
    ):
        main(config, None, "aviti", False)

    mock_fc.assert_called_once_with(
        name="FC2",
        bclPath="/data/AV251009/FC2",
        logFile=logDir / "AV251009" / "FC2.log",
        config=config,
        sequencer="aviti",
        forceLaneSplit=False,
    )
    flowcell.demux_aviti.assert_called_once()
    flowcell.demux.assert_not_called()


def test_createFlowcell_no_logfile_defaults_to_stdout():
    flowcell = MagicMock()
    with (
        patch("dissectBCL.dissect.getConf", return_value="theconfig") as mock_getConf,
        patch(
            "dissectBCL.dissect.getNewFlowCell",
            return_value=("FC1", "/data/FC1", "illumina"),
        ),
        patch("dissectBCL.dissect.flowCellClass", return_value=flowcell) as mock_fc,
    ):
        result = createFlowcell("conf.ini", "/data/FC1", "illumina")

    mock_getConf.assert_called_once_with("conf.ini", sequencer="illumina")
    mock_fc.assert_called_once_with(
        name="FC1",
        bclPath="/data/FC1",
        logFile="STDOUT",
        config="theconfig",
        sequencer="illumina",
        forceLaneSplit=False,
    )
    assert result is flowcell


def test_createFlowcell_with_logfile_creates_parent_dir(tmp_path):
    logFile = tmp_path / "sub" / "fc.log"
    flowcell = MagicMock()
    with (
        patch("dissectBCL.dissect.getConf", return_value="theconfig"),
        patch(
            "dissectBCL.dissect.getNewFlowCell",
            return_value=("FC1", "/data/FC1", "illumina"),
        ),
        patch("dissectBCL.dissect.flowCellClass", return_value=flowcell) as mock_fc,
    ):
        createFlowcell("conf.ini", "/data/FC1", "illumina", logFile=str(logFile))

    assert logFile.parent.is_dir()
    mock_fc.assert_called_once_with(
        name="FC1",
        bclPath="/data/FC1",
        logFile=str(logFile),
        config="theconfig",
        sequencer="illumina",
        forceLaneSplit=False,
    )
