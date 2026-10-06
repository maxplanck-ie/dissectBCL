import configparser
from unittest.mock import patch

from click.testing import CliRunner

from wd40.reset import collect_targets, detect_mode, reset
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


def _make_illumina_outlane(base):
    (base / "demuxSheet.csv").write_text("sheet")
    (base / "Reports").mkdir()
    (base / "Reports" / "Demultiplex_Stats.csv").write_text("x")
    (base / "Logs").mkdir()
    (base / "Logs" / "Info.log").write_text("x")
    (base / "Project_1_Foo").mkdir()
    (base / "Project_1_Foo" / "sample.fastq.gz").write_text("x")
    (base / "FASTQC_Project_1_Foo").mkdir()
    (base / "Undetermined_S0_L001_R1_001.fastq.gz").write_text("x")
    (base / "bclconvert.done").write_text("")
    (base / "communication.done").write_text("")
    (base / "renamed.done").write_text("")
    (base / "postmux.done").write_text("")
    (base / "analysis.done").write_text("")
    (base / "run.failed").write_text("")
    (base / "fastq.made").write_text("")
    return base


def _make_aviti_outlane(base):
    (base / "manifest").mkdir()
    (base / "manifest" / "RunManifest.csv").write_text("manifest")
    (base / "RunManifest.csv").write_text("b2f-copy")
    (base / "Samples").mkdir()
    (base / "Samples" / "Project_1_Foo").mkdir()
    (base / "info").write_text("x")
    (base / "RunManifest.json").write_text("x")
    (base / "RunParameters.json").write_text("x")
    (base / "IndexAssignment.csv").write_text("x")
    (base / "UnassignedSequences.csv").write_text("x")
    (base / "Metrics.csv").write_text("x")
    (base / "multiqc_report.html").write_text("x")
    (base / "multiqc_data").mkdir()
    (base / "Logs").mkdir()
    (base / "Project_1_Foo").mkdir()
    (base / "FASTQC_Project_1_Foo").mkdir()
    (base / "RunStats.json").write_text("x")
    (base / "bases2fastq.done").write_text("")
    (base / "communication.done").write_text("")
    (base / "analysis.done").write_text("")
    (base / ".1_Foo.renamed.done").write_text("")
    (base / ".1_Foo.postmux.done").write_text("")
    (base / "run.failed").write_text("")
    (base / "fastq.made").write_text("")
    return base


class Test_detect_mode:
    def test_illumina(self, tmp_path):
        _make_illumina_outlane(tmp_path)
        assert detect_mode(tmp_path) == "illumina"

    def test_aviti(self, tmp_path):
        _make_aviti_outlane(tmp_path)
        assert detect_mode(tmp_path) == "aviti"

    def test_neither(self, tmp_path):
        assert detect_mode(tmp_path) is None


class Test_collect_targets:
    def test_illumina_targets(self, tmp_path):
        _make_illumina_outlane(tmp_path)

        mode, targets = collect_targets(tmp_path)

        assert mode == "illumina"
        names = {p.name for p in targets}
        assert names == {
            "Reports",
            "Logs",
            "Project_1_Foo",
            "FASTQC_Project_1_Foo",
            "Undetermined_S0_L001_R1_001.fastq.gz",
            "bclconvert.done",
            "communication.done",
            "renamed.done",
            "postmux.done",
            "analysis.done",
            "run.failed",
            "fastq.made",
        }
        assert (tmp_path / "demuxSheet.csv") not in targets

    def test_aviti_targets(self, tmp_path):
        _make_aviti_outlane(tmp_path)

        mode, targets = collect_targets(tmp_path)

        assert mode == "aviti"
        names = {p.name for p in targets}
        assert names == {
            "Samples",
            "info",
            "RunManifest.json",
            "RunParameters.json",
            "IndexAssignment.csv",
            "UnassignedSequences.csv",
            "Metrics.csv",
            "multiqc_report.html",
            "multiqc_data",
            "Logs",
            "Project_1_Foo",
            "FASTQC_Project_1_Foo",
            "RunStats.json",
            "bases2fastq.done",
            "communication.done",
            "analysis.done",
            ".1_Foo.renamed.done",
            ".1_Foo.postmux.done",
            "run.failed",
            "fastq.made",
            "RunManifest.csv",
        }
        assert (tmp_path / "manifest" / "RunManifest.csv") not in targets

    def test_no_manifest_returns_empty(self, tmp_path):
        mode, targets = collect_targets(tmp_path)

        assert mode is None
        assert targets == []


class Test_reset:
    def test_aborts_on_unrecognized_dir(self, tmp_path, capsys):
        reset(tmp_path)

        out = capsys.readouterr().out.replace("\n", "")
        assert "doesn't look like" in out
        assert "outLane dir" in out
        assert list(tmp_path.iterdir()) == []

    def test_declining_confirmation_deletes_nothing(self, tmp_path):
        _make_illumina_outlane(tmp_path)
        before = sorted(p.name for p in tmp_path.iterdir())

        with patch("wd40.reset.click.confirm", return_value=False):
            reset(tmp_path)

        after = sorted(p.name for p in tmp_path.iterdir())
        assert before == after

    def test_confirming_deletes_illumina_targets_keeps_manifest(self, tmp_path):
        _make_illumina_outlane(tmp_path)

        with patch("wd40.reset.click.confirm", return_value=True):
            reset(tmp_path)

        remaining = {p.name for p in tmp_path.iterdir()}
        assert remaining == {"demuxSheet.csv"}

    def test_confirming_deletes_aviti_targets_keeps_manifest(self, tmp_path):
        _make_aviti_outlane(tmp_path)

        with patch("wd40.reset.click.confirm", return_value=True):
            reset(tmp_path)

        remaining = {p.name for p in tmp_path.iterdir()}
        assert remaining == {"manifest"}
        manifest_remaining = {p.name for p in (tmp_path / "manifest").iterdir()}
        assert manifest_remaining == {"RunManifest.csv"}
        assert (tmp_path / "manifest" / "RunManifest.csv").read_text() == "manifest"

    def test_empty_but_recognized_dir_reports_nothing_to_do(self, tmp_path, capsys):
        (tmp_path / "demuxSheet.csv").write_text("sheet")

        reset(tmp_path)

        out = capsys.readouterr().out
        assert "Nothing to reset" in out
        assert list(tmp_path.iterdir()) == [tmp_path / "demuxSheet.csv"]


def test_cli_reset_invokes_reset_outLane(tmp_path):
    configfile = tmp_path / "conf.ini"
    configfile.write_text("[dummy]\nkey=val\n")

    with (
        patch("wd40.wd40.getConf", return_value=_config()),
        patch("wd40.wd40.reset_outLane") as mock_reset,
    ):
        runner = CliRunner()
        result = runner.invoke(
            cli, ["--configpath", str(configfile), "reset", str(tmp_path)]
        )

    assert result.exit_code == 0, result.output
    mock_reset.assert_called_once_with(str(tmp_path))
