import configparser
from unittest.mock import patch

import pandas as pd

from dissectBCL.postmux import kraken, runPlusPF


def _config():
    c = configparser.ConfigParser()
    c["misc"] = {"threads": "10"}
    c["screening"] = {"plusPFdb": "/fake/pluspf"}
    return c


def _make_sample(laneFolder, project, sampleID):
    sampleFolder = laneFolder / f"Project_{project}" / f"Sample_{sampleID}"
    sampleFolder.mkdir(parents=True)
    (sampleFolder / f"{sampleID}_R1.fastq.gz").write_bytes(b"fake")
    (laneFolder / f"FASTQC_Project_{project}" / f"Sample_{sampleID}").mkdir(
        parents=True
    )
    return sampleFolder


class _SyncPool:
    """
    Stand-in for multiprocessing.Pool that runs .map() synchronously in the
    current process (no fork), so mocks patched in the test are visible to
    the code under test and their call history is inspectable afterwards.
    """

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, fn, iterable):
        return [fn(x) for x in iterable]


class Test_runPlusPF:
    @patch("dissectBCL.postmux.Pool", _SyncPool)
    @patch("dissectBCL.postmux.Popen")
    def test_writes_a_plusPF_report_per_flagged_sample(self, mock_popen, tmp_path):
        mock_popen.return_value.wait.return_value = 0
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")

        runPlusPF("1_proj", laneFolder, ["S1"], _config())

        assert mock_popen.called
        cmd = mock_popen.call_args[0][0]
        assert cmd[0] == "kraken2"
        assert "/fake/pluspf" in cmd
        assert any(arg.endswith(".plusPF.krakenreport") for arg in cmd)

    @patch("dissectBCL.postmux.Pool", _SyncPool)
    @patch("dissectBCL.postmux.mailHome")
    @patch("dissectBCL.postmux.Popen")
    def test_failed_run_emails_core_and_cleans_up_partial_report(
        self, mock_popen, mock_mailHome, tmp_path
    ):
        mock_popen.return_value.wait.return_value = 1
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        # kraken2 can write a partial --report file before later failing;
        # simulate that here to check runPlusPF cleans it up.
        partialReport = (
            laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.plusPF.krakenreport"
        )
        partialReport.write_text("truncated")

        runPlusPF("1_proj", laneFolder, ["S1"], _config())

        mock_mailHome.assert_called_once()
        assert not partialReport.exists()

    @patch("dissectBCL.postmux.Pool", _SyncPool)
    @patch("dissectBCL.postmux.Popen")
    def test_no_flagged_samples_does_not_call_kraken2(self, mock_popen, tmp_path):
        laneFolder = tmp_path / "lane"
        runPlusPF("1_proj", laneFolder, [], _config())
        mock_popen.assert_not_called()


class Test_kraken_escalation:
    def _config(self):
        c = configparser.ConfigParser()
        c["misc"] = {"threads": "10"}
        c["software"] = {"kraken2db": "/fake/krakendb"}
        c["screening"] = {
            "plusPFdb": "/fake/pluspf",
            "unclassified_threshold": "10",
            "relaxed_library_types": "ATAC-Seq",
            "relaxed_threshold": "20",
        }
        return c

    def _ssdf(self, sampleID, libraryType):
        return pd.DataFrame(
            {"Sample_ID": [sampleID], "Library_Type": [libraryType]}
        )

    @patch("dissectBCL.postmux.runPlusPF")
    def test_flags_sample_over_default_threshold_for_escalation(
        self, mock_runPlusPF, tmp_path
    ):
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        reportPath = (
            laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep"
        )
        reportPath.write_text("15.0\t100\t100\tU\t0\tunclassified\n")

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config())

        mock_runPlusPF.assert_called_once()
        called_ids = mock_runPlusPF.call_args[0][2]
        assert called_ids == ["S1"]

    @patch("dissectBCL.postmux.runPlusPF")
    def test_does_not_flag_atac_sample_under_relaxed_threshold(
        self, mock_runPlusPF, tmp_path
    ):
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        reportPath = (
            laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep"
        )
        reportPath.write_text("15.0\t100\t100\tU\t0\tunclassified\n")

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ATAC-Seq"), self._config())

        mock_runPlusPF.assert_not_called()

    @patch("dissectBCL.postmux.runPlusPF")
    def test_skips_already_escalated_sample(self, mock_runPlusPF, tmp_path):
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        fqcSample = laneFolder / "FASTQC_Project_1_proj" / "Sample_S1"
        (fqcSample / "S1.rep").write_text("15.0\t100\t100\tU\t0\tunclassified\n")
        (fqcSample / "S1.plusPF.krakenreport").write_text("2.0\t100\t100\tU\t0\tunclassified\n")

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config())

        mock_runPlusPF.assert_not_called()

    @patch("dissectBCL.postmux.runPlusPF")
    @patch("dissectBCL.postmux.Pool", _SyncPool)
    @patch("dissectBCL.postmux.Popen")
    def test_fresh_run_writes_report_then_evaluates_it_for_escalation(
        self, mock_popen, mock_runPlusPF, tmp_path
    ):
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        reportPath = (
            laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep"
        )

        def _fake_kraken2(*args, **kwargs):
            # simulate kraken2 --report writing its output file
            reportPath.write_text("15.0\t100\t100\tU\t0\tunclassified\n")
            mock_wait = mock_popen.return_value.wait
            mock_wait.return_value = 0
            return mock_popen.return_value

        mock_popen.side_effect = _fake_kraken2

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config())

        assert reportPath.exists()
        mock_runPlusPF.assert_called_once()

    @patch("dissectBCL.postmux.runPlusPF")
    def test_missing_screening_section_skips_escalation_without_raising(
        self, mock_runPlusPF, tmp_path
    ):
        # Deployed dissectBCL.ini files that predate this feature won't
        # have a [screening] section -- the escalation check must degrade
        # to a no-op, not crash the whole postmux() run.
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        (laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep").write_text(
            "15.0\t100\t100\tU\t0\tunclassified\n"
        )
        config = configparser.ConfigParser()
        config["misc"] = {"threads": "10"}  # no [screening] section

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), config)

        mock_runPlusPF.assert_not_called()
