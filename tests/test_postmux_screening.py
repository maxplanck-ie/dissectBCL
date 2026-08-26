import configparser
from unittest.mock import patch

from dissectBCL.postmux import runPlusPF


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
