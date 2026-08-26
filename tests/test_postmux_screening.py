import configparser
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from dissectBCL.fakeNews import buildContaminationDic
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
    def _config(self, tmp_path):
        # plusPFdb must exist on disk (kraken() now bails out early if it
        # doesn't -- see the guard added for finding #2), so point it at a
        # real (if empty) directory rather than a literal fake path.
        plusPFdb = tmp_path / "pluspf"
        plusPFdb.mkdir(exist_ok=True)
        c = configparser.ConfigParser()
        c["misc"] = {"threads": "10"}
        c["software"] = {"kraken2db": "/fake/krakendb"}
        c["screening"] = {
            "plusPFdb": str(plusPFdb),
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

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config(tmp_path))

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

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ATAC-Seq"), self._config(tmp_path))

        mock_runPlusPF.assert_not_called()

    @patch("dissectBCL.postmux.runPlusPF")
    def test_skips_already_escalated_sample(self, mock_runPlusPF, tmp_path):
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        fqcSample = laneFolder / "FASTQC_Project_1_proj" / "Sample_S1"
        (fqcSample / "S1.rep").write_text("15.0\t100\t100\tU\t0\tunclassified\n")
        (fqcSample / "S1.plusPF.krakenreport").write_text("2.0\t100\t100\tU\t0\tunclassified\n")

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config(tmp_path))

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

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config(tmp_path))

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

    @patch("dissectBCL.postmux.runPlusPF")
    def test_missing_library_type_column_skips_escalation_without_raising(
        self, mock_runPlusPF, tmp_path
    ):
        # ssdf can lack a Library_Type column entirely (the parkourDF.empty
        # path -- see flowcell.py around lines 790/828). The escalation
        # block must degrade gracefully rather than raise KeyError.
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        (laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep").write_text(
            "15.0\t100\t100\tU\t0\tunclassified\n"
        )
        ssdf = pd.DataFrame({"Sample_ID": ["S1"]})  # no Library_Type column

        kraken("1_proj", laneFolder, ["S1"], ssdf, self._config(tmp_path))

        mock_runPlusPF.assert_not_called()

    @patch("dissectBCL.postmux.runPlusPF")
    def test_sample_folder_with_no_fastqs_skips_escalation_without_raising(
        self, mock_runPlusPF, tmp_path
    ):
        # krakenfqs() indexes into an empty fastq list (IndexError) when a
        # sample folder has zero matching fastq files, rather than
        # returning None like it does for the "too many fastqs" case. The
        # escalation loop must treat both the same way: skip, don't crash.
        laneFolder = tmp_path / "lane"
        sampleFolder = laneFolder / "Project_1_proj" / "Sample_S1"
        sampleFolder.mkdir(parents=True)  # no fastq.gz files inside
        fqcSample = laneFolder / "FASTQC_Project_1_proj" / "Sample_S1"
        fqcSample.mkdir(parents=True)
        (fqcSample / "S1.rep").write_text("15.0\t100\t100\tU\t0\tunclassified\n")

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), self._config(tmp_path))

        mock_runPlusPF.assert_not_called()

    @patch("dissectBCL.postmux.runPlusPF")
    def test_missing_plusPFdb_skips_escalation_without_raising(
        self, mock_runPlusPF, tmp_path
    ):
        # A [screening] section with plusPFdb missing/unset, or pointing at
        # a path that doesn't exist (e.g. the shipped template's literal
        # placeholder), must not fire a doomed kraken2 run.
        laneFolder = tmp_path / "lane"
        _make_sample(laneFolder, "1_proj", "S1")
        (laneFolder / "FASTQC_Project_1_proj" / "Sample_S1" / "S1.rep").write_text(
            "15.0\t100\t100\tU\t0\tunclassified\n"
        )
        config = configparser.ConfigParser()
        config["misc"] = {"threads": "10"}
        config["screening"] = {
            "plusPFdb": "/path/to/kraken2_contaminome/pluspf",
            "unclassified_threshold": "10",
        }

        kraken("1_proj", laneFolder, ["S1"], self._ssdf("S1", "ChIP-Seq"), config)

        mock_runPlusPF.assert_not_called()


class Test_runPlusPF_buildContaminationDic_handoff:
    """
    Integration test locking the filename handoff between runPlusPF()
    (writer) and buildContaminationDic() (reader): both currently agree on
    "<sample>.plusPF.krakenreport" only via a hardcoded literal in each
    side's own tests -- nothing exercises the two functions together.
    """

    @patch("dissectBCL.postmux.Pool", _SyncPool)
    @patch("dissectBCL.postmux.Popen")
    def test_runPlusPF_report_is_picked_up_by_buildContaminationDic(
        self, mock_popen, tmp_path
    ):
        laneFolder = tmp_path / "lane"
        sampleFolder = _make_sample(laneFolder, "1_proj", "S1")
        # Primary (routine) kraken report -- buildContaminationDic reads
        # this too, and expects it to already exist.
        fqcSampleDir = laneFolder / "FASTQC_Project_1_proj" / "Sample_S1"
        (fqcSampleDir / "S1.rep").write_text(
            "94.0\t940\t940\tU\t0\tunclassified\n6.0\t60\t60\tS\t10090\tmouse\n"
        )

        def _fake_kraken2(cmd, *args, **kwargs):
            # Locate the --report path from the real call args (as built by
            # runPlusPF) and write a minimal kraken2 report to it, exactly
            # as kraken2 itself would.
            reportPath = Path(cmd[cmd.index("--report") + 1])
            reportPath.write_text(
                "5.0\t50\t50\tU\t0\tunclassified\n"
                "95.0\t950\t950\tS\t3702\tarabidopsis\n"
            )
            mock_popen.return_value.wait.return_value = 0
            return mock_popen.return_value

        mock_popen.side_effect = _fake_kraken2

        config = configparser.ConfigParser()
        config["misc"] = {"threads": "10"}
        config["screening"] = {"plusPFdb": str(tmp_path / "pluspf")}

        runPlusPF("1_proj", laneFolder, ["S1"], config)

        assert (fqcSampleDir / "S1.plusPF.krakenreport").exists()

        ssdf = pd.DataFrame(
            {"Sample_ID": ["S1"], "Organism": [["mouse (GRCm39)"]]}
        )
        result = buildContaminationDic(laneFolder, ssdf)

        assert result["S1"][3] == "arabidopsis"
        assert sampleFolder.exists()  # sanity: fixture actually made the sample
