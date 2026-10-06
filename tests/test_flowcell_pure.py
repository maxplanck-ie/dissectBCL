import configparser
import datetime
import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from dissectBCL.flowcell import flowCellClass, sampleSheetClass

MINIMAL_RUNINFO_XML = """<?xml version="1.0"?>
<RunInfo Version="4">
  <Run Id="Flowcell" Number="1">
    <Flowcell>HHHHHHHHH</Flowcell>
    <Instrument>NB000000</Instrument>
    <Date>000000</Date>
    <Reads>
      <Read Number="1" NumCycles="150" IsIndexedRead="N" />
      <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="8" IsIndexedRead="Y" />
      <Read Number="4" NumCycles="150" IsIndexedRead="N" />
    </Reads>
    <FlowcellLayout LaneCount="4" SurfaceCount="2" SwathCount="3" TileCount="12" />
  </Run>
</RunInfo>
"""


class Test_parseRunInfo:
    def test_parses_reads_lanes_instrument_flowcell(self, tmp_path):
        runInfo = tmp_path / "RunInfo.xml"
        runInfo.write_text(MINIMAL_RUNINFO_XML)
        fake_self = SimpleNamespace(runInfo=runInfo)

        seqRecipe, lanes, instrument, flowcellID = flowCellClass.parseRunInfo(fake_self)

        assert seqRecipe == {
            "Read1": ["Y", 150],
            "Index1": ["I", 8],
            "Index2": ["I", 8],
            "Read2": ["Y", 150],
        }
        assert lanes == 4
        assert instrument == "NB000000"
        assert flowcellID == "HHHHHHHHH"


class Test_parseRunInfoAviti:
    def _write(self, tmp_path, cycles, lanes="1+2"):
        runInfo = tmp_path / "RunParameters.json"
        runInfo.write_text(
            json.dumps(
                {
                    "Cycles": cycles,
                    "AnalysisLanes": lanes,
                    "InstrumentName": "AV1",
                    "FlowcellID": "AVFC1",
                }
            )
        )
        return runInfo

    def test_parses_single_index_two_lanes(self, tmp_path):
        runInfo = self._write(tmp_path, {"R1": 150, "R2": 150, "I1": 8})
        fake_self = SimpleNamespace(runInfo=runInfo)

        seqRecipe, lanes, instrument, flowcellID = flowCellClass.parseRunInfoAviti(
            fake_self
        )

        assert seqRecipe == {
            "Read1": ["Y", 150],
            "Read2": ["Y", 150],
            "Index1": ["I", 8],
        }
        assert lanes == 2
        assert instrument == "AV1"
        assert flowcellID == "AVFC1"

    def test_dual_index_single_lane(self, tmp_path):
        runInfo = self._write(
            tmp_path, {"R1": 150, "R2": 150, "I1": 8, "I2": 8}, lanes="1"
        )
        fake_self = SimpleNamespace(runInfo=runInfo)

        seqRecipe, lanes, instrument, flowcellID = flowCellClass.parseRunInfoAviti(
            fake_self
        )

        assert seqRecipe["Index2"] == ["I", 8]
        assert lanes == 1


class Test_validateRunCompletion:
    def test_novaseq_assumes_success(self):
        fake_self = SimpleNamespace(sequencer="NovaSeq")

        assert flowCellClass.validateRunCompletion(fake_self) == "SuccessfullyCompleted"

    def test_miseq_reads_completion_status_xml(self, tmp_path):
        status = tmp_path / "RunCompletionStatus.xml"
        status.write_text(
            "<Run><CompletionStatus>SuccessfullyCompleted</CompletionStatus></Run>"
        )
        fake_self = SimpleNamespace(sequencer="MiSeq", runCompletionStatus=status)

        assert flowCellClass.validateRunCompletion(fake_self) == "SuccessfullyCompleted"

    def test_miseq_reports_failure_status(self, tmp_path):
        status = tmp_path / "RunCompletionStatus.xml"
        status.write_text("<Run><CompletionStatus>Failed</CompletionStatus></Run>")
        fake_self = SimpleNamespace(sequencer="MiSeq", runCompletionStatus=status)

        assert flowCellClass.validateRunCompletion(fake_self) == "Failed"


class Test_filesExist:
    def _fake_self(self, tmp_path, missing=None):
        paths = {
            "bclPath": tmp_path / "bcl",
            "origSS": tmp_path / "SampleSheet.csv",
            "runInfo": tmp_path / "RunInfo.xml",
            "inBaseDir": tmp_path / "in",
            "outBaseDir": tmp_path / "out",
        }
        for key, p in paths.items():
            if key == missing:
                continue
            if key in ("bclPath", "inBaseDir", "outBaseDir"):
                p.mkdir()
            else:
                p.write_text("x")
        return SimpleNamespace(name="FC1", config={}, **paths)

    def test_all_paths_exist_no_exit(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        flowCellClass.filesExist(fake_self)

    def test_missing_path_mails_home_and_exits(self, tmp_path):
        fake_self = self._fake_self(tmp_path, missing="runInfo")

        with (
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
            pytest.raises(SystemExit),
        ):
            flowCellClass.filesExist(fake_self)

        mock_mail.assert_called_once()


class Test_asdict:
    def test_returns_expected_keys_and_values(self, tmp_path):
        now = datetime.datetime(2026, 1, 1, 12, 0, 0)
        fake_self = SimpleNamespace(
            name="FC1",
            sequencer="NovaSeq",
            bclPath=tmp_path / "bcl",
            origSS=tmp_path / "SampleSheet.csv",
            runInfo=tmp_path / "RunInfo.xml",
            runCompletionStatus=tmp_path / "RunCompletionStatus.xml",
            successfulrun="SuccessfullyCompleted",
            inBaseDir=tmp_path / "in",
            outBaseDir=tmp_path / "out",
            logFile=tmp_path / "fc.log",
            seqRecipe={"Read1": ["Y", 150]},
            lanes=4,
            instrument="NB000000",
            flowcellID="HHHHHHHHH",
            startTime=now,
            config={"a": 1},
        )

        result = flowCellClass.asdict(fake_self)

        assert result["name"] == "FC1"
        assert result["flowcellID"] == "HHHHHHHHH"
        assert result["lanes"] == 4
        assert result["Time initiated"] == "01/01/2026, 12:00:00"
        assert result["config"] == {"a": 1}

    def test_missing_runCompletionStatus_defaults_empty_string(self, tmp_path):
        now = datetime.datetime(2026, 1, 1, 12, 0, 0)
        fake_self = SimpleNamespace(
            name="FC1",
            sequencer="aviti",
            bclPath=tmp_path / "bcl",
            origSS=tmp_path / "RunManifest.csv",
            runInfo=tmp_path / "RunParameters.json",
            successfulrun="SuccessfullyCompleted",
            inBaseDir=tmp_path / "in",
            outBaseDir=tmp_path / "out",
            logFile=tmp_path / "fc.log",
            seqRecipe={},
            lanes=2,
            instrument="AV1",
            flowcellID="FC1",
            startTime=now,
            config={},
        )

        result = flowCellClass.asdict(fake_self)

        assert result["runCompletionStatus"] == ""


class Test_decideSplit:
    def _fake_self(self, fullSS, runInfoLanes=2, forceLaneSplit=False):
        return SimpleNamespace(
            fullSS=fullSS, runInfoLanes=runInfoLanes, forceLaneSplit=forceLaneSplit
        )

    def test_splits_by_default_single_lane_samples(self):
        fullSS = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "Sample_Project": ["P1", "P2"],
                "Lane": ["1", "2"],
                "index": ["AAAA", "CCCC"],
                "index2": ["GGGG", "TTTT"],
            }
        )
        fake_self = self._fake_self(fullSS)

        assert sampleSheetClass.decideSplit(fake_self, aviti=False) is True

    def test_sample_on_multiple_lanes_disables_split(self):
        fullSS = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1"],
                "Sample_Project": ["P1", "P1"],
                "Lane": ["1", "2"],
                "index": ["AAAA", "AAAA"],
                "index2": ["GGGG", "GGGG"],
            }
        )
        fake_self = self._fake_self(fullSS)

        assert sampleSheetClass.decideSplit(fake_self, aviti=False) is False

    def test_project_on_multiple_lanes_disables_split(self):
        fullSS = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "Sample_Project": ["P1", "P1"],
                "Lane": ["1", "2"],
                "index": ["AAAA", "CCCC"],
                "index2": ["GGGG", "TTTT"],
            }
        )
        fake_self = self._fake_self(fullSS, runInfoLanes=2)

        assert sampleSheetClass.decideSplit(fake_self, aviti=False) is False

    def test_forceLaneSplit_overrides_to_true(self):
        fullSS = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1"],
                "Sample_Project": ["P1", "P1"],
                "Lane": ["1", "2"],
                "index": ["AAAA", "AAAA"],
                "index2": ["GGGG", "GGGG"],
            }
        )
        fake_self = self._fake_self(fullSS, forceLaneSplit=True)

        assert sampleSheetClass.decideSplit(fake_self, aviti=False) is True

    def test_index_clash_forces_split_back_on(self):
        # S1 sits on 2 lanes, which alone disables lane splitting. But
        # S2 and S3 share an identical index+index2 combo -- if lane
        # splitting stayed off, they'd be indistinguishable in the
        # combined sheet, so decideSplit must override back to True.
        fullSS = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1", "S2", "S3"],
                "Sample_Project": ["P1", "P1", "P1", "P1"],
                "Lane": ["1", "2", "1", "1"],
                "index": ["AAAA", "AAAA", "CCCC", "CCCC"],
                "index2": ["GGGG", "GGGG", "TTTT", "TTTT"],
            }
        )
        fake_self = self._fake_self(fullSS, runInfoLanes=2)

        assert sampleSheetClass.decideSplit(fake_self, aviti=False) is True

    def test_aviti_colnames_used_when_aviti_true(self):
        fullSS = pd.DataFrame(
            {
                "SampleName": ["S1", "S2"],
                "Project": ["P1", "P2"],
                "Lane": ["1", "2"],
                "Index1": ["AAAA", "CCCC"],
                "Index2": ["GGGG", "TTTT"],
            }
        )
        fake_self = self._fake_self(fullSS)

        assert sampleSheetClass.decideSplit(fake_self, aviti=True) is True


class Test_prepConvert:
    def _fake_self(self, ss, sequencer="NovaSeq"):
        ssDic = {"lane1": {"sampleSheet": ss}}
        return SimpleNamespace(
            sampleSheet=SimpleNamespace(ssDic=ssDic),
            seqRecipe={"Read1": ["Y", 150]},
            sequencer=sequencer,
            exitStats={},
        )

    def test_populates_mask_dualix_pe_convertopts_mismatch(self):
        ss = pd.DataFrame({"index": ["AAAAAAAA", "CCCCCCCC"]})
        fake_self = self._fake_self(ss)

        with (
            patch(
                "dissectBCL.flowcell.detMask",
                return_value=("maskstr", False, True, ["opt"], np.nan, np.nan),
            ) as mock_detmask,
            patch(
                "dissectBCL.flowcell.misMatcher", return_value="1,1"
            ) as mock_mismatch,
            patch("dissectBCL.flowcell.P5Seriesret", return_value="P5series"),
        ):
            flowCellClass.prepConvert(fake_self)

        ss_dict = fake_self.sampleSheet.ssDic["lane1"]
        assert ss_dict["mask"] == "maskstr"
        assert ss_dict["dualIx"] is False
        assert ss_dict["PE"] is True
        assert ss_dict["convertOpts"] == ["opt"]
        assert ss_dict["mismatch"] == "1,1"
        assert fake_self.exitStats["premux"] == 0
        mock_detmask.assert_called_once_with(
            fake_self.seqRecipe, ss, "lane1", "NovaSeq"
        )
        mock_mismatch.assert_called_once()

    def test_truncates_index_to_minP7_when_single_index(self):
        ss = pd.DataFrame({"index": ["AAAAAAAA", "CCCCCCCC"]})
        fake_self = self._fake_self(ss)

        with (
            patch(
                "dissectBCL.flowcell.detMask",
                return_value=("maskstr", False, True, [], np.nan, 4),
            ),
            patch("dissectBCL.flowcell.misMatcher", return_value="1"),
            patch("dissectBCL.flowcell.P5Seriesret", return_value="P5series"),
        ):
            flowCellClass.prepConvert(fake_self)

        truncated = fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"]["index"]
        assert list(truncated) == ["AAAA", "CCCC"]

    def test_truncates_both_indices_when_dualix(self):
        ss = pd.DataFrame(
            {"index": ["AAAAAAAA", "CCCCCCCC"], "index2": ["GGGGGGGG", "TTTTTTTT"]}
        )
        fake_self = self._fake_self(ss)

        with (
            patch(
                "dissectBCL.flowcell.detMask",
                return_value=("maskstr", True, True, [], 3, 5),
            ),
            patch("dissectBCL.flowcell.misMatcher", return_value="1,1"),
            patch("dissectBCL.flowcell.P5Seriesret", return_value="P5series"),
        ):
            flowCellClass.prepConvert(fake_self)

        updated = fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"]
        assert list(updated["index"]) == ["AAAAA", "CCCCC"]
        assert list(updated["index2"]) == ["GGG", "TTT"]

    def test_aviti_uses_Index1_Index2_colnames(self):
        ss = pd.DataFrame({"Index1": ["AAAAAAAA"], "Index2": ["GGGGGGGG"]})
        fake_self = self._fake_self(ss, sequencer="aviti")

        with (
            patch(
                "dissectBCL.flowcell.detMask",
                return_value=("maskstr", True, True, [], np.nan, np.nan),
            ) as mock_detmask,
            patch(
                "dissectBCL.flowcell.misMatcher", return_value="1,1"
            ) as mock_mismatch,
            patch("dissectBCL.flowcell.P5Seriesret", return_value="P5series"),
        ):
            flowCellClass.prepConvert(fake_self)

        mock_mismatch.assert_called_once()
        called_index_series = mock_mismatch.call_args[0][0]
        assert list(called_index_series) == ["AAAAAAAA"]
        mock_detmask.assert_called_once_with(fake_self.seqRecipe, ss, "lane1", "aviti")


class Test_organiseLogs:
    def _fake_self(self, tmp_path, ss):
        outBaseDir = tmp_path / "out"
        (outBaseDir / "lane1").mkdir(parents=True)
        ssDic = {"lane1": {"sampleSheet": ss, "lane": 1}}
        config = configparser.ConfigParser()
        config["software"] = {"bclconvert": "/bin/bclconvert"}
        fake_self = SimpleNamespace(
            sampleSheet=SimpleNamespace(ssDic=ssDic),
            outBaseDir=outBaseDir,
        )
        fake_self.asdict = lambda: {"name": "FC1", "config": config}
        return fake_self

    def test_writes_ssdf_and_yaml_and_config_files(self, tmp_path):
        ss = pd.DataFrame({"Sample_ID": ["S1"], "Sample_Project": ["P1"]})
        fake_self = self._fake_self(tmp_path, ss)

        flowCellClass.organiseLogs(fake_self)

        logDir = tmp_path / "out" / "lane1" / "Logs"
        assert (logDir / "sampleSheetdf.tsv").exists()
        assert (logDir / "outLaneInfo.yaml").exists()
        assert (logDir / "config.ini").exists()
        assert (logDir / "flowcellInfo.yaml").exists()

        ssdf_content = (logDir / "sampleSheetdf.tsv").read_text()
        assert "S1" in ssdf_content
        assert "P1" in ssdf_content

        yaml_content = (logDir / "outLaneInfo.yaml").read_text()
        assert "lane: 1" in yaml_content
        assert "sampleSheet" not in yaml_content

        config_content = (logDir / "config.ini").read_text()
        assert "bclconvert = /bin/bclconvert" in config_content

        flowcellinfo_content = (logDir / "flowcellInfo.yaml").read_text()
        assert "name: FC1" in flowcellinfo_content
        assert "config" not in flowcellinfo_content

    def test_removes_sampleSheet_key_from_ssDic_in_place(self, tmp_path):
        ss = pd.DataFrame({"Sample_ID": ["S1"], "Sample_Project": ["P1"]})
        fake_self = self._fake_self(tmp_path, ss)

        flowCellClass.organiseLogs(fake_self)

        assert "sampleSheet" not in fake_self.sampleSheet.ssDic["lane1"]


class Test_fakenews:
    def _fake_self(self, tmp_path):
        outBaseDir = tmp_path / "out"
        (outBaseDir / "lane1").mkdir(parents=True)
        return SimpleNamespace(
            sampleSheet=SimpleNamespace(ssDic={"lane1": {}}),
            outBaseDir=outBaseDir,
            config="theconfig",
            flowcellID="FC1",
            bclPath="/data/FC1",
            sequencer="NovaSeq",
            exitStats={},
        )

    def test_skips_lane_with_existing_communication_flag(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        (fake_self.outBaseDir / "lane1" / "communication.done").touch()

        with (
            patch("dissectBCL.flowcell.shipFiles") as mock_ship,
            patch("dissectBCL.flowcell.pushParkour") as mock_push,
            patch("dissectBCL.flowcell.gatherFinalMetrics") as mock_gather,
            patch("dissectBCL.flowcell.drHouseClass") as mock_house,
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
        ):
            flowCellClass.fakenews(fake_self)

        mock_ship.assert_not_called()
        mock_push.assert_not_called()
        mock_gather.assert_not_called()
        mock_house.assert_not_called()
        mock_mail.assert_not_called()

    def test_ships_pushes_mails_and_marks_done_on_success(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch("dissectBCL.flowcell.shipFiles", return_value={}) as mock_ship,
            patch("dissectBCL.flowcell.pushParkour", return_value=True) as mock_push,
            patch(
                "dissectBCL.flowcell.gatherFinalMetrics", return_value="metrics"
            ) as mock_gather,
            patch("dissectBCL.flowcell.drHouseClass") as mock_house,
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
        ):
            mock_house.return_value.prepMail.return_value = ("subj", "html")
            flowCellClass.fakenews(fake_self)

        mock_ship.assert_called_once_with(
            fake_self.outBaseDir / "lane1", fake_self.config
        )
        mock_push.assert_called_once_with(
            fake_self.flowcellID,
            fake_self.sampleSheet,
            fake_self.config,
            fake_self.bclPath,
            fake_self.sequencer,
            outBaseDir=fake_self.outBaseDir,
        )
        mock_gather.assert_called_once_with("lane1", fake_self)
        mock_house.assert_called_once_with("metrics")
        mock_mail.assert_called_once_with("subj", "html", fake_self.config)
        assert fake_self.exitStats["lane1"]["pushParkour"] is True
        assert (fake_self.outBaseDir / "lane1" / "communication.done").exists()

    def test_does_not_mark_done_when_shipping_failed(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch(
                "dissectBCL.flowcell.shipFiles",
                return_value={"failedProjects": ["P1"]},
            ),
            patch("dissectBCL.flowcell.pushParkour", return_value=True),
            patch("dissectBCL.flowcell.gatherFinalMetrics", return_value="metrics"),
            patch("dissectBCL.flowcell.drHouseClass") as mock_house,
            patch("dissectBCL.flowcell.mailHome"),
        ):
            mock_house.return_value.prepMail.return_value = ("subj", "html")
            flowCellClass.fakenews(fake_self)

        assert not (fake_self.outBaseDir / "lane1" / "communication.done").exists()


class Test_postmux:
    def _fake_self(self, tmp_path, sequencer="NovaSeq"):
        outBaseDir = tmp_path / "out"
        (outBaseDir / "lane1").mkdir(parents=True)
        ss = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "Sample_Project": ["P1", "P1"],
            }
        )
        return SimpleNamespace(
            sampleSheet=SimpleNamespace(
                ssDic={"lane1": {"sampleSheet": ss, "PE": True}},
                laneSplitStatus=True,
            ),
            outBaseDir=outBaseDir,
            config="theconfig",
            sequencer=sequencer,
            exitStats={},
        )

    def test_illumina_runs_full_pipeline_and_touches_flags(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        laneFolder = fake_self.outBaseDir / "lane1"

        with (
            patch("dissectBCL.flowcell.renameProject") as mock_rename,
            patch("dissectBCL.flowcell.validateFqEnds") as mock_validate,
            patch("dissectBCL.flowcell.qcs") as mock_qcs,
            patch("dissectBCL.flowcell.clumper") as mock_clump,
            patch("dissectBCL.flowcell.kraken") as mock_kraken,
            patch("dissectBCL.flowcell.md5_multiqc") as mock_multiqc,
            patch("dissectBCL.flowcell.moveOptDup") as mock_movedup,
        ):
            flowCellClass.postmux(fake_self)

        ss = fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"]
        mock_rename.assert_called_once_with(
            laneFolder / "P1", ss, fake_self.sampleSheet.laneSplitStatus
        )
        mock_validate.assert_called_once_with(laneFolder / "P1", fake_self)
        mock_qcs.assert_called_once_with(
            "P1", laneFolder, {"S1", "S2"}, fake_self.config
        )
        mock_clump.assert_called_once_with(
            "P1", laneFolder, {"S1", "S2"}, fake_self.config, True, "NovaSeq"
        )
        mock_kraken.assert_called_once_with(
            "P1", laneFolder, {"S1", "S2"}, ss, fake_self.config
        )
        mock_multiqc.assert_called_once_with("P1", laneFolder, fake_self)
        mock_movedup.assert_called_once_with(laneFolder)

        assert (laneFolder / ".P1.renamed.done").exists()
        assert (laneFolder / ".P1.postmux.done").exists()
        assert fake_self.exitStats["postmux"] == 0

    def test_aviti_renames_under_samples_subdir(self, tmp_path):
        fake_self = self._fake_self(tmp_path, sequencer="aviti")
        laneFolder = fake_self.outBaseDir / "lane1"

        with (
            patch("dissectBCL.flowcell.renameProject") as mock_rename,
            patch("dissectBCL.flowcell.validateFqEnds"),
            patch("dissectBCL.flowcell.qcs"),
            patch("dissectBCL.flowcell.clumper"),
            patch("dissectBCL.flowcell.kraken"),
            patch("dissectBCL.flowcell.md5_multiqc"),
            patch("dissectBCL.flowcell.moveOptDup"),
        ):
            flowCellClass.postmux(fake_self)

        ss = fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"]
        mock_rename.assert_called_once_with(
            laneFolder / "Samples" / "P1", ss, fake_self.sampleSheet.laneSplitStatus
        )

    def test_skips_rename_and_postmux_when_flags_already_present(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        laneFolder = fake_self.outBaseDir / "lane1"
        (laneFolder / ".P1.renamed.done").touch()
        (laneFolder / ".P1.postmux.done").touch()

        with (
            patch("dissectBCL.flowcell.renameProject") as mock_rename,
            patch("dissectBCL.flowcell.validateFqEnds") as mock_validate,
            patch("dissectBCL.flowcell.qcs") as mock_qcs,
            patch("dissectBCL.flowcell.clumper") as mock_clump,
            patch("dissectBCL.flowcell.kraken") as mock_kraken,
            patch("dissectBCL.flowcell.md5_multiqc") as mock_multiqc,
            patch("dissectBCL.flowcell.moveOptDup") as mock_movedup,
        ):
            flowCellClass.postmux(fake_self)

        mock_rename.assert_not_called()
        mock_validate.assert_not_called()
        mock_qcs.assert_not_called()
        mock_clump.assert_not_called()
        mock_kraken.assert_not_called()
        mock_multiqc.assert_not_called()
        mock_movedup.assert_not_called()
        assert fake_self.exitStats["postmux"] == 0


class Test_demux:
    def _fake_self(
        self, tmp_path, sequencer="NovaSeq", successfulrun="SuccessfullyCompleted"
    ):
        outBaseDir = tmp_path / "out"
        outBaseDir.mkdir()
        ss = pd.DataFrame({"Sample_ID": ["S1"], "Sample_Project": ["P1"]})
        return SimpleNamespace(
            successfulrun=successfulrun,
            sampleSheet=SimpleNamespace(
                ssDic={"lane1": {"sampleSheet": ss, "dualIx": False}},
                laneSplitStatus=True,
            ),
            outBaseDir=outBaseDir,
            config="theconfig",
            bclconvert_path="/bin/bclconvert",
            bclPath="/data/FC1",
            num_threads=4,
            sequencer=sequencer,
            name="FC1",
            exitStats={},
        )

    def _popen(self, returncode=0, stderr=b""):
        proc = patch("dissectBCL.flowcell.Popen").start()
        proc.return_value.communicate.return_value = (b"", stderr)
        proc.return_value.returncode = returncode
        return proc

    def test_failed_run_marks_lanes_failed_and_mails_no_demux(self, tmp_path):
        fake_self = self._fake_self(tmp_path, successfulrun="Failed")

        with (
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
            patch("dissectBCL.flowcell.writeDemuxSheet") as mock_write,
            patch("dissectBCL.flowcell.Popen") as mock_popen,
        ):
            flowCellClass.demux(fake_self)

        assert (fake_self.outBaseDir / "lane1" / "run.failed").exists()
        mock_mail.assert_called_once()
        assert mock_mail.call_args.kwargs["toCore"] is True
        mock_write.assert_not_called()
        mock_popen.assert_not_called()

    def test_success_writes_sheet_runs_bclconvert_and_parses_stats(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch("dissectBCL.flowcell.writeDemuxSheet") as mock_write,
            patch("dissectBCL.flowcell.compareDemuxSheet") as mock_compare,
            patch(
                "dissectBCL.flowcell.parseStats", return_value="parsed"
            ) as mock_parse,
        ):
            self._popen(returncode=0)
            try:
                flowCellClass.demux(fake_self)
            finally:
                patch.stopall()

        outputFolder = fake_self.outBaseDir / "lane1"
        mock_write.assert_called_once_with(
            outputFolder / "demuxSheet.csv",
            fake_self.sampleSheet.ssDic["lane1"],
            True,
        )
        mock_compare.assert_not_called()
        assert (outputFolder / "bclconvert.done").exists()
        mock_parse.assert_called_once()
        assert fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"] == "parsed"
        assert fake_self.sampleSheet.ssDic["lane1"]["P5RC"] is False
        assert fake_self.exitStats["demux"] == 0

    def test_existing_demux_sheet_compares_instead_of_writing(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        outputFolder = fake_self.outBaseDir / "lane1"
        outputFolder.mkdir()
        (outputFolder / "demuxSheet.csv").write_text("x")

        with (
            patch("dissectBCL.flowcell.writeDemuxSheet") as mock_write,
            patch("dissectBCL.flowcell.compareDemuxSheet") as mock_compare,
            patch("dissectBCL.flowcell.parseStats", return_value="parsed"),
        ):
            self._popen(returncode=0)
            try:
                flowCellClass.demux(fake_self)
            finally:
                patch.stopall()

        mock_write.assert_not_called()
        mock_compare.assert_called_once_with(
            fake_self.sampleSheet.ssDic["lane1"], outputFolder / "demuxSheet.csv"
        )

    def test_existing_bclconvert_done_skips_popen(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        outputFolder = fake_self.outBaseDir / "lane1"
        outputFolder.mkdir()
        (outputFolder / "bclconvert.done").touch()

        with (
            patch("dissectBCL.flowcell.writeDemuxSheet"),
            patch("dissectBCL.flowcell.Popen") as mock_popen,
            patch("dissectBCL.flowcell.parseStats", return_value="parsed"),
        ):
            flowCellClass.demux(fake_self)

        mock_popen.assert_not_called()
        assert fake_self.exitStats["demux"] == 0

    def test_nonzero_exitcode_mails_and_exits(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch("dissectBCL.flowcell.writeDemuxSheet"),
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
        ):
            self._popen(returncode=1, stderr=b"boom")
            try:
                with pytest.raises(SystemExit):
                    flowCellClass.demux(fake_self)
            finally:
                patch.stopall()

        mock_mail.assert_called_once()
        assert "boom" in mock_mail.call_args.args[1]
        outputFolder = fake_self.outBaseDir / "lane1"
        assert not (outputFolder / "bclconvert.done").exists()

    def test_miseq_p5rc_triggers_rerun_and_matching_sheets(self, tmp_path):
        fake_self = self._fake_self(tmp_path, sequencer="MiSeq")
        outputFolder = fake_self.outBaseDir / "lane1"
        origSheet = fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"]

        with (
            patch("dissectBCL.flowcell.writeDemuxSheet"),
            patch("dissectBCL.flowcell.evalMiSeqP5", return_value=True),
            patch("dissectBCL.flowcell.readDemuxSheet", return_value="demuxdf"),
            patch(
                "dissectBCL.flowcell.matchingSheets", return_value="matched"
            ) as mock_match,
            patch("dissectBCL.flowcell.parseStats", return_value="parsed"),
        ):
            proc = patch("dissectBCL.flowcell.Popen").start()
            proc.return_value.communicate.return_value = (b"", b"")
            proc.return_value.returncode = 0
            outputFolder.mkdir()

            def side_effect(*args, **kwargs):
                (outputFolder / "Reports").mkdir(exist_ok=True)
                (outputFolder / "Logs").mkdir(exist_ok=True)
                return proc.return_value

            proc.side_effect = side_effect
            try:
                flowCellClass.demux(fake_self)
            finally:
                patch.stopall()

        mock_match.assert_called_once()
        called_args = mock_match.call_args.args
        assert called_args[0] is origSheet
        assert called_args[1] == "demuxdf"
        assert fake_self.sampleSheet.ssDic["lane1"]["P5RC"] is True
        assert proc.call_count == 2


class Test_demux_aviti:
    def _fake_self(self, tmp_path, successfulrun="SuccessfullyCompleted"):
        outBaseDir = tmp_path / "out"
        outBaseDir.mkdir()
        ss = pd.DataFrame({"Sample_ID": ["S1"], "Sample_Project": ["P1"]})
        return SimpleNamespace(
            successfulrun=successfulrun,
            sampleSheet=SimpleNamespace(
                ssDic={"lane1": {"sampleSheet": ss}},
                laneSplitStatus=True,
            ),
            outBaseDir=outBaseDir,
            config="theconfig",
            bases2fastq_path="/bin/bases2fastq",
            bclPath="/data/FC1",
            num_threads=4,
            name="FC1",
            exitStats={},
        )

    def _popen(self, returncode=0, stderr=b""):
        proc = patch("dissectBCL.flowcell.Popen").start()
        proc.return_value.communicate.return_value = (b"", stderr)
        proc.return_value.returncode = returncode
        return proc

    def test_failed_run_marks_lanes_failed_and_mails_no_demux(self, tmp_path):
        fake_self = self._fake_self(tmp_path, successfulrun="Failed")

        with patch("dissectBCL.flowcell.mailHome") as mock_mail:
            flowCellClass.demux_aviti(fake_self)

        assert (fake_self.outBaseDir / "lane1" / "run.failed").exists()
        mock_mail.assert_called_once()
        assert mock_mail.call_args.kwargs["toCore"] is True
        assert "demux" not in fake_self.exitStats

    def test_success_writes_manifest_runs_bases2fastq_and_parses_stats(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch("dissectBCL.flowcell.writeDemuxSheetAviti") as mock_write,
            patch(
                "dissectBCL.flowcell.parseStats", return_value="parsed"
            ) as mock_parse,
        ):
            self._popen(returncode=0)
            try:
                flowCellClass.demux_aviti(fake_self)
            finally:
                patch.stopall()

        outputFolder = fake_self.outBaseDir / "lane1"
        mock_write.assert_called_once_with(
            outputFolder / "manifest" / "RunManifest.csv",
            fake_self.sampleSheet.ssDic["lane1"],
            True,
        )
        assert (outputFolder / "bases2fastq.done").exists()
        mock_parse.assert_called_once()
        assert mock_parse.call_args.kwargs.get("mode") == "aviti"
        assert fake_self.sampleSheet.ssDic["lane1"]["P5RC"] is False
        assert fake_self.sampleSheet.ssDic["lane1"]["sampleSheet"] == "parsed"
        assert fake_self.exitStats["demux"] == 0

    def test_existing_manifest_skips_write(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        outputFolder = fake_self.outBaseDir / "lane1"
        (outputFolder / "manifest").mkdir(parents=True)
        (outputFolder / "manifest" / "RunManifest.csv").touch()

        with (
            patch("dissectBCL.flowcell.writeDemuxSheetAviti") as mock_write,
            patch("dissectBCL.flowcell.parseStats", return_value="parsed"),
        ):
            self._popen(returncode=0)
            try:
                flowCellClass.demux_aviti(fake_self)
            finally:
                patch.stopall()

        mock_write.assert_not_called()

    def test_existing_bases2fastq_done_skips_popen(self, tmp_path):
        fake_self = self._fake_self(tmp_path)
        outputFolder = fake_self.outBaseDir / "lane1"
        outputFolder.mkdir()
        (outputFolder / "bases2fastq.done").touch()

        with (
            patch("dissectBCL.flowcell.writeDemuxSheetAviti"),
            patch("dissectBCL.flowcell.Popen") as mock_popen,
            patch("dissectBCL.flowcell.parseStats", return_value="parsed"),
        ):
            flowCellClass.demux_aviti(fake_self)

        mock_popen.assert_not_called()
        assert fake_self.exitStats["demux"] == 0

    def test_nonzero_exitcode_mails_and_exits(self, tmp_path):
        fake_self = self._fake_self(tmp_path)

        with (
            patch("dissectBCL.flowcell.writeDemuxSheetAviti"),
            patch("dissectBCL.flowcell.mailHome") as mock_mail,
        ):
            self._popen(returncode=1, stderr=b"boom")
            try:
                with pytest.raises(SystemExit):
                    flowCellClass.demux_aviti(fake_self)
            finally:
                patch.stopall()

        mock_mail.assert_called_once()
        assert "boom" in mock_mail.call_args.args[1]
        outputFolder = fake_self.outBaseDir / "lane1"
        assert not (outputFolder / "bases2fastq.done").exists()
