import datetime
import json
from types import SimpleNamespace
from unittest.mock import patch

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
