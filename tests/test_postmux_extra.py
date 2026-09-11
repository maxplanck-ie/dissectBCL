from unittest.mock import patch

import pandas as pd
import pytest

from dissectBCL.postmux import (
    matchIDtoName,
    renamefq,
    renameProject,
    validateFqEnds,
)


def _ssdf():
    return pd.DataFrame(
        {
            "Sample_ID": ["S1", "S2"],
            "Sample_Name": ["SampleOne", "SampleTwo"],
        }
    )


class Test_matchIDtoName:
    def test_returns_name_for_known_id(self):
        assert matchIDtoName("S1", _ssdf()) == "SampleOne"

    def test_unknown_id_exits(self):
        with pytest.raises(SystemExit):
            matchIDtoName("UNKNOWN", _ssdf())

    def test_id_split_across_lanes_with_consistent_name_returns_name(self):
        ssdf = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1"],
                "Sample_Name": ["SampleOne", "SampleOne"],
            }
        )
        assert matchIDtoName("S1", ssdf) == "SampleOne"

    def test_id_split_across_lanes_with_conflicting_names_exits(self):
        ssdf = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1"],
                "Sample_Name": ["SampleOne", "SampleOneDifferent"],
            }
        )
        with pytest.raises(SystemExit):
            matchIDtoName("S1", ssdf)

    def test_nan_name_exits(self):
        ssdf = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S1"],
                "Sample_Name": [None, None],
            }
        )
        with pytest.raises(SystemExit):
            matchIDtoName("S1", ssdf)


class Test_renamefq:
    def test_lane_split_strips_lane_and_setnum(self, tmp_path):
        projectFolder = tmp_path / "1906_Hein_B03_Hein"
        projectFolder.mkdir()
        fqFile = projectFolder / "S1_S63_L001_R2_001.fastq.gz"
        fqFile.write_bytes(b"")

        result = renamefq(fqFile, projectFolder, _ssdf(), laneSplitStatus=True)

        assert result == projectFolder / "Sample_S1" / "SampleOne_R2.fastq.gz"

    def test_no_lane_split_strips_setnum_only(self, tmp_path):
        projectFolder = tmp_path / "1906_Hein_B03_Hein"
        projectFolder.mkdir()
        fqFile = projectFolder / "S1_S63_R2_001.fastq.gz"
        fqFile.write_bytes(b"")

        result = renamefq(fqFile, projectFolder, _ssdf(), laneSplitStatus=False)

        assert result == projectFolder / "Sample_S1" / "SampleOne_R2.fastq.gz"

    def test_moves_matching_stats_json_into_sample_dir(self, tmp_path):
        projectFolder = tmp_path / "1906_Hein_B03_Hein"
        projectFolder.mkdir()
        fqFile = projectFolder / "S1_S63_R2_001.fastq.gz"
        fqFile.write_bytes(b"")
        (projectFolder / "S1_stats.json").write_text("{}")

        renamefq(fqFile, projectFolder, _ssdf(), laneSplitStatus=False)

        assert (projectFolder / "Sample_S1" / "S1_stats.json").exists()
        assert not (projectFolder / "S1_stats.json").exists()


class Test_renameProject:
    def test_moves_fastqs_and_renames_project_folder(self, tmp_path):
        projectFolder = tmp_path / "1906_Hein_B03_Hein"
        projectFolder.mkdir()
        (projectFolder / "S1_S63_R2_001.fastq.gz").write_bytes(b"")

        renameProject(projectFolder, _ssdf(), laneSplitStatus=False)

        renamed = tmp_path / "Project_1906_Hein_B03_Hein"
        assert renamed.exists()
        assert (renamed / "Sample_S1" / "SampleOne_R2.fastq.gz").exists()
        assert not projectFolder.exists()


class Test_validateFqEnds:
    class _FakeFlowcell:
        name = "FC123"
        config = {}

    def test_all_proper_endings_do_not_raise(self, tmp_path):
        (tmp_path / "S1_R1.fastq.gz").write_bytes(b"")
        (tmp_path / "S1_R2.fastq.gz").write_bytes(b"")

        validateFqEnds(tmp_path, self._FakeFlowcell())

    @patch("dissectBCL.postmux.mailHome")
    def test_malformed_ending_mails_and_exits(self, mock_mailHome, tmp_path):
        (tmp_path / "S1_weird.fastq.gz").write_bytes(b"")

        with pytest.raises(SystemExit):
            validateFqEnds(tmp_path, self._FakeFlowcell())
        mock_mailHome.assert_called_once()

    def test_undetermined_is_ignored(self, tmp_path):
        (tmp_path / "Undetermined_weird.fastq.gz").write_bytes(b"")

        validateFqEnds(tmp_path, self._FakeFlowcell())
