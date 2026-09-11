import os
from pathlib import Path

import numpy as np
import pandas as pd

from dissectBCL.demux import (
    compareDemuxSheet,
    matchingSheets,
    misMatcher,
    parseStats,
    readDemuxSheet,
)


class Test_misMatcher:
    def test_empty_series_return_empty_dict(self):
        assert misMatcher(pd.Series([], dtype=object), pd.Series([], dtype=object), "illumina") == {}

    def test_all_null_series_return_empty_dict(self):
        P7s = pd.Series([np.nan, np.nan])
        P5s = pd.Series([np.nan, np.nan])
        assert misMatcher(P7s, P5s, "illumina") == {}

    def test_illumina_keys_and_mismatch_from_hamming_distance(self):
        # "AAAAA" vs "AAAAT" differ in 1 position -> hamming2Mismatch(1) == 0
        P7s = pd.Series(["AAAAA", "AAAAT"])
        # "AAAAA" vs "TTTTT" differ in 5 positions -> hamming2Mismatch(5) == 2
        P5s = pd.Series(["AAAAA", "TTTTT"])

        result = misMatcher(P7s, P5s, "illumina")

        assert result == {
            "BarcodeMismatchesIndex1": 0,
            "BarcodeMismatchesIndex2": 2,
        }

    def test_aviti_uses_mismatch_threshold_keys(self):
        P7s = pd.Series(["AAAAA", "TTTTT"])
        P5s = pd.Series([], dtype=object)

        result = misMatcher(P7s, P5s, "aviti")

        assert result == {"I1MismatchThreshold": 2}


class Test_matchingSheets:
    def _autodf(self):
        return pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "index": ["AAAA", "CCCC"],
                "index2": ["GGGG", "TTTT"],
                "I7_Index_ID": ["i7-A", "i7-C"],
                "I5_Index_ID": ["i5-G", "i5-T"],
            }
        )

    def test_no_index_column_returns_autodf_unchanged(self):
        autodf = pd.DataFrame({"Sample_ID": ["S1"]})
        mandf = pd.DataFrame({"Sample_ID": ["S1"], "index": ["AAAA"]})

        result = matchingSheets(autodf, mandf)

        assert result is autodf

    def test_dual_index_overwrite_clears_index_id_columns(self):
        autodf = self._autodf()
        mandf = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "index": ["AAAA", "NEWC"],
                "index2": ["GGGG", "NEWT"],
            }
        )

        result = matchingSheets(autodf, mandf)

        s2 = result[result["Sample_ID"] == "S2"].iloc[0]
        assert s2["index"] == "NEWC"
        assert s2["index2"] == "NEWT"
        assert pd.isna(s2["I7_Index_ID"])
        assert pd.isna(s2["I5_Index_ID"])
        # untouched sample keeps its original index-derived columns.
        s1 = result[result["Sample_ID"] == "S1"].iloc[0]
        assert s1["I7_Index_ID"] == "i7-A"

    def test_single_index_overwrite_nans_index2(self):
        autodf = self._autodf()
        mandf = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "index": ["AAAA", "NEWC"],
            }
        )

        result = matchingSheets(autodf, mandf)

        s2 = result[result["Sample_ID"] == "S2"].iloc[0]
        assert s2["index"] == "NEWC"
        assert pd.isna(s2["index2"])
        assert pd.isna(s2["I5_Index_ID"])


class Test_compareDemuxSheet:
    def test_updates_ssDic_from_manually_edited_demuxSheet(self, tmp_path):
        fixture = Path(
            os.path.dirname(os.path.realpath(__file__)),
            "test_demux",
            "demuxSheet.csv",
        )
        # demuxSheet.csv defines mask 'Y101;I8N2;I8N16;Y101',
        # BarcodeMismatchesIndex1/2 == 1, and Sample_ID '22L002048' with
        # index CTCGAATA / index2 TCGACATC. matchingSheets requires every
        # mandf Sample_ID to already exist in autodf, so start from the
        # fixture's own data and only perturb the one sample under test.
        _, man_df, _, _ = readDemuxSheet(fixture)
        ssdf = man_df.copy()
        ssdf["I7_Index_ID"] = "orig-i7"
        ssdf["I5_Index_ID"] = "orig-i5"
        ssdf.loc[ssdf["Sample_ID"] == "22L002048", "index"] = "OLDOLDOL"
        ssdf.loc[ssdf["Sample_ID"] == "22L002048", "index2"] = "OLDOLDOL"
        ssDic = {
            "mismatch": {"BarcodeMismatchesIndex1": 2, "BarcodeMismatchesIndex2": 2},
            "mask": "Y100;I7;I7;Y100",
            "dualIx": False,
            "sampleSheet": ssdf,
        }

        compareDemuxSheet(ssDic, fixture)

        assert ssDic["mismatch"] == {
            "BarcodeMismatchesIndex1": 1,
            "BarcodeMismatchesIndex2": 1,
        }
        assert ssDic["mask"] == "Y101;I8N2;I8N16;Y101"
        assert ssDic["dualIx"] is True
        assert ssDic["P5RC"] is False
        updated = ssDic["sampleSheet"]
        row = updated[updated["Sample_ID"] == "22L002048"].iloc[0]
        assert row["index"] == "CTCGAATA"
        assert row["index2"] == "TCGACATC"


class Test_parseStats:
    def test_illumina_mode_computes_qc_and_depth_per_sample(self, tmp_path):
        outputFolder = tmp_path / "outLane"
        reportsDir = outputFolder / "Reports"
        reportsDir.mkdir(parents=True)
        pd.DataFrame(
            {
                "SampleID": ["S1", "S1", "Undetermined"],
                "ReadNumber": [1, 2, 1],
                "Mean Quality Score (PF)": [35.0, 34.0, 20.0],
                "% Q30": [0.9, 0.85, 0.5],
            }
        ).to_csv(reportsDir / "Quality_Metrics.csv", index=False)
        pd.DataFrame(
            {
                "SampleID": ["S1", "Undetermined"],
                "# Reads": [1000, 50],
            }
        ).to_csv(reportsDir / "Demultiplex_Stats.csv", index=False)
        ssdf = pd.DataFrame({"Sample_ID": ["S1"], "Sample_Project": ["P1"]})

        result = parseStats(outputFolder, ssdf, mode="illumina")

        assert list(result["Sample_ID"]) == ["S1"]
        row = result.iloc[0]
        assert row["meanQ"] == "1:35.0,2:34.0"
        assert row["percQ30"] == "1:90,2:85"
        assert row["gotDepth"] == 1000

    def test_unsupported_mode_logs_and_returns_none(self, tmp_path):
        ssdf = pd.DataFrame({"Sample_ID": ["S1"]})
        assert parseStats(tmp_path, ssdf, mode="nope") is None
