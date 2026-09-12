import configparser

import numpy as np
import pandas as pd
import pytest

from dissectBCL.misc import fetchLatestSeqDir, matchOptdupsReqs


def _config(piDir, seqDir="sequencing_data"):
    config = configparser.ConfigParser()
    config["Dirs"] = {"piDir": str(piDir)}
    config["Internals"] = {"seqDir": seqDir}
    return config


class Test_fetchLatestSeqDir:
    def test_no_suffixed_dir_returns_bare_sequencing_data(self, tmp_path):
        (tmp_path / "goodpi" / "sequencing_data").mkdir(parents=True)

        result = fetchLatestSeqDir(_config(tmp_path), "goodpi")

        assert result == tmp_path / "goodpi" / "sequencing_data"

    def test_picks_highest_numbered_suffix(self, tmp_path):
        (tmp_path / "goodpi" / "sequencing_data").mkdir(parents=True)
        (tmp_path / "goodpi" / "sequencing_data2").mkdir()
        (tmp_path / "goodpi" / "sequencing_data10").mkdir()

        result = fetchLatestSeqDir(_config(tmp_path), "goodpi")

        assert result == tmp_path / "goodpi" / "sequencing_data10"


class Test_matchOptdupsReqs:
    def _ssdf(self, **cols):
        return pd.DataFrame(cols)

    def test_normal_sample_computes_req_v_got_ratio(self):
        ssdf = self._ssdf(
            Sample_ID=["S1"], reqDepth=[1000000], gotDepth=[800000]
        )
        optDups = [["P1", "S1", "SampleOne", 1.5]]

        result = matchOptdupsReqs(optDups, ssdf)

        assert result == [["P1", "S1", "SampleOne", 1.5, 0.8, 800000]]

    def test_sample_missing_gotDepth_fills_zeroes(self):
        ssdf = self._ssdf(
            Sample_ID=["S1"], reqDepth=[1000000], gotDepth=[np.nan]
        )
        optDups = [["P1", "S1", "SampleOne", "NA"]]

        result = matchOptdupsReqs(optDups, ssdf)

        assert result == [["P1", "S1", "SampleOne", "NA", 0, 0]]

    def test_sorted_by_sampleID(self):
        ssdf = self._ssdf(
            Sample_ID=["S2", "S1"],
            reqDepth=[1000000, 1000000],
            gotDepth=[500000, 900000],
        )
        optDups = [
            ["P1", "S2", "SampleTwo", "NA"],
            ["P1", "S1", "SampleOne", "NA"],
        ]

        result = matchOptdupsReqs(optDups, ssdf)

        assert [r[1] for r in result] == ["S1", "S2"]

    def test_phix_multiple_depths_exits_when_inconsistent(self):
        # A sample (e.g. PhiX spiked into multiple lanes) can yield >1
        # gotDepth row; if those depths disagree, that's unexpected data
        # and matchOptdupsReqs must fail loudly rather than pick one.
        ssdf = self._ssdf(
            Sample_ID=["PhiX", "PhiX"],
            reqDepth=[1000000, 1000000],
            gotDepth=[500000, 600000],
        )
        optDups = [["P1", "PhiX", "PhiX", "NA"]]

        with pytest.raises(SystemExit):
            matchOptdupsReqs(optDups, ssdf)
