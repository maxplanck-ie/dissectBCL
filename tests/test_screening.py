import configparser

import pytest

from dissectBCL.screening import needsEscalation, parseUnclassifiedPct, pickThreshold


def _write_report(path, unclassified_pct=None, rows=None):
    """
    Writes a minimal kraken2 --report TSV. kraken2's report columns are:
    pct, count_clade, count_direct, rank_code, taxid, name (indented
    with leading spaces per taxonomic depth).
    """
    lines = []
    if unclassified_pct is not None:
        lines.append(f"{unclassified_pct}\t100\t100\tU\t0\tunclassified")
    if rows:
        lines.extend(rows)
    else:
        lines.append(f"{100 - (unclassified_pct or 0)}\t900\t900\tR\t1\troot")
    path.write_text("\n".join(lines) + "\n")


class Test_parseUnclassifiedPct:
    def test_returns_the_unclassified_row_value(self, tmp_path):
        report = tmp_path / "sample.rep"
        _write_report(report, unclassified_pct=12.5)
        assert parseUnclassifiedPct(report) == 12.5

    def test_returns_none_when_no_unclassified_row(self, tmp_path):
        report = tmp_path / "sample.rep"
        report.write_text("100.0\t900\t900\tR\t1\troot\n")
        assert parseUnclassifiedPct(report) is None

    def test_returns_none_for_empty_report(self, tmp_path):
        report = tmp_path / "sample.rep"
        report.write_text("")
        assert parseUnclassifiedPct(report) is None


class Test_pickThreshold:
    @pytest.fixture
    def config(self):
        c = configparser.ConfigParser()
        c["screening"] = {
            "plusPFdb": "/dev/null",
            "unclassified_threshold": "10",
            "relaxed_library_types": "ATAC-Seq",
            "relaxed_threshold": "20",
        }
        return c

    def test_default_threshold_for_unlisted_library_type(self, config):
        assert pickThreshold("ChIP-Seq", config) == 10.0

    def test_relaxed_threshold_for_listed_library_type(self, config):
        assert pickThreshold("ATAC-Seq", config) == 20.0

    def test_relaxed_match_is_case_insensitive(self, config):
        assert pickThreshold("atac-seq", config) == 20.0

    def test_default_threshold_when_library_type_is_none(self, config):
        assert pickThreshold(None, config) == 10.0


class Test_needsEscalation:
    @pytest.fixture
    def config(self):
        c = configparser.ConfigParser()
        c["screening"] = {
            "plusPFdb": "/dev/null",
            "unclassified_threshold": "10",
            "relaxed_library_types": "ATAC-Seq",
            "relaxed_threshold": "20",
        }
        return c

    def test_escalates_over_default_threshold(self, tmp_path, config):
        report = tmp_path / "sample.rep"
        _write_report(report, unclassified_pct=15)
        assert needsEscalation(report, "ChIP-Seq", config) is True

    def test_does_not_escalate_under_default_threshold(self, tmp_path, config):
        report = tmp_path / "sample.rep"
        _write_report(report, unclassified_pct=5)
        assert needsEscalation(report, "ChIP-Seq", config) is False

    def test_atac_at_15pct_does_not_escalate_relaxed_threshold(self, tmp_path, config):
        report = tmp_path / "sample.rep"
        _write_report(report, unclassified_pct=15)
        assert needsEscalation(report, "ATAC-Seq", config) is False

    def test_atac_at_25pct_escalates_past_relaxed_threshold(self, tmp_path, config):
        report = tmp_path / "sample.rep"
        _write_report(report, unclassified_pct=25)
        assert needsEscalation(report, "ATAC-Seq", config) is True

    def test_empty_report_never_escalates(self, tmp_path, config):
        report = tmp_path / "sample.rep"
        report.write_text("")
        assert needsEscalation(report, "ChIP-Seq", config) is False
