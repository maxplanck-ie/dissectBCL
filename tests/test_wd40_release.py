import configparser
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from wd40.release import (
    checkBRBDone,
    fetchLatestSeqDir,
    forceShip,
    release_rights,
)


def _make_tree(base):
    """
    multiqc_data/-shaped layout: most files "belong" to the project group,
    a couple land with the wrong group -- mirroring the real-world case
    where multiqc_sources.txt/multiqc_citations.txt/multiqc_data.json are
    written with the operator's primary group instead of the PI group.
    """
    d = base / "multiqc_data"
    d.mkdir()
    right = d / "deeptools_frag_size_table.txt"
    right.write_text("ok")
    wrong1 = d / "multiqc_sources.txt"
    wrong1.write_text("wrong")
    wrong2 = d / "multiqc_citations.txt"
    wrong2.write_text("wrong")
    return d, right, [wrong1, wrong2]


def _patch_group(
    monkeypatch, wrong_paths, wrong_group="bioinfo", right_group="akhtargrp"
):
    wrong_names = {p.name for p in wrong_paths}

    def fake_group(self):
        return wrong_group if self.name in wrong_names else right_group

    monkeypatch.setattr(Path, "group", fake_group)


def test_release_rights_chgrps_mismatched_files(tmp_path, monkeypatch):
    d, right, wrong = _make_tree(tmp_path)
    _patch_group(monkeypatch, wrong)

    chowned = []

    def fake_chown(path, uid, gid):
        chowned.append((path, uid, gid))

    with (
        patch(
            "wd40.release.grp_module.getgrnam",
            return_value=SimpleNamespace(gr_gid=4242),
        ),
        patch("wd40.release.os.chown", side_effect=fake_chown),
    ):
        successRate = release_rights(str(d), "akhtargrp")

    assert successRate == 1.0
    assert sorted(str(p) for p, _, _ in chowned) == sorted(str(p) for p in wrong)
    assert all(uid == -1 for _, uid, _ in chowned)
    assert all(gid == 4242 for _, _, gid in chowned)


def test_release_rights_reports_files_it_could_not_chgrp(tmp_path, monkeypatch, capsys):
    d, right, wrong = _make_tree(tmp_path)
    _patch_group(monkeypatch, wrong)

    with (
        patch(
            "wd40.release.grp_module.getgrnam",
            return_value=SimpleNamespace(gr_gid=4242),
        ),
        patch("wd40.release.os.chown", side_effect=PermissionError),
    ):
        successRate = release_rights(str(d), "akhtargrp")

    assert successRate == 1.0
    out = capsys.readouterr().out.replace("\n", "")
    assert "wrong grp" in out
    for p in wrong:
        assert str(p) in out


def test_release_rights_no_group_mismatch_skips_chown(tmp_path, monkeypatch):
    d, right, wrong = _make_tree(tmp_path)
    # Nothing is "wrong" this time -- everything reports the target group.
    monkeypatch.setattr(Path, "group", lambda self: "akhtargrp")

    with (
        patch(
            "wd40.release.grp_module.getgrnam",
            return_value=SimpleNamespace(gr_gid=4242),
        ),
        patch("wd40.release.os.chown") as mock_chown,
    ):
        successRate = release_rights(str(d), "akhtargrp")

    assert successRate == 1.0
    mock_chown.assert_not_called()


def test_release_rights_unknown_group_does_not_attempt_chown(
    tmp_path, monkeypatch, capsys
):
    d, right, wrong = _make_tree(tmp_path)
    _patch_group(monkeypatch, wrong)

    with (
        patch("wd40.release.grp_module.getgrnam", side_effect=KeyError),
        patch("wd40.release.os.chown") as mock_chown,
    ):
        successRate = release_rights(str(d), "akhtargrp")

    assert successRate == 1.0
    mock_chown.assert_not_called()
    out = capsys.readouterr().out.replace("\n", "")
    assert "wrong grp" in out


class Test_fetchLatestSeqDir:
    def test_single_match_returns_it_directly(self, tmp_path):
        (tmp_path / "goodpi" / "sequencing_data").mkdir(parents=True)

        result = fetchLatestSeqDir(str(tmp_path), "goodpi", "sequencing_data")

        assert result == str(tmp_path / "goodpi" / "sequencing_data")

    def test_multiple_matches_picks_highest_suffix(self, tmp_path):
        (tmp_path / "goodpi" / "sequencing_data2024").mkdir(parents=True)
        (tmp_path / "goodpi" / "sequencing_data2025").mkdir(parents=True)

        result = fetchLatestSeqDir(str(tmp_path), "goodpi", "sequencing_data")

        assert result == str(tmp_path / "goodpi" / "sequencing_data2025")

    def test_uses_configured_prefix(self, tmp_path):
        (tmp_path / "goodpi" / "seqfolderstr").mkdir(parents=True)
        (tmp_path / "goodpi" / "seqfolderstr2").mkdir()

        result = fetchLatestSeqDir(str(tmp_path), "goodpi", "seqfolderstr")

        assert result == str(tmp_path / "goodpi" / "seqfolderstr2")


class Test_checkBRBDone:
    def test_flag_present_prints_nothing(self, tmp_path, capsys):
        (tmp_path / "analysis.done").write_text("")

        checkBRBDone(tmp_path)

        assert capsys.readouterr().out == ""

    def test_flag_missing_warns(self, tmp_path, capsys):
        checkBRBDone(tmp_path)

        out = capsys.readouterr().out
        assert "analysis.done" in out
        assert str(tmp_path) in out


def test_force_ship_copies_and_releases_requested_project(tmp_path, monkeypatch):
    lane = tmp_path / "20260922_AV261103_2605514357_lanes_2"
    project = lane / "Project_4070_Hummel_Domschke"
    fastqc = lane / "FASTQC_Project_4070_Hummel_Domschke"
    sibling = lane / "Project_4071_other_user"
    project.mkdir(parents=True)
    fastqc.mkdir()
    sibling.mkdir()
    (project / "sample_R1.fastq.gz").write_bytes(b"project")
    (fastqc / "multiqc_report.html").write_text("report")

    config = configparser.ConfigParser()
    config["Dirs"] = {
        "piDir": str(tmp_path / "data"),
        "bioinfoCoreDir": str(tmp_path / "bioinfo"),
        "seqFacDir": str(tmp_path / "seqfac"),
    }
    config["Internals"] = {
        "PIs": "someone_else",
        "seqDir": "sequencing_data",
        "fex": "False",
    }
    config["communication"] = {"fromAddress": "from@example.com"}
    (tmp_path / "data" / "iovino" / "sequencing_data2").mkdir(parents=True)
    preexistingSibling = (
        tmp_path
        / "data"
        / "iovino"
        / "sequencing_data2"
        / lane.name
        / "Project_4071_existing"
    )
    preexistingSibling.mkdir(parents=True)
    preexistingSibling.chmod(0o700)
    monkeypatch.chdir(lane)

    with patch("dissectBCL.fakeNews.sendMqcReports") as mock_send_mqc:
        result = forceShip(".", "4070,iovino", config)

    mock_send_mqc.assert_called_once()
    assert mock_send_mqc.call_args.args[2] == project.name

    destination = (
        tmp_path
        / "data"
        / "iovino"
        / "sequencing_data2"
        / lane.name
    )
    copiedProject = destination / project.name
    copiedFastqc = destination / fastqc.name
    assert list(result["shipDic"]) == [project.name]
    assert result["shipDic"][project.name][0] == "Copied"
    assert (copiedProject / "sample_R1.fastq.gz").read_bytes() == b"project"
    assert (copiedFastqc / "multiqc_report.html").read_text() == "report"
    assert preexistingSibling.stat().st_mode & 0o777 == 0o700
    assert not (destination / sibling.name).exists()
    for path in (
        destination,
        copiedProject,
        copiedFastqc,
        copiedProject / "sample_R1.fastq.gz",
        copiedFastqc / "multiqc_report.html",
    ):
        assert path.stat().st_mode & 0o777 == 0o750
