from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from wd40.release import release_rights


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
