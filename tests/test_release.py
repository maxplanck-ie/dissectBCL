from unittest.mock import Mock, patch

from wd40.release import fetchFolders, rel


@patch("wd40.release.requests.post")
@patch("wd40.release.release_folder")
@patch("wd40.release.checkBRBDone")
@patch("wd40.release.fetchFolders")
class Test_rel_deliverTo:
    """
    rel()'s put_filepaths step re-derives the PI token from the released
    project path and must honour the same deliverTo override fetchFolders
    already applies (parkour2 #317/#294), instead of a hardcoded per-PI
    special case (previously only "cabezas-wallscheid" -> "cabezas").
    """

    def _projDic(self, project_pi_segment):
        proj = f"Project_1234_jdoe_{project_pi_segment}"
        return {
            proj: [
                f"{project_pi_segment}grp",
                [
                    "/flowcellF",
                    f"/data/pidir/{proj}",
                    "/fqcF",
                    "/analysisF",
                ],
            ]
        }

    def _call_rel(self, piList, deliverTo, project_pi_segment):
        rel(
            "/flowcell",
            piList,
            "/prefix",
            "postfix",
            "https://parkour.tld",
            ("u", "p"),
            "/cert.pem",
            True,
            "from@test.io",
            deliverTo=deliverTo,
        )

    def test_deliver_to_override_resolves_full_pi_name_to_its_short_dir(
        self, mock_fetchFolders, mock_checkBRBDone, mock_release_folder, mock_post
    ):
        # On-disk project name still carries the PI's full name; only the
        # sequencing_data directory uses the short deliver_to token. The
        # membership check must still recognize it once resolved.
        mock_fetchFolders.return_value = self._projDic("Cabezas-Wallscheid")
        mock_release_folder.return_value = [1.0, 1.0]
        mock_post.return_value = Mock(status_code=200)

        self._call_rel(
            piList="cabezas-wallscheid,manke",
            deliverTo={"cabezas-wallscheid": "cabezas"},
            project_pi_segment="Cabezas-Wallscheid",
        )

        assert mock_post.called

    def test_pi_with_no_override_falls_back_to_its_own_name(
        self, mock_fetchFolders, mock_checkBRBDone, mock_release_folder, mock_post
    ):
        mock_fetchFolders.return_value = self._projDic("manke")
        mock_release_folder.return_value = [1.0, 1.0]
        mock_post.return_value = Mock(status_code=200)

        self._call_rel(piList="manke", deliverTo={}, project_pi_segment="manke")

        assert mock_post.called

    def test_no_deliver_to_dict_defaults_to_empty(
        self, mock_fetchFolders, mock_checkBRBDone, mock_release_folder, mock_post
    ):
        # deliverTo defaults to None on the CLI path when config has no
        # overrides at all - must not raise (e.g. AttributeError on None).
        mock_fetchFolders.return_value = self._projDic("manke")
        mock_release_folder.return_value = [1.0, 1.0]
        mock_post.return_value = Mock(status_code=200)

        self._call_rel(piList="manke", deliverTo=None, project_pi_segment="manke")

        assert mock_post.called

    def test_pi_in_piList_is_exact_not_substring(
        self, mock_fetchFolders, mock_checkBRBDone, mock_release_folder, mock_post
    ):
        # "be" must not match because it's a substring of "alhajabed" in the
        # piList string - membership must be exact, per-name.
        mock_fetchFolders.return_value = self._projDic("be")
        mock_release_folder.return_value = [1.0, 1.0]
        mock_post.return_value = Mock(status_code=200)

        self._call_rel(piList="alhajabed,manke", deliverTo={}, project_pi_segment="be")

        assert not mock_post.called


class Test_fetchFolders_deliverTo:
    """
    fetchFolders()'s internal/external membership check must recognize a PI
    via a deliver_to key too, not just Parkour's current PI list - this
    covers a PI whose Parkour name changed after some of their projects were
    already created (e.g. "Cabezas-Wallscheid" -> "Cabezas"), so those older
    projects still carry the PI's previous name and would otherwise be
    treated as external and fex'ed.
    """

    def test_membership_recognizes_a_deliver_to_key_not_in_institute_pis(
        self, tmp_path
    ):
        flowcellPath = tmp_path / "260818_A00931_0932_BHNHMNDRX7_lanes_2"
        proj = "Project_4035_Demollin_Cabezas-Wallscheid"
        (flowcellPath / proj).mkdir(parents=True)

        seqDir = tmp_path / "data" / "cabezas" / "sequencing_data"
        (seqDir / flowcellPath.name).mkdir(parents=True)

        projDic = fetchFolders(
            str(flowcellPath),
            "manke",  # institute_PIs has neither cabezas-wallscheid nor cabezas
            str(tmp_path / "data"),
            "sequencing_data",
            False,
            ("https://parkour.tld", ("u", "p"), "/cert.pem", "from@test.io"),
            deliverTo={"cabezas-wallscheid": "cabezas"},
        )

        assert proj in projDic
        grp, paths = projDic[proj]
        # The group and every path must use the resolved "cabezas" token,
        # not the raw "cabezas-wallscheid" one.
        assert grp == "cabezasgrp"
        assert paths[0] == str(seqDir / flowcellPath.name)

    def test_no_override_and_not_in_institute_pis_falls_to_fex(self, tmp_path):
        flowcellPath = tmp_path / "260101_M001_0001_AAAA"
        proj = "Project_1_jdoe_unknownpi"
        (flowcellPath / proj).mkdir(parents=True)

        with patch("wd40.release.check_output") as mock_check_output:
            mock_check_output.return_value = b""
            projDic = fetchFolders(
                str(flowcellPath),
                "manke",
                str(tmp_path / "data"),
                "sequencing_data",
                False,
                ("https://parkour.tld", ("u", "p"), "/cert.pem", "from@test.io"),
                deliverTo={},
            )

        assert proj not in projDic
