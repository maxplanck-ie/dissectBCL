from unittest.mock import Mock, patch

from wd40.release import rel


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
