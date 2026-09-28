import configparser
import logging
import subprocess as sp
from pathlib import Path
from unittest.mock import Mock, patch

import click
import pytest

from wd40.fex import (
    _add_fastq_file_entities,
    _fetch_ro_crate_metadata,
    _fex_archive_exists,
    _fex_delete,
    _project_size,
    _sample_barcodes,
    _warn_if_parkour_found_no_records,
    fex,
)

ARCHIVE = "Project_42_jdoe_manke_rocrate.zip"


def _config(url="https://parkour.example.org"):
    cfg = configparser.ConfigParser()
    cfg.read_dict(
        {
            "parkour": {
                "URL": url,
                "user": "u",
                "password": "p",
                "cert": "/nonexistent.pem",
            }
        }
    )
    return cfg


def _project(tmp_path):
    project_dir = tmp_path / "Project_42_jdoe_manke"
    sample_dir = project_dir / "Sample_24L000001"
    sample_dir.mkdir(parents=True)
    (sample_dir / "s_R1.fastq.gz").write_bytes(b"fake-r1")
    return project_dir


class Test_fetch_ro_crate_metadata:
    @patch("wd40.fex.requests.get")
    def test_queries_parkour_by_barcodes(self, mock_get):
        mock_get.return_value = Mock(
            raise_for_status=Mock(),
            json=Mock(return_value={"ro_crate": {"@graph": []}}),
        )

        result = _fetch_ro_crate_metadata(["26L000001", "26L000002"], _config())

        assert result == {"@graph": []}
        assert mock_get.call_args.kwargs["params"] == {
            "barcodes": "26L000001,26L000002",
            "preview": "true",
        }
        assert mock_get.call_args.kwargs["timeout"] == 60

    @patch("wd40.fex.requests.get")
    def test_uses_parkour_url_override(self, mock_get):
        mock_get.return_value = Mock(
            raise_for_status=Mock(),
            json=Mock(return_value={"ro_crate": {"@graph": []}}),
        )
        _fetch_ro_crate_metadata(
            ["X"],
            _config(url="https://default.example.org"),
            parkour_url="https://test.example.org",
        )
        url = mock_get.call_args.args[0]
        assert url.startswith("https://test.example.org/api/generate_ro_crate/")

    @patch("wd40.fex.requests.get")
    def test_returns_none_and_warns_on_failure(self, mock_get, caplog):
        import requests as real_requests

        mock_get.side_effect = real_requests.ConnectionError("nope")

        with caplog.at_level(logging.WARNING):
            result = _fetch_ro_crate_metadata(["26L000001"], _config())

        assert result is None
        assert "Shipping without ro-crate-metadata.json" in caplog.text

    @patch("wd40.fex.requests.get")
    def test_warns_when_parkour_returns_no_records_and_keeps_crate(
        self, mock_get, caplog
    ):
        hollow = {
            "@graph": [
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "description": "No matching barcodes or requests were found.",
                }
            ]
        }
        mock_get.return_value = Mock(
            raise_for_status=Mock(),
            json=Mock(return_value={"ro_crate": hollow}),
        )

        with caplog.at_level(logging.WARNING):
            result = _fetch_ro_crate_metadata(["26L004591"], _config())

        assert "Parkour returned no records for barcodes 26L004591" in caplog.text
        # still return the (hollow) crate so stubs can be synthesized from disk
        assert result is hollow


class Test_warn_if_parkour_found_no_records:
    def test_silent_for_a_normal_crate(self, caplog):
        crate = {"@graph": [{"@id": "./", "@type": "Dataset", "name": "r1"}]}
        with caplog.at_level(logging.WARNING):
            _warn_if_parkour_found_no_records(crate, ["26L000001"])
        assert caplog.text == ""


class Test_add_fastq_file_entities:
    def test_synthesizes_missing_stub_and_links_files(self, tmp_path, caplog):
        project_dir = _project(tmp_path)
        sample_dir = project_dir / "Sample_24L000001"
        (sample_dir / "my-sample_R2.fastq.gz").write_bytes(b"fake-r2")
        (project_dir / "md5sums.txt").write_text(
            "s_R1.fastq.gz\tabc111\nmy-sample_R2.fastq.gz\tabc222\n"
        )
        ro_crate_metadata = {
            "@graph": [
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "hasPart": [{"@id": "other-thing"}],
                }
            ]
        }

        with caplog.at_level(logging.WARNING):
            _add_fastq_file_entities(ro_crate_metadata, project_dir)

        # the missing stub is synthesized (old Parkour did not emit it) so
        # fastq files are still described in the delivered RO-Crate
        assert "Parkour did not provide 1 #fastq-data-* stub entities" in caplog.text
        stub = next(
            e
            for e in ro_crate_metadata["@graph"]
            if e["@id"] == "#fastq-data-24L000001"
        )
        assert stub["@type"] == "Dataset"
        assert stub["identifier"] == "urn:parkour:fastq-data:24L000001"
        stub_file_ids = {r["@id"] for r in stub["hasPart"]}
        assert "#fastq-file-24L000001-s_R1.fastq.gz" in stub_file_ids
        assert "#fastq-file-24L000001-my-sample_R2.fastq.gz" in stub_file_ids
        root = next(e for e in ro_crate_metadata["@graph"] if e["@id"] == "./")
        root_part_ids = {r["@id"] for r in root["hasPart"]}
        assert "#fastq-data-24L000001" in root_part_ids
        assert "other-thing" in root_part_ids

    def test_uses_parkour_stub_when_present(self, tmp_path, caplog):
        project_dir = _project(tmp_path)
        ro_crate_metadata = {
            "@graph": [
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "hasPart": [{"@id": "#fastq-data-24L000001"}],
                },
                {
                    "@id": "#fastq-data-24L000001",
                    "@type": "Dataset",
                    "hasPart": [],
                },
            ]
        }

        with caplog.at_level(logging.WARNING):
            _add_fastq_file_entities(ro_crate_metadata, project_dir)

        assert caplog.text == ""
        stub = next(
            e
            for e in ro_crate_metadata["@graph"]
            if e["@id"] == "#fastq-data-24L000001"
        )
        stub_file_ids = {r["@id"] for r in stub["hasPart"]}
        assert len(stub_file_ids) == 1
        # root hasPart is not duplicated by synthesis
        root = next(e for e in ro_crate_metadata["@graph"] if e["@id"] == "./")
        assert root["hasPart"].count({"@id": "#fastq-data-24L000001"}) == 1

    def test_missing_md5sums_file_does_not_raise(self, tmp_path):
        project_dir = _project(tmp_path)
        ro_crate_metadata = {"@graph": [{"@id": "./", "@type": "Dataset"}]}

        _add_fastq_file_entities(ro_crate_metadata, project_dir)

        r1_entry = next(
            e
            for e in ro_crate_metadata["@graph"]
            if e["@id"] == "#fastq-file-24L000001-s_R1.fastq.gz"
        )
        assert "additionalProperty" not in r1_entry


class Test_sample_barcodes:
    def test_returns_sorted_barcodes(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        (project_dir / "Sample_26L000003").mkdir(parents=True)
        (project_dir / "Sample_26L000001").mkdir()
        (project_dir / "Sample_26L000002").mkdir()
        (project_dir / "README.txt").write_text("not a sample")

        assert _sample_barcodes(project_dir) == [
            "26L000001",
            "26L000002",
            "26L000003",
        ]

    def test_empty_with_no_sample_folders(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        project_dir.mkdir()
        assert _sample_barcodes(project_dir) == []


class Test_project_size:
    def test_sums_all_files_recursively(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        (project_dir / "Sample_26L000001").mkdir(parents=True)
        (project_dir / "Sample_26L000001" / "a.fastq.gz").write_bytes(b"12345")
        (project_dir / "Sample_26L000001" / "b.fastq.gz").write_bytes(b"123")
        (project_dir / "md5sums.txt").write_bytes(b"md5")

        assert _project_size(project_dir) == 11

    def test_zero_for_empty_project(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        project_dir.mkdir()
        assert _project_size(project_dir) == 0


class Test_fex_delete:
    @patch("wd40.fex.sp.run")
    def test_deletes_archive_when_present(self, mock_run, tmp_path):
        with patch("wd40.fex._fex_archive_exists", return_value=True) as exists:
            _fex_delete(ARCHIVE, "someone@example.com")
        exists.assert_called_once_with(ARCHIVE, "someone@example.com")
        mock_run.assert_called_once_with(
            ["fexsend", "-d", ARCHIVE, "someone@example.com"], check=False
        )

    @patch("wd40.fex.sp.run")
    def test_does_not_delete_when_absent(self, mock_run, tmp_path):
        with patch("wd40.fex._fex_archive_exists", return_value=False):
            _fex_delete(ARCHIVE, "someone@example.com")
        mock_run.assert_not_called()


class Test_fex_archive_exists:
    @patch("wd40.fex.sp.check_output", return_value=b"#1) 100 MB [10 d] some_other.zip\n")
    def test_matches_full_archive_name(self, mock_check):
        assert _fex_archive_exists(ARCHIVE, "someone@example.com") is False
        assert _fex_archive_exists("some_other.zip", "someone@example.com") is True


class Test_fex:
    @patch("wd40.fex.shutil.which", return_value="/usr/bin/fexsend")
    @patch("wd40.fex.sp.Popen")
    @patch("wd40.fex._build_ro_crate_archive")
    @patch("wd40.fex._fetch_ro_crate_metadata")
    @patch("wd40.fex._fex_delete")
    def test_uploads_and_removes_stale_archive_first(
        self, mock_delete, mock_fetch, mock_build, mock_popen, mock_which, tmp_path
    ):
        project_dir = _project(tmp_path)
        mock_fetch.return_value = {"@graph": []}
        fake_proc = Mock()
        fake_proc.wait.return_value = 0
        mock_popen.return_value = fake_proc

        fex(str(project_dir), _config(), "someone@example.com")

        assert mock_popen.call_args.args[0] == [
            "/usr/bin/fexsend",
            "-s",
            ARCHIVE,
            "someone@example.com",
        ]
        assert mock_popen.call_args.kwargs == {"stdin": sp.PIPE}
        fake_proc.stdin.close.assert_called_once()
        fake_proc.wait.assert_called_once()
        # stale archive pre-deleted exactly once, no delete after success
        mock_delete.assert_called_once_with(ARCHIVE, "someone@example.com")

    @patch("wd40.fex.shutil.which", return_value="/usr/bin/fexsend")
    @patch("wd40.fex.sp.run")
    @patch("wd40.fex.sp.Popen")
    @patch("wd40.fex._build_ro_crate_archive")
    @patch("wd40.fex._fetch_ro_crate_metadata")
    @patch("wd40.fex._fex_delete")
    def test_uploads_large_project_as_regular_file_and_removes_temp_copy(
        self, mock_delete, mock_fetch, mock_build, mock_popen, mock_run, mock_which, tmp_path
    ):
        project_dir = _project(tmp_path)
        with open(project_dir / "Sample_24L000001" / "big.fastq.gz", "wb") as f:
            f.truncate(5 * 2**30)
        mock_fetch.return_value = {"@graph": []}
        mock_run.return_value = Mock(returncode=0)

        fex(str(project_dir), _config(), "someone@example.com")

        # streaming (-s) is not used for multi-GiB projects
        mock_popen.assert_not_called()
        args = mock_run.call_args.args[0]
        assert args[0] == "/usr/bin/fexsend"
        assert args[2] == "someone@example.com"
        tmp_zip = Path(args[1])
        assert tmp_zip.name == ARCHIVE
        assert tmp_zip.parent.name.startswith(".fextmp_")
        assert tmp_zip.parent.parent == tmp_path
        # the temp copy is removed; the archive lives only on the FEX server
        assert not tmp_zip.exists()
        assert not tmp_zip.parent.exists()
        # the zip was built into that temp file on disk
        built_into = mock_build.call_args.args[2]
        assert built_into.name == str(tmp_zip)

    @patch("wd40.fex.shutil.which", return_value="/usr/bin/fexsend")
    @patch("wd40.fex.sp.Popen")
    @patch("wd40.fex._build_ro_crate_archive")
    @patch("wd40.fex._fetch_ro_crate_metadata")
    @patch("wd40.fex._fex_delete")
    def test_retries_after_byte_mismatch_and_deletes_corrupt_upload(
        self, mock_delete, mock_fetch, mock_build, mock_popen, mock_which, tmp_path
    ):
        project_dir = _project(tmp_path)
        mock_fetch.return_value = {"@graph": []}
        proc1 = Mock()
        proc1.wait.return_value = 29
        proc2 = Mock()
        proc2.wait.return_value = 0
        mock_popen.side_effect = [proc1, proc2]

        fex(str(project_dir), _config(), "someone@example.com")

        assert mock_popen.call_count == 2
        assert mock_build.call_count == 2
        # pre-delete + deletion of the corrupt upload between attempts
        mock_delete.assert_called_with(ARCHIVE, "someone@example.com")
        assert mock_delete.call_count == 2

    @patch("wd40.fex.shutil.which", return_value="/usr/bin/fexsend")
    @patch("wd40.fex.sp.Popen")
    @patch("wd40.fex._build_ro_crate_archive")
    @patch("wd40.fex._fetch_ro_crate_metadata")
    @patch("wd40.fex._fex_delete")
    def test_aborts_if_both_attempts_fail(
        self, mock_delete, mock_fetch, mock_build, mock_popen, mock_which, tmp_path
    ):
        project_dir = _project(tmp_path)
        mock_fetch.return_value = {"@graph": []}
        mock_popen.side_effect = [
            Mock(wait=Mock(return_value=29)),
            Mock(wait=Mock(return_value=29)),
        ]

        with pytest.raises(click.Abort):
            fex(str(project_dir), _config(), "someone@example.com")

        assert mock_popen.call_count == 2
        # pre-delete + one delete after each failed attempt
        assert mock_delete.call_count == 3

    @patch("wd40.fex._fetch_ro_crate_metadata")
    @patch("wd40.fex._fex_delete")
    @patch("wd40.fex._build_ro_crate_archive")
    @patch("wd40.fex.sp.Popen")
    def test_kills_fexsend_and_cleans_up_if_archive_build_raises(
        self, mock_popen, mock_build, mock_delete, mock_fetch, tmp_path
    ):
        project_dir = _project(tmp_path)
        mock_fetch.return_value = {"@graph": []}
        mock_build.side_effect = RuntimeError("boom")
        fake_proc = Mock()
        mock_popen.return_value = fake_proc

        with pytest.raises(RuntimeError):
            fex(str(project_dir), _config(), "someone@example.com")

        fake_proc.kill.assert_called_once()
        # stale archive + partial upload cleaned up
        mock_delete.assert_called_with(ARCHIVE, "someone@example.com")
        assert mock_delete.call_count == 2

    @patch("wd40.fex.sp.Popen")
    def test_returns_without_uploading_for_unparseable_project_name(
        self, mock_popen, tmp_path
    ):
        project_dir = tmp_path / "NotAProject_Dir"
        project_dir.mkdir()

        fex(str(project_dir), _config(), "someone@example.com")

        mock_popen.assert_not_called()