import configparser
import subprocess as sp
from unittest.mock import Mock, patch

import pandas as pd
import os
import pytest
from dissectBCL.misc import joinLis
from dissectBCL.misc import hamming
from dissectBCL.misc import lenMask
from dissectBCL.misc import P5Seriesret
from dissectBCL.misc import retBCstr
from dissectBCL.misc import retIxtype
from dissectBCL.misc import retMean_perc_Q
from dissectBCL.misc import formatSeqRecipe
from dissectBCL.misc import formatMisMatches
from dissectBCL.misc import umlautDestroyer
from dissectBCL.misc import parseRunInfo
from dissectBCL.misc import getConf
from dissectBCL.misc import getNewFlowCell
from dissectBCL.misc import projectPI
from dissectBCL.misc import _resolve_internal_pis
from dissectBCL.misc import _fetch_ro_crate_metadata
from dissectBCL.misc import _build_ro_crate_archive
from dissectBCL.misc import _add_fastq_file_entities
from dissectBCL.misc import fexUpload
from zipfile import ZipFile, ZIP_STORED


def _write_test_ini(tmp_path, organizations="MPI-IE"):
    ini_path = tmp_path / "test.ini"
    ini_path.write_text(
        "[Internals]\n"
        f"Organizations={organizations}\n"
        "seqDir=seqfolderstr\n"
        "fex=False\n"
        "\n"
        "[parkour]\n"
        "user=parkourUser\n"
        "password=parkourPw\n"
        "cert=/path/to/cert.pem\n"
        "URL=https://parkour.domain.tld\n"
    )
    return ini_path


class Test_projectPI:
    def test_simple_surname(self):
        assert projectPI("Project_1234_jdoe_manke") == "manke"

    def test_compound_surname_is_kept_intact(self):
        # Parkour's internal_pis endpoint returns compound surnames in full
        # (verified against prod: 'cabezas-wallscheid'), so we must not truncate.
        assert (
            projectPI("Project_1234_jdoe_cabezas-wallscheid")
            == "cabezas-wallscheid"
        )

    def test_lowercases(self):
        assert projectPI("Project_1234_jdoe_Manke") == "manke"


class Test_getConf_internal_pis:
    @patch("dissectBCL.misc.requests.get")
    def test_resolves_pi_list_from_parkour(self, mock_get, tmp_path):
        mock_get.return_value = Mock(
            status_code=200,
            json=lambda: {"pis": ["Manke", "Cabezas"]},
        )
        ini_path = _write_test_ini(tmp_path)

        config = getConf(str(ini_path), quickload=True)

        mock_get.assert_called_once_with(
            "https://parkour.domain.tld/api/internal_pis/",
            params={"organizations": "MPI-IE"},
            auth=("parkourUser", "parkourPw"),
            verify="/path/to/cert.pem",
        )
        # Names are lowercased so they match the lowercased PI tokens the
        # shipping code compares against.
        assert config["Internals"]["PIs"] == "cabezas,manke"

    @patch("dissectBCL.misc.requests.get")
    def test_raises_loudly_on_parkour_failure(self, mock_get, tmp_path):
        # Parkour being unreachable must crash rather than degrade: an empty PI
        # list would misroute internal PIs to external FEX shipment. Restarting
        # once Parkour is back is the intended recovery.
        mock_get.side_effect = ConnectionError("network down")
        ini_path = _write_test_ini(tmp_path)

        with pytest.raises(RuntimeError):
            getConf(str(ini_path), quickload=True)

    @patch("dissectBCL.misc.requests.get")
    def test_raises_loudly_on_non_200_response(self, mock_get, tmp_path):
        mock_get.return_value = Mock(status_code=500, text="server error")
        ini_path = _write_test_ini(tmp_path)

        with pytest.raises(RuntimeError):
            getConf(str(ini_path), quickload=True)

    @patch("dissectBCL.misc.requests.get")
    def test_empty_pi_list_raises_runtimeerror(self, mock_get):
        # A 200 with an empty list (e.g. a misconfigured/bracketed Organizations
        # value) must crash rather than treat every PI as external.
        mock_get.return_value = Mock(status_code=200, json=lambda: {"pis": []})
        config = configparser.ConfigParser()
        config["Internals"] = {"Organizations": "[MPI-IE]"}
        config["parkour"] = {
            "URL": "https://parkour.domain.tld",
            "user": "u",
            "password": "p",
            "cert": "/cert.pem",
        }

        with pytest.raises(RuntimeError):
            _resolve_internal_pis(config)

    @patch("dissectBCL.misc.requests.get")
    def test_malformed_200_body_raises_runtimeerror(self, mock_get):
        # A 200 whose JSON lacks "pis" must surface as the descriptive
        # RuntimeError, not a bare KeyError.
        mock_get.return_value = Mock(
            status_code=200,
            json=lambda: {"unexpected": "shape"},
        )
        config = configparser.ConfigParser()
        config["Internals"] = {"Organizations": "MPI-IE"}
        config["parkour"] = {
            "URL": "https://parkour.domain.tld",
            "user": "u",
            "password": "p",
            "cert": "/cert.pem",
        }

        with pytest.raises(RuntimeError):
            _resolve_internal_pis(config)


class Test_getNewFlowCell_sequencer_gating:
    # A sequencer-restricted call must never touch the other platform's
    # Dirs keys - proven here by simply not defining them, so any read
    # would raise KeyError.

    def test_illumina_only_never_reads_aviti_keys(self, tmp_path):
        illumina_out = tmp_path / "out_illumina"
        illumina_base = tmp_path / "base_illumina"
        illumina_out.mkdir()
        illumina_base.mkdir()
        config = {
            "Dirs": {
                "outputDir_illumina": str(illumina_out),
                "baseDir_illumina": str(illumina_base),
            }
        }

        result = getNewFlowCell(config, None, "illumina")

        assert result == (None, None, None)

    def test_aviti_only_never_reads_illumina_keys(self, tmp_path):
        aviti_out = tmp_path / "out_aviti"
        aviti_base = tmp_path / "base_aviti"
        aviti_out.mkdir()
        aviti_base.mkdir()
        config = {
            "Dirs": {
                "outputDir_aviti": str(aviti_out),
                "baseDir_aviti": str(aviti_base),
            }
        }

        result = getNewFlowCell(config, None, "aviti")

        assert result == (None, None, None)


class Test_getNewFlowCell_aviti_serial_id_nesting:
    # baseDir_aviti holds one subdir per sequencer serial ID, each holding
    # its own flowcells - a flowcell is found two levels down, not one.
    def test_finds_flowcell_nested_under_serial_id_dir(self, tmp_path):
        out = tmp_path / "out"
        base = tmp_path / "base"
        out.mkdir()
        flowcellName = "20260804_AV261103_installpvrun-sideb-av261103"
        flowcellDir = base / "AV261103" / flowcellName
        flowcellDir.mkdir(parents=True)
        (flowcellDir / "RunUploaded.json").write_text('{"outcome": "OK"}')
        config = {
            "Dirs": {
                "outputDir_aviti": str(out),
                "baseDir_aviti": str(base),
            }
        }

        result = getNewFlowCell(config, None, "aviti")

        assert result == (flowcellName, flowcellDir, "aviti")

    def test_scans_multiple_serial_id_dirs(self, tmp_path):
        out = tmp_path / "out"
        base = tmp_path / "base"
        out.mkdir()
        # An empty serial-ID dir (e.g. a freshly added machine) must not
        # stop a real flowcell in a sibling serial-ID dir from being found.
        (base / "AV251009").mkdir(parents=True)
        flowcellName = "20260804_AV261103_installpvrun-sidea-av261103"
        flowcellDir = base / "AV261103" / flowcellName
        flowcellDir.mkdir(parents=True)
        (flowcellDir / "RunUploaded.json").write_text('{"outcome": "OK"}')
        config = {
            "Dirs": {
                "outputDir_aviti": str(out),
                "baseDir_aviti": str(base),
            }
        }

        result = getNewFlowCell(config, None, "aviti")

        assert result == (flowcellName, flowcellDir, "aviti")

    def test_completion_check_looks_under_matching_serial_id_output_dir(
        self, tmp_path
    ):
        # Output mirrors baseDir_aviti's serial-ID nesting: a flowcell
        # already marked done lives at outputDir_aviti/<serialID>/<name>*,
        # not flat under outputDir_aviti. The completion check must look
        # there, or it will re-offer an already-finished flowcell forever.
        out = tmp_path / "out"
        base = tmp_path / "base"
        flowcellName = "20260804_AV251009_run1"
        flowcellDir = base / "AV251009" / flowcellName
        flowcellDir.mkdir(parents=True)
        (flowcellDir / "RunUploaded.json").write_text('{"outcome": "OK"}')

        doneDir = out / "AV251009" / f"{flowcellName}_lanes_1"
        doneDir.mkdir(parents=True)
        (doneDir / "communication.done").touch()

        config = {
            "Dirs": {
                "outputDir_aviti": str(out),
                "baseDir_aviti": str(base),
            }
        }

        result = getNewFlowCell(config, None, "aviti")

        assert result == (None, None, None)


class Test_getConf_sequencer_gating:
    def _write_full_ini(self, tmp_path):
        adapters = tmp_path / "adapters.txt"
        adapters.write_text("")
        ini_path = tmp_path / "full.ini"
        ini_path.write_text(
            "[Internals]\n"
            "Organizations=MPI-IE\n"
            "seqDir=seqfolderstr\n"
            "fex=False\n"
            "\n"
            "[parkour]\n"
            "user=parkourUser\n"
            "password=parkourPw\n"
            "cert=/path/to/cert.pem\n"
            "URL=https://parkour.domain.tld\n"
            "\n"
            "[software]\n"
            "bclconvert=/usr/bin/bclconvert\n"
            "bases2fastq=/usr/bin/bases2fastq\n"
            f"fastqc_adapters={adapters}\n"
            "splitFastq=/usr/bin/splitFastq\n"
        )
        return ini_path

    @staticmethod
    def _run_side_effect(cmd, **kwargs):
        prog = cmd[0]
        if "bclconvert" in prog:
            return Mock(stdout=b"", stderr=b"bcl-convert Version 00.000.000.4.4.4\n")
        if "bases2fastq" in prog:
            return Mock(stdout="bases2fastq Version 2.1.0,\n", stderr="")
        if prog == "fastqc":
            return Mock(stdout=b"FastQC v0.12.1\n", stderr=b"")
        if prog == "kraken2":
            return Mock(stdout=b"Kraken version 2.1.3\n", stderr=b"")
        if prog == "clumpify.sh":
            return Mock(stdout=b"", stderr=b"BBTools version 39.01\n")
        raise AssertionError(f"unexpected command probed: {cmd}")

    @patch("dissectBCL.misc.requests.get")
    @patch("dissectBCL.misc.sp.run")
    @patch("dissectBCL.misc.version", return_value="1.0")
    def test_illumina_only_skips_bases2fastq_probe(
        self, mock_version, mock_run, mock_get, tmp_path
    ):
        mock_get.return_value = Mock(status_code=200, json=lambda: {"pis": ["Manke"]})
        mock_run.side_effect = self._run_side_effect
        ini_path = self._write_full_ini(tmp_path)

        config = getConf(str(ini_path), quickload=False, sequencer="illumina")

        probed = [call.args[0][0] for call in mock_run.call_args_list]
        assert "/usr/bin/bclconvert" in probed
        assert "/usr/bin/bases2fastq" not in probed
        assert "bclconvert" in config["softwareVers"]
        assert "bases2fastq" not in config["softwareVers"]
        assert "splitFastq" in config["softwareVers"]

    @patch("dissectBCL.misc.requests.get")
    @patch("dissectBCL.misc.sp.run")
    @patch("dissectBCL.misc.version", return_value="1.0")
    def test_aviti_only_skips_bclconvert_probe(
        self, mock_version, mock_run, mock_get, tmp_path
    ):
        mock_get.return_value = Mock(status_code=200, json=lambda: {"pis": ["Manke"]})
        mock_run.side_effect = self._run_side_effect
        ini_path = self._write_full_ini(tmp_path)

        config = getConf(str(ini_path), quickload=False, sequencer="aviti")

        probed = [call.args[0][0] for call in mock_run.call_args_list]
        assert "/usr/bin/bases2fastq" in probed
        assert "/usr/bin/bclconvert" not in probed
        assert "bases2fastq" in config["softwareVers"]
        assert "bclconvert" not in config["softwareVers"]
        assert "splitFastq" in config["softwareVers"]


class Test_ro_crate_archive:
    @patch("dissectBCL.misc.requests.get")
    def test_fetch_ro_crate_metadata_returns_graph_on_success(self, mock_get):
        mock_get.return_value = Mock(
            status_code=200,
            json=lambda: {"ro_crate": {"@graph": []}, "skipped_records": []},
        )
        config = configparser.ConfigParser()
        config["parkour"] = {
            "URL": "https://parkour.domain.tld",
            "user": "u",
            "password": "p",
            "cert": "/cert.pem",
        }

        result = _fetch_ro_crate_metadata("42", config)

        assert result == {"@graph": []}
        mock_get.assert_called_once_with(
            "https://parkour.domain.tld/api/generate_ro_crate/",
            params={"requests": "42", "preview": "true"},
            auth=("u", "p"),
            verify="/cert.pem",
        )

    @patch("dissectBCL.misc.requests.get")
    def test_fetch_ro_crate_metadata_returns_none_on_failure(self, mock_get):
        mock_get.side_effect = ConnectionError("down")
        config = configparser.ConfigParser()
        config["parkour"] = {
            "URL": "https://parkour.domain.tld",
            "user": "u",
            "password": "p",
            "cert": "/cert.pem",
        }

        result = _fetch_ro_crate_metadata("42", config)

        assert result is None

    def test_build_ro_crate_archive_contains_expected_files(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        fastqc_dir = tmp_path / "FASTQC_Project_42_jdoe_manke"
        project_dir.mkdir()
        fastqc_dir.mkdir()
        (project_dir / "sample_R1.fastq.gz").write_bytes(b"fake-gzip-bytes")
        (project_dir / "md5sums.txt").write_text("sample_R1.fastq.gz\tabc123\n")
        (fastqc_dir / "report.html").write_text("<html></html>")

        archive_path = _build_ro_crate_archive(
            "250101_M001_0001_AAAA",
            "Project_42_jdoe_manke",
            (project_dir, fastqc_dir),
            {"@graph": []},
        )

        try:
            with ZipFile(archive_path) as zf:
                names = set(zf.namelist())
                assert "Project_42_jdoe_manke/sample_R1.fastq.gz" in names
                assert "Project_42_jdoe_manke/md5sums.txt" in names
                assert "FASTQC_Project_42_jdoe_manke/report.html" in names
                assert "ro-crate-metadata.json" in names
                for info in zf.infolist():
                    assert info.compress_type == ZIP_STORED
        finally:
            archive_path.unlink()

    def test_build_ro_crate_archive_omits_metadata_file_when_none(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        fastqc_dir = tmp_path / "FASTQC_Project_42_jdoe_manke"
        project_dir.mkdir()
        fastqc_dir.mkdir()
        (project_dir / "sample_R1.fastq.gz").write_bytes(b"fake-gzip-bytes")

        archive_path = _build_ro_crate_archive(
            "250101_M001_0001_AAAA",
            "Project_42_jdoe_manke",
            (project_dir, fastqc_dir),
            None,
        )

        try:
            with ZipFile(archive_path) as zf:
                assert "ro-crate-metadata.json" not in zf.namelist()
        finally:
            archive_path.unlink()

    def test_build_ro_crate_archive_streams_to_fileobj_without_touching_disk(
        self, tmp_path
    ):
        import io

        project_dir = tmp_path / "Project_42_jdoe_manke"
        fastqc_dir = tmp_path / "FASTQC_Project_42_jdoe_manke"
        project_dir.mkdir()
        fastqc_dir.mkdir()
        (project_dir / "sample_R1.fastq.gz").write_bytes(b"fake-gzip-bytes")

        buffer = io.BytesIO()
        result = _build_ro_crate_archive(
            "250101_M001_0001_AAAA",
            "Project_42_jdoe_manke",
            (project_dir, fastqc_dir),
            {"@graph": []},
            fileobj=buffer,
        )

        assert result is None
        expected_archive = (
            tmp_path / "250101_M001_0001_AAAA_Project_42_jdoe_manke_ro_crate.zip"
        )
        assert not expected_archive.exists()

        buffer.seek(0)
        with ZipFile(buffer) as zf:
            names = set(zf.namelist())
            assert "Project_42_jdoe_manke/sample_R1.fastq.gz" in names
            assert "ro-crate-metadata.json" in names


class Test_fexUpload:
    @patch("dissectBCL.misc._build_ro_crate_archive")
    @patch("dissectBCL.misc._fetch_ro_crate_metadata")
    @patch("dissectBCL.misc.sp.Popen")
    @patch("dissectBCL.misc.sp.check_output")
    def test_streams_archive_into_fexsend_stdin_without_writing_to_disk(
        self, mock_check_output, mock_popen, mock_fetch, mock_build_archive, tmp_path
    ):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        fastqc_dir = tmp_path / "FASTQC_Project_42_jdoe_manke"
        project_dir.mkdir()
        fastqc_dir.mkdir()
        (project_dir / "sample_R1.fastq.gz").write_bytes(b"fake-gzip-bytes")

        mock_check_output.return_value = b""
        mock_fetch.return_value = None
        fake_proc = Mock()
        mock_popen.return_value = fake_proc

        config = configparser.ConfigParser()

        result = fexUpload(
            "250101_M001_0001_AAAA",
            "Project_42_jdoe_manke",
            "someone@example.com",
            (project_dir, fastqc_dir),
            config,
        )

        assert result == "Uploaded"
        mock_popen.assert_called_once_with(
            [
                "fexsend",
                "-s",
                "250101_M001_0001_AAAA_Project_42_jdoe_manke_ro_crate.zip",
                "someone@example.com",
            ],
            stdin=sp.PIPE,
        )
        mock_build_archive.assert_called_once_with(
            "250101_M001_0001_AAAA",
            "Project_42_jdoe_manke",
            (project_dir, fastqc_dir),
            None,
            fileobj=fake_proc.stdin,
        )
        fake_proc.stdin.close.assert_called_once()
        fake_proc.wait.assert_called_once()

        expected_archive = (
            tmp_path / "250101_M001_0001_AAAA_Project_42_jdoe_manke_ro_crate.zip"
        )
        assert not expected_archive.exists()

    @patch("dissectBCL.misc._build_ro_crate_archive")
    @patch("dissectBCL.misc._fetch_ro_crate_metadata")
    @patch("dissectBCL.misc.sp.Popen")
    @patch("dissectBCL.misc.sp.check_output")
    def test_kills_fexsend_and_reaps_it_if_archive_build_raises(
        self, mock_check_output, mock_popen, mock_fetch, mock_build_archive, tmp_path
    ):
        mock_check_output.return_value = b""
        mock_fetch.return_value = None
        mock_build_archive.side_effect = RuntimeError("boom")
        fake_proc = Mock()
        mock_popen.return_value = fake_proc

        config = configparser.ConfigParser()

        with pytest.raises(RuntimeError):
            fexUpload(
                "250101_M001_0001_AAAA",
                "Project_42_jdoe_manke",
                "someone@example.com",
                (tmp_path, tmp_path),
                config,
            )

        # fexsend is killed (rather than sent a clean EOF that would let it
        # finalize a truncated upload), stdin is closed, and the process is
        # reaped so it isn't left as a zombie.
        fake_proc.kill.assert_called_once()
        fake_proc.stdin.close.assert_called_once()
        fake_proc.wait.assert_called_once()


class Test_add_fastq_file_entities:
    def test_adds_file_entities_linked_to_matching_stub(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        sample_dir = project_dir / "Sample_24L000001"
        sample_dir.mkdir(parents=True)
        (sample_dir / "my-sample_R1.fastq.gz").write_bytes(b"fake-r1")
        (sample_dir / "my-sample_R2.fastq.gz").write_bytes(b"fake-r2")
        (project_dir / "md5sums.txt").write_text(
            "my-sample_R1.fastq.gz\tabc111\nmy-sample_R2.fastq.gz\tabc222\n"
        )
        ro_crate_metadata = {
            "@graph": [
                {"@id": "#fastq-data-24L000001", "@type": "Dataset", "hasPart": []},
            ]
        }

        _add_fastq_file_entities(ro_crate_metadata, project_dir)

        graph_ids = {entry["@id"] for entry in ro_crate_metadata["@graph"]}
        r1_id = "#fastq-file-24L000001-my-sample_R1.fastq.gz"
        r2_id = "#fastq-file-24L000001-my-sample_R2.fastq.gz"
        assert r1_id in graph_ids
        assert r2_id in graph_ids

        stub = next(
            e
            for e in ro_crate_metadata["@graph"]
            if e["@id"] == "#fastq-data-24L000001"
        )
        stub_part_ids = {ref["@id"] for ref in stub["hasPart"]}
        assert r1_id in stub_part_ids
        assert r2_id in stub_part_ids

        r1_entry = next(
            e for e in ro_crate_metadata["@graph"] if e["@id"] == r1_id
        )
        assert r1_entry["@type"] == ["File", "MediaObject"]
        assert r1_entry["name"] == "my-sample_R1.fastq.gz"
        assert (
            r1_entry["contentUrl"]
            == "Project_42_jdoe_manke/Sample_24L000001/my-sample_R1.fastq.gz"
        )
        assert r1_entry["encodingFormat"] == "application/gzip"
        md5_property_id = r1_entry["additionalProperty"][0]["@id"]
        md5_entry = next(
            e for e in ro_crate_metadata["@graph"] if e["@id"] == md5_property_id
        )
        assert md5_entry == {
            "@id": md5_property_id,
            "@type": "PropertyValue",
            "name": "md5",
            "value": "abc111",
        }

    def test_skips_files_without_a_matching_stub_and_does_not_raise(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        sample_dir = project_dir / "Sample_24L999999"
        sample_dir.mkdir(parents=True)
        (sample_dir / "orphan_R1.fastq.gz").write_bytes(b"fake-r1")
        ro_crate_metadata = {"@graph": []}

        _add_fastq_file_entities(ro_crate_metadata, project_dir)

        assert ro_crate_metadata["@graph"] == []

    def test_missing_md5sums_file_does_not_raise(self, tmp_path):
        project_dir = tmp_path / "Project_42_jdoe_manke"
        sample_dir = project_dir / "Sample_24L000001"
        sample_dir.mkdir(parents=True)
        (sample_dir / "my-sample_R1.fastq.gz").write_bytes(b"fake-r1")
        ro_crate_metadata = {
            "@graph": [
                {"@id": "#fastq-data-24L000001", "@type": "Dataset", "hasPart": []},
            ]
        }

        _add_fastq_file_entities(ro_crate_metadata, project_dir)

        r1_id = "#fastq-file-24L000001-my-sample_R1.fastq.gz"
        r1_entry = next(
            e for e in ro_crate_metadata["@graph"] if e["@id"] == r1_id
        )
        assert "additionalProperty" not in r1_entry


class Test_misc_data():
    def test_joinLis(self):
        assert joinLis([1, 2, 'A']) == "12A"
        assert joinLis([1, 2, 3]) == "123"
        assert joinLis(['a', 2, 'b'], joinStr=',') == 'a,2,b'

    def test_hamming(self):
        assert hamming(float(1), float(2)) == 0
        assert hamming('aabb', 'aaba') == 1
        assert hamming('aaaa', 'bbbb') == 4

    def test_lenMask(self):
        assert lenMask(8, 4, aviti=False) == "I4N4"
        assert lenMask(10, 10, aviti=False) == "I10"

    def test_P5Seriesret(self):
        adf = pd.DataFrame(
            {
                'index': [1, 2, 3, 4, 5],
                'index2': [1, 2, 3, 4, 5]
            }
        )

        bdf = pd.DataFrame(
            {
                'index': [1, 2, 3, 4, 5]
            }
        )
        for i, j in zip(P5Seriesret(adf), adf['index2']):
            assert i == j
        assert P5Seriesret(bdf).empty

    def test_retBCstr(self):
        _a = pd.Series(
            data=[1, 1],
            index=['index', 'index2']
        )
        _b = pd.Series(
            data=[1, 'bar'],
            index=['index', 'foo']
        )
        _c = pd.Series(
            data=1,
            index=['foo']
        )
        assert retBCstr(_a) == '1\t1'
        assert retBCstr(_b) == '1'
        assert retBCstr(_c) == 'nan'

    def test_retIxtype(self):
        _a = pd.Series(
            data=['I7type', 'I5type'],
            index=['I7_Index_ID', 'I5_Index_ID']
        )
        _b = pd.Series(
            data=['I7type'],
            index=['I7_Index_ID']
        )
        _c = pd.Series(
            data=['foo'],
            index=['bar']
        )
        assert retIxtype(_a) == 'I7type\tI5type'
        assert retIxtype(_b) == 'I7type'
        assert retIxtype(_c) == 'NA'

    def test_retMean_perc_Q(self):
        _a = pd.Series(
            data=['1:26.1,I1:22.8,2:23.0'],
            index=['meanQ']
        )
        _b = pd.Series(
            data=['1:48.2'],
            index=['meanQ']
        )
        assert retMean_perc_Q(_a) == '26.0\t23.0\t23.0'
        assert retMean_perc_Q(_b) == '48.0'
        assert retMean_perc_Q(_a, returnHeader=True) == (
            'R1_meanQ\tI1_meanQ\tR2_meanQ',
            '26.0\t23.0\t23.0'
        )

    def test_formatSeqRecipe(self):
        _a = {
            'Read1': ['Y', 100],
            'Read2': ['Y', 100]
        }
        _b = {
            'Read1': ['Y', 51],
            'Index1': ['I', 10]
        }
        assert formatSeqRecipe(_a) == "Read1:100; Read2:100"
        assert formatSeqRecipe(_b) == "Read1:51; Index1:10"

    def test_formatMisMatches(self):
        _a = {
            'BarcodeMismatchesIndex1': 2,
            'BarcodeMismatchesIndex2': 1
        }
        _b = {
            'BarcodeMismatchesIndex1': 1
        }
        assert formatMisMatches(
            _a
        ) == "BarcodeMismatchesIndex1:2, BarcodeMismatchesIndex2:1"
        assert formatMisMatches(_b) == "BarcodeMismatchesIndex1:1"

    def test_umlautDestroyer(self):
        _a = "ö"
        _b = "ä"
        _c = "ß"
        assert umlautDestroyer(_a) == "o"
        assert umlautDestroyer(_b) == "a"
        assert umlautDestroyer(_c) == "ss"


class Test_misc_Files():
    def RTF(self, testFile):
        return os.path.join(
            os.path.dirname(
                os.path.realpath(__file__)
            ),
            'test_misc',
            testFile
        )

    def test_parseRunInfo(self):
        _runInfo = parseRunInfo(
            self.RTF("RunInfo.xml")
        )
        _readDic = {
            'Read1': ['150', 'Read'],
            'Read2': ['8', 'Index'],
            'Read3': ['8', 'Index'],
            'Read4': ['150', 'Read']
        }
        assert _runInfo['instrument'] == 'NB000000'
        assert _runInfo['readDic'] == _readDic
        assert _runInfo['lanes'] == 4
        assert _runInfo['flowcellID'] == 'HHHHHHHHH'
