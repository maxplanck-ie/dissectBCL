import configparser
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
from dissectBCL.misc import _fetch_ro_crate_metadata
from dissectBCL.misc import _build_ro_crate_archive
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


class Test_getConf_internal_pis:
    @patch("dissectBCL.misc.requests.get")
    def test_resolves_pi_list_from_parkour(self, mock_get, tmp_path):
        mock_get.return_value = Mock(
            status_code=200,
            json=lambda: {"pis": ["manke", "cabezas"]},
        )
        ini_path = _write_test_ini(tmp_path)

        config = getConf(str(ini_path), quickload=True)

        mock_get.assert_called_once_with(
            "https://parkour.domain.tld/api/internal_pis/",
            params={"organizations": "MPI-IE"},
            auth=("parkourUser", "parkourPw"),
            verify="/path/to/cert.pem",
        )
        assert config["Internals"]["PIs"] == "cabezas,manke"

    @patch("dissectBCL.misc.requests.get")
    def test_raises_loudly_on_parkour_failure(self, mock_get, tmp_path):
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
