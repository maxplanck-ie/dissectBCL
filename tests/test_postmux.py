import hashlib
import os
from unittest.mock import patch

from dissectBCL.postmux import md5Runner
from dissectBCL.postmux import md5_multiqc


class TestMd5Runner:
    def test_matches_whole_file_digest(self, tmp_path):
        fqfile = tmp_path / "sample.fastq.gz"
        data = os.urandom(2 * 1024 * 1024 + 137)
        fqfile.write_bytes(data)

        expected = hashlib.md5(data).hexdigest()
        name, digest = md5Runner(fqfile)

        assert name == "sample.fastq.gz"
        assert digest == expected

    def test_handles_files_larger_than_chunk_size(self, tmp_path):
        fqfile = tmp_path / "big.fastq.gz"
        data = os.urandom(3 * 1024 * 1024)
        fqfile.write_bytes(data)

        _, digest = md5Runner(fqfile)

        assert digest == hashlib.md5(data).hexdigest()

    def test_empty_file(self, tmp_path):
        fqfile = tmp_path / "empty.fastq.gz"
        fqfile.write_bytes(b"")

        _, digest = md5Runner(fqfile)

        assert digest == hashlib.md5(b"").hexdigest()


class TestMd5Multiqc:
    @patch("dissectBCL.postmux.multiQC_yaml")
    @patch("dissectBCL.postmux.Popen")
    def test_writes_sorted_md5sums(self, mock_popen, mock_multiqc_yaml, tmp_path):
        mock_popen.return_value.wait.return_value = 0
        mock_multiqc_yaml.return_value = ({}, "", "", "", "")

        laneFolder = tmp_path
        project = "TestProject"
        projectFolder = laneFolder / f"Project_{project}"
        qcFolder = laneFolder / f"FASTQC_Project_{project}"
        sampleDir = projectFolder / "Sample_1"
        sampleDir.mkdir(parents=True)
        qcFolder.mkdir(parents=True)

        contents = {
            "b_R2.fastq.gz": b"second file contents",
            "a_R1.fastq.gz": b"first file contents",
        }
        for name, data in contents.items():
            (sampleDir / name).write_bytes(data)

        md5_multiqc(project, laneFolder, flowcell=object())

        md5out = projectFolder / "md5sums.txt"
        assert md5out.exists()

        lines = md5out.read_text().splitlines()
        names = [line.split("\t")[0] for line in lines]
        assert names == sorted(names)

        written = dict(line.split("\t") for line in lines)
        for name, data in contents.items():
            assert written[name] == hashlib.md5(data).hexdigest()

    @patch("dissectBCL.postmux.multiQC_yaml")
    @patch("dissectBCL.postmux.Popen")
    def test_skips_recompute_if_md5sums_exists(
        self, mock_popen, mock_multiqc_yaml, tmp_path
    ):
        mock_popen.return_value.wait.return_value = 0
        mock_multiqc_yaml.return_value = ({}, "", "", "", "")

        laneFolder = tmp_path
        project = "TestProject"
        projectFolder = laneFolder / f"Project_{project}"
        qcFolder = laneFolder / f"FASTQC_Project_{project}"
        sampleDir = projectFolder / "Sample_1"
        sampleDir.mkdir(parents=True)
        qcFolder.mkdir(parents=True)
        (sampleDir / "a.fastq.gz").write_bytes(b"contents")

        md5out = projectFolder / "md5sums.txt"
        md5out.write_text("preexisting\tdeadbeef\n")

        md5_multiqc(project, laneFolder, flowcell=object())

        assert md5out.read_text() == "preexisting\tdeadbeef\n"
