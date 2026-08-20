import configparser
from pathlib import Path
from unittest.mock import patch

from dissectBCL.fakeNews import shipFiles


def _write_test_config(bioinfo_dir, seqfac_dir):
    config = configparser.ConfigParser()
    config["Internals"] = {
        "PIs": "goodpi,brokenpi",
        "seqDir": "sequencing_data",
        "fex": "False",
    }
    config["communication"] = {"fromAddress": "someone@example.com"}
    config["Dirs"] = {
        "bioinfoCoreDir": str(bioinfo_dir),
        "seqFacDir": str(seqfac_dir),
    }
    return config


def _make_project(outPath, project, pi_seqdata_base):
    """
    Create a Project_*/FASTQC_Project_* pair under outPath, mirroring what
    postmux() leaves behind for shipFiles() to pick up.
    """
    projectPath = outPath / project
    fqcPath = outPath / project.replace("Project_", "FASTQC_Project_")
    projectPath.mkdir()
    fqcPath.mkdir()
    (projectPath / "sample_R1.fastq.gz").write_bytes(b"fake-gzip-bytes")
    (fqcPath / "multiqc_report.html").write_text("<html></html>")
    pi_seqdata_base.mkdir(parents=True, exist_ok=True)
    return projectPath, fqcPath


def _fake_copytree(src, dst):
    """
    Stand-in for shutil.copytree that fails for the "broken" project's
    trees (simulating e.g. a disk error), and does a trivial real copy
    otherwise, without recursing back into the patched shutil.copytree.
    """
    src = Path(src)
    if "brokenpi" in src.name:
        raise OSError("disk full (simulated)")
    dst = Path(dst)
    dst.mkdir(parents=True)
    for f in src.iterdir():
        (dst / f.name).write_bytes(f.read_bytes())


class Test_shipFiles_per_project_isolation:
    @patch("dissectBCL.fakeNews.mailHome")
    @patch("dissectBCL.fakeNews.fetchLatestSeqDir")
    @patch("dissectBCL.fakeNews.shutil.copytree", side_effect=_fake_copytree)
    def test_one_project_failing_does_not_block_siblings(
        self, mock_copytree, mock_fetchLatestSeqDir, mock_mailHome, tmp_path
    ):
        outLane = "250101_M001_0001_AAAA_lanes_1"
        outPath = tmp_path / outLane
        outPath.mkdir()

        good_base = tmp_path / "data" / "goodpi" / "sequencing_data"
        broken_base = tmp_path / "data" / "brokenpi" / "sequencing_data"
        mock_fetchLatestSeqDir.side_effect = lambda config, PI: (
            good_base if PI == "goodpi" else broken_base
        )

        _make_project(outPath, "Project_1_jdoe_goodpi", good_base)
        _make_project(outPath, "Project_2_jdoe_brokenpi", broken_base)

        bioinfo_dir = tmp_path / "bioinfo"
        bioinfo_dir.mkdir()
        config = _write_test_config(bioinfo_dir, tmp_path / "seqfac")

        result = shipFiles(outPath, config)

        # The broken project is reported as failed, and doesn't take the
        # whole run down with it.
        assert result["failedProjects"] == ["Project_2_jdoe_brokenpi"]
        failedEntry = result["shipDic"]["Project_2_jdoe_brokenpi"]
        assert failedEntry["status"] == "FAILED"
        assert "disk full" in failedEntry["error"]

        # The sibling project still shipped successfully.
        goodEntry = result["shipDic"]["Project_1_jdoe_goodpi"]
        assert goodEntry[0] == "Copied"
        assert (good_base / outLane / "Project_1_jdoe_goodpi").exists()
        assert (good_base / outLane / "FASTQC_Project_1_jdoe_goodpi").exists()

        # The broken project never got anything copied into place.
        assert not (broken_base / outLane / "Project_2_jdoe_brokenpi").exists()
        assert not (broken_base / outLane / "FASTQC_Project_2_jdoe_brokenpi").exists()

        # Failure is loud: a dedicated email was sent for the failed project.
        mock_mailHome.assert_called_once()
        subject = mock_mailHome.call_args.args[0]
        assert "SHIPPING FAILED" in subject
        assert "Project_2_jdoe_brokenpi" in subject


class Test_shipFiles_deliverTo_membership:
    @patch("dissectBCL.fakeNews.fetchLatestSeqDir")
    def test_deliver_to_key_ships_internally_even_when_absent_from_pi_list(
        self, mock_fetchLatestSeqDir, tmp_path
    ):
        # Project string still carries the PI's previous Parkour name
        # ("cabezas-wallscheid"); only a deliverTo override maps it to the
        # current internal PI dir - the raw token is not in [Internals] PIs.
        outLane = "260818_A00931_0932_BHNHMNDRX7_lanes_2"
        outPath = tmp_path / outLane
        outPath.mkdir()

        seq_base = tmp_path / "data" / "cabezas" / "sequencing_data"
        mock_fetchLatestSeqDir.return_value = seq_base

        _make_project(outPath, "Project_4035_Demollin_Cabezas-Wallscheid", seq_base)

        bioinfo_dir = tmp_path / "bioinfo"
        bioinfo_dir.mkdir()
        config = configparser.ConfigParser()
        config["Internals"] = {
            "PIs": "cabezas,akhtar",
            "seqDir": "sequencing_data",
            "fex": "False",
            "deliverTo": '{"cabezas-wallscheid": "cabezas"}',
        }
        config["communication"] = {"fromAddress": "someone@example.com"}
        config["Dirs"] = {
            "bioinfoCoreDir": str(bioinfo_dir),
            "seqFacDir": str(tmp_path / "seqfac"),
        }

        result = shipFiles(outPath, config)

        assert result["failedProjects"] == []
        entry = result["shipDic"]["Project_4035_Demollin_Cabezas-Wallscheid"]
        assert entry[0] == "Copied"
        assert (
            seq_base / outLane / "Project_4035_Demollin_Cabezas-Wallscheid"
        ).exists()
