import configparser
import json
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from dissectBCL.fakeNews import buildContaminationDic, pushParkour, shipFiles


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


class _FakeSampleSheet:
    def __init__(self, ssDic):
        self.ssDic = ssDic


def _parkour_config(url="https://parkour.example.com"):
    config = configparser.ConfigParser()
    config["parkour"] = {
        "URL": url,
        "user": "u",
        "password": "p",
        "cert": "",
    }
    return config


class Test_pushParkour_aviti_outBaseDir:
    @patch("dissectBCL.fakeNews.requests.post")
    def test_reads_RunStats_from_outBaseDir_not_config_outputDir(
        self, mock_post, tmp_path
    ):
        """
        Aviti output nests under a serial-ID subdir (e.g. AV251009), which
        outputDir_aviti alone does not include. pushParkour must read
        RunStats.json from the caller-supplied outBaseDir (the actual,
        nested location), not reconstruct a flat path from config.
        """
        mock_post.return_value.status_code = 200
        outLane = "20260804_AV251009_run1_lanes_1"
        flatOutputDir = tmp_path / "flat_outputDir_aviti"
        flatOutputDir.mkdir()
        nestedOutBaseDir = tmp_path / "nested" / "AV251009"
        nestedOutBaseDir.mkdir(parents=True)
        (nestedOutBaseDir / outLane).mkdir()
        (nestedOutBaseDir / outLane / "RunStats.json").write_text(
            json.dumps(
                {
                    "Lanes": [
                        {
                            "Lane": 1,
                            "NumPolonies": 1000,
                            "Reads": [{"PercentQ30": 90.0}],
                            "PercentQ30": 95.0,
                            "PercentAssignedReads": 80.0,
                        }
                    ]
                }
            )
        )

        config = configparser.ConfigParser()
        config["Dirs"] = {"outputDir_aviti": str(flatOutputDir)}
        config["parkour"] = {
            "URL": "https://parkour.example.com",
            "user": "u",
            "password": "p",
            "cert": "",
        }
        sampleSheet = _FakeSampleSheet({outLane: {}})

        pushParkour(
            "20260804_AV251009_run1",
            sampleSheet,
            config,
            None,
            "aviti",
            outBaseDir=nestedOutBaseDir,
        )

        mock_post.assert_called_once()

    @patch("dissectBCL.fakeNews.requests.post")
    def test_payload_flowcell_id_and_matrix_content(self, mock_post, tmp_path):
        """
        The POST body must carry the *last* underscore-delimited token of
        the flowcell folder name as flowcell_id (Aviti has no leading-letter
        or hyphen stripping, unlike Illumina), and a matrix whose fields are
        derived correctly from RunStats.json: reads_pf from NumPolonies,
        read_1/read_2 from the two Reads entries, undetermined_indices as
        the complement of PercentAssignedReads.
        """
        mock_post.return_value.status_code = 200
        outLane = "run1_lanes_1"
        outBaseDir = tmp_path / "aviti_out"
        outBaseDir.mkdir()
        (outBaseDir / outLane).mkdir()
        (outBaseDir / outLane / "RunStats.json").write_text(
            json.dumps(
                {
                    "Lanes": [
                        {
                            "Lane": 1,
                            "NumPolonies": 1234,
                            "Reads": [
                                {"PercentQ30": 91.5},
                                {"PercentQ30": 88.25},
                            ],
                            "PercentQ30": 90.0,
                            "PercentAssignedReads": 97.3,
                        }
                    ]
                }
            )
        )
        config = _parkour_config()
        sampleSheet = _FakeSampleSheet({outLane: {}})

        pushParkour(
            "20260804_AV251009_myrun",
            sampleSheet,
            config,
            None,
            "aviti",
            outBaseDir=outBaseDir,
        )

        _, kwargs = mock_post.call_args
        assert kwargs["data"]["flowcell_id"] == "myrun"
        matrix = json.loads(kwargs["data"]["matrix"])
        assert matrix == [
            {
                "reads_pf": 1234,
                "read_1": 91.5,
                "read_2": 88.25,
                "cluster_pf": 90.0,
                "undetermined_indices": 2.7,
                "name": "Lane 1",
            }
        ]

    @patch("dissectBCL.fakeNews.requests.post")
    def test_single_read_leaves_read_2_none(self, mock_post, tmp_path):
        """Single-end Aviti runs only report one Reads entry."""
        mock_post.return_value.status_code = 200
        outLane = "run1_lanes_1"
        outBaseDir = tmp_path / "aviti_out"
        outBaseDir.mkdir()
        (outBaseDir / outLane).mkdir()
        (outBaseDir / outLane / "RunStats.json").write_text(
            json.dumps(
                {
                    "Lanes": [
                        {
                            "Lane": 1,
                            "NumPolonies": 500,
                            "Reads": [{"PercentQ30": 80.0}],
                            "PercentQ30": 80.0,
                            "PercentAssignedReads": 90.0,
                        }
                    ]
                }
            )
        )
        config = _parkour_config()
        sampleSheet = _FakeSampleSheet({outLane: {}})

        pushParkour(
            "20260804_AV251009_myrun",
            sampleSheet,
            config,
            None,
            "aviti",
            outBaseDir=outBaseDir,
        )

        matrix = json.loads(mock_post.call_args.kwargs["data"]["matrix"])
        assert matrix[0]["read_2"] is None

    @patch("dissectBCL.fakeNews.requests.post")
    def test_return_value_is_passed_through(self, mock_post, tmp_path):
        """Callers (flowCellClass.fakenews) rely on the response being returned."""
        outLane = "run1_lanes_1"
        outBaseDir = tmp_path / "aviti_out"
        outBaseDir.mkdir()
        (outBaseDir / outLane).mkdir()
        (outBaseDir / outLane / "RunStats.json").write_text(
            json.dumps(
                {
                    "Lanes": [
                        {
                            "Lane": 1,
                            "NumPolonies": 500,
                            "Reads": [{"PercentQ30": 80.0}],
                            "PercentQ30": 80.0,
                            "PercentAssignedReads": 90.0,
                        }
                    ]
                }
            )
        )
        sentinel = Mock(status_code=200)
        mock_post.return_value = sentinel
        config = _parkour_config()
        sampleSheet = _FakeSampleSheet({outLane: {}})

        result = pushParkour(
            "20260804_AV251009_myrun",
            sampleSheet,
            config,
            None,
            "aviti",
            outBaseDir=outBaseDir,
        )

        assert result is sentinel

    @patch("dissectBCL.fakeNews.mailHome")
    @patch("dissectBCL.fakeNews.requests.post")
    def test_non_200_response_mails_but_does_not_abort(
        self, mock_post, mock_mailHome, tmp_path
    ):
        """A non-200 Parkour response must not be swallowed silently, but
        pushing stats is best-effort: unlike pullParkour it must not raise
        or otherwise take the pipeline down."""
        mock_post.return_value.status_code = 500
        outLane = "run1_lanes_1"
        outBaseDir = tmp_path / "aviti_out"
        outBaseDir.mkdir()
        (outBaseDir / outLane).mkdir()
        (outBaseDir / outLane / "RunStats.json").write_text(
            json.dumps(
                {
                    "Lanes": [
                        {
                            "Lane": 1,
                            "NumPolonies": 500,
                            "Reads": [{"PercentQ30": 80.0}],
                            "PercentQ30": 80.0,
                            "PercentAssignedReads": 90.0,
                        }
                    ]
                }
            )
        )
        config = _parkour_config()
        sampleSheet = _FakeSampleSheet({outLane: {}})

        result = pushParkour(
            "20260804_AV251009_myrun",
            sampleSheet,
            config,
            None,
            "aviti",
            outBaseDir=outBaseDir,
        )

        assert result.status_code == 500
        mock_mailHome.assert_called_once()
        assert "500" in mock_mailHome.call_args.args[1]


class Test_pushParkour_illumina:
    def _write_quality_metrics(self, path):
        pd.DataFrame(
            {
                "Lane": [1, 1, 1, 1],
                "SampleID": ["Undetermined", "Undetermined", "Sample1", "Sample1"],
                "YieldQ30": [100, 100, 900, 900],
                "ReadNumber": [1, 2, 1, 2],
                "% Q30": [0.5, 0.5, 0.9, 0.85],
                "Yield": [200, 200, 1000, 1000],
            }
        ).to_csv(path, index=False)

    @patch("dissectBCL.fakeNews.requests.post")
    @patch("dissectBCL.fakeNews.interop.summary")
    @patch("dissectBCL.fakeNews.interop.read")
    def test_computed_lane_metrics_and_flowcell_id_hyphen_stripping(
        self, mock_iop_read, mock_iop_summary, mock_post, tmp_path
    ):
        mock_post.return_value.status_code = 200
        outLane = "run1_lanes_1"
        outputDir = tmp_path / "illumina_out"
        reportsDir = outputDir / outLane / "Reports"
        reportsDir.mkdir(parents=True)
        self._write_quality_metrics(reportsDir / "Quality_Metrics.csv")

        mock_iop_summary.return_value = [
            {"ReadNumber": 1, "Lane": 1, "Reads Pf": 5000000.0},
        ]

        config = _parkour_config()
        config["Dirs"] = {"outputDir_illumina": str(outputDir)}
        sampleSheet = _FakeSampleSheet({outLane: {}})

        pushParkour(
            "some-HCCMWDRXY",
            sampleSheet,
            config,
            tmp_path / "flowcellBase",
            "illumina",
        )

        _, kwargs = mock_post.call_args
        # FID is split on "-" and the second half is taken.
        assert kwargs["data"]["flowcell_id"] == "HCCMWDRXY"
        matrix = json.loads(kwargs["data"]["matrix"])
        assert matrix == [
            {
                "reads_pf": 5000000.0,
                "undetermined_indices": 10.0,
                "read_1": 90.0,
                "read_2": 85.0,
                "cluster_pf": 90.0,
                "name": "Lane 1",
            }
        ]


class Test_buildContaminationDic:
    def _ssdf(self, sampleID, organism="mouse (GRCm39)"):
        return pd.DataFrame({"Sample_ID": [sampleID], "Organism": [[organism]]})

    def test_reads_primary_report_only_when_no_escalation(self, tmp_path):
        outPath = tmp_path / "lane"
        sampleDir = outPath / "FASTQC_Project_1_proj" / "Sample_S1"
        sampleDir.mkdir(parents=True)
        (sampleDir / "S1.rep").write_text(
            "5.0\t50\t50\tU\t0\tunclassified\n95.0\t950\t950\tS\t10090\tmouse\n"
        )

        result = buildContaminationDic(outPath, self._ssdf("S1"))

        assert result["S1"][0] == 0.95  # fraction: top hit (950) / total (1000)
        assert result["S1"][1] == "mouse"
        assert result["S1"][2] == "mouse (GRCm39)"
        assert result["S1"][3] == ""  # no extended screening happened

    def test_includes_extended_top_hit_when_escalated(self, tmp_path):
        outPath = tmp_path / "lane"
        sampleDir = outPath / "FASTQC_Project_1_proj" / "Sample_S1"
        sampleDir.mkdir(parents=True)
        (sampleDir / "S1.rep").write_text(
            "94.0\t940\t940\tU\t0\tunclassified\n6.0\t60\t60\tS\t10090\tmouse\n"
        )
        (sampleDir / "S1.extended.krakenreport").write_text(
            "5.0\t50\t50\tU\t0\tunclassified\n95.0\t950\t950\tS\t3702\tarabidopsis\n"
        )

        result = buildContaminationDic(outPath, self._ssdf("S1"))

        assert result["S1"][3] == "arabidopsis"

    def test_empty_primary_report_yields_NA_row(self, tmp_path):
        outPath = tmp_path / "lane"
        sampleDir = outPath / "FASTQC_Project_1_proj" / "Sample_S1"
        sampleDir.mkdir(parents=True)
        (sampleDir / "S1.rep").write_text("")

        result = buildContaminationDic(outPath, self._ssdf("S1"))

        assert result["S1"] == ["NA", "None", "mouse (GRCm39)", ""]
