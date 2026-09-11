from tools.prep_contaminome import dumpname, dumpnode


class Test_dumpnode:
    def test_writes_ncbi_nodes_dmp_format(self, tmp_path):
        taxmap = {
            "root": [1, 1, "no rank"],
            "human": [9606, 9, "species"],
        }
        ofile = tmp_path / "nodes.dmp"

        dumpnode(taxmap, ofile)

        lines = ofile.read_text().splitlines()
        assert lines == [
            "1\t|\t1\t|\tno rank\t|\t-\t|",
            "9606\t|\t9\t|\tspecies\t|\t-\t|",
        ]


class Test_dumpname:
    def test_writes_ncbi_names_dmp_format(self, tmp_path):
        taxmap = {
            "root": [1, 1, "no rank"],
            "human": [9606, 9, "species"],
        }
        ofile = tmp_path / "names.dmp"

        dumpname(taxmap, ofile)

        lines = ofile.read_text().splitlines()
        assert lines == [
            "1\t|\troot\t|\t-\t|\tscientific name\t|",
            "9606\t|\thuman\t|\t-\t|\tscientific name\t|",
        ]
