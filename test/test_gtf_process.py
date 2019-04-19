from operon.gtf_process import GtfProcess
import os.path


def test_parse_line():
    cwd = os.path.dirname(__file__)
    filename = os.path.join(cwd, "data/test_file.gtf")
    gtf_process = GtfProcess(filename)
    lines = gtf_process.lines()
    gene_1 = lines.__next__()
    assert gene_1['gene_id'] == "Rv0166"
    gene_2 = lines.__next__()
    assert gene_2['gene_id'] == "Rv0167"
