from operon.gene import Gene


def test_gene_csv_row():
    start = 1
    end = 100
    gene_id = 'i001x'
    cov = 75.23
    igr_cov = 102.31

    gene = Gene(start, end, gene_id, cov, igr_cov)
    row = gene.gene_csv_row()
    assert len(row) == 10
    for x in range(7):
        assert row[x] == ''
    assert row[7] == 'i001x'
    assert row[8] == str(75.23)
    assert row[9] == str(102.31)
