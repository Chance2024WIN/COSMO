from operon.operon import Operon


def test_operon_constructor():
    start = 1
    end = 100
    strand = '+'
    gene_id = 'i001x'
    cov = 75.23
    op = Operon.operon(start, end, strand, gene_id, cov)
    assert len(op.genes) == 1
    assert op.genes[0].gene_id == 'i001x'


def test_add_gene():
    op = Operon.operon(1, 100, '+', 'i001x', 75.23)
    op.add_gene(101, 200, 'i002x', 63.47, 102.03)
    assert len(op.genes) == 2
    assert op.genes[1].gene_id == 'i002x'


def test_info_csv_row():
    op = Operon.operon(1, 100, '+', 'i001x', 75.23)
    op.add_gene(101, 200, 'i002x', 63.47, 102.03)
    row = op.info_csv_row()
    assert row[0] == 'i001x - i002x'
    assert int(row[1]) == 1
    assert int(row[2]) == 200
    assert int(row[3]) == 199
    assert row[4] == '+'
    assert round(float(row[5]),2) == 80.24
    assert int(row[6]) == 2


def test_gene_csv_rows():
    op = Operon.operon(1, 100, '+', 'i001x', 75.23)
    op.add_gene(101, 200, 'i002x', 63.47, 102.03)
    rows = op.gene_csv_rows()
    assert len(rows) == 2
    assert len(rows[0]) == 10
    assert rows[0][0] == ''
    assert rows[0][7] == 'i001x'
    assert rows[0][8] == str(75.23)
    assert rows[0][9] == ''
    assert rows[1][7] == 'i002x'
    assert rows[1][8] == str(63.47)
    assert rows[1][9] == str(102.03)


def test_last_gene_cov():
    op = Operon.operon(1, 100, '+', 'i001x', 75.23)
    op.add_gene(101, 200, 'i002x', 63.47, 102.03)
    last_gene = op.last_gene_cov()
    assert last_gene == 63.47


def test_genes_covs():
    op = Operon.operon(1, 100, '+', 'i001x', 75.23)
    op.add_gene(101, 200, 'i002x', 63.47, 102.03)
    covs = op.genes_covs()
    assert len(covs) == 2
    assert covs[0] == 75.23
    assert covs[1] == 63.47
