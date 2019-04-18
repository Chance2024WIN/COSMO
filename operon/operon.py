from dataclasses import dataclass, field
from typing import List

from .gene import Gene

@dataclass
class Operon:
    start: int
    end: int
    strand: str
    cov: float
    first: str
    last: str
    genes: List[Gene] = field(default_factory=list)
    num_genes: int = field(default=1)

    @classmethod
    def operon(cls, start, end, strand, gene_id, cov):
        op = cls(start, end, strand, cov, gene_id, '')
        op.genes.append(Gene(start, end, gene_id, cov, 0))

        return op

    def add_gene(self, start, end, gene_id, cov, igr):
        self.genes.append(Gene(start, end, gene_id, cov, igr))
        self.end = end
        self.cov = round((self.cov + cov + igr) / 3, 2)
        self.num_genes += 1
        self.last = gene_id

    def info_csv_row(self):
        row = [
            str(self.start),
            str(self.end),
            str(self.end - self.start),
            self.strand,
            str(self.cov),
            str(self.num_genes)
        ]

        if self.last == '':
            row.insert(0, self.first)
        else:
            row.insert(0, self.first + ' - ' + self.last)

        return row

    def gene_csv_rows(self):
        rows = []
        for gene in self.genes:
            rows.append(gene.gene_csv_row())

        return rows

    def last_gene_cov(self):
        """
        return the last gene coverage of the current operon
        """
        last_gene: Gene = self.genes[-1]
        return last_gene.cov

    def genes_covs(self):
        """
        Return a list of coverages of the operon genes
        """
        covs = []
        for gene in self.genes:
            covs.append(gene.cov)
        return covs