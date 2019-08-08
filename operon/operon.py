# Copyright (C) 2019  Hocine Bendou <hocine@sanbi.ac.za>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>
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
            str(self.op_avg_cov()),
            #str(self.cov),
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

    def op_avg_cov(self):
        count = 0
        j = 0
        eq = ""
        for i, gene in enumerate(self.genes):
            if i == 0:
                j += 1
                count += gene.cov
                eq += str(gene.cov)
            if i > 0:
                if gene.igr_cov > 0:
                    eq += ' + ' + str(gene.igr_cov)
                    count += gene.igr_cov
                    j+=1
                count += gene.cov
                eq += ' + ' + str(gene.cov)
                j += 1

        return  round(count / j, 2)
