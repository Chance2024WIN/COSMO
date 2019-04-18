from dataclasses import dataclass

@dataclass
class Gene:
    start: int
    end: int
    gene_id: str
    cov: float
    igr_cov: float

    def gene_csv_row(self):
        row = ['' for _ in range(7)]
        row.append(self.gene_id)
        row.append(str(self.cov))
        if self.igr_cov > 0:
            row.append(str(self.igr_cov))
        else:
            row.append('')

        return row
