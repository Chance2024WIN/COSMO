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
