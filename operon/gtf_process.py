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
# Adapted from: https://gist.github.com/slowkow/8101481
from dataclasses import dataclass, field
import gzip
import re

GTF_HEADER = ['seq_name', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA = re.compile(r'\s*,\s*')
R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')

@dataclass
class GtfProcess:
    """
    Read GTF file and return a generator of processed lines from the GTF rows.
    """
    filename: str = field(repr=False)

    def lines(self):
        """
        Open an optionally gzipped GTF file and generate a dict for each line.
        """
        fn_open = gzip.open if self.filename.endswith('.gz') else open

        with fn_open(self.filename) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                else:
                    yield self.parse(line)

    def parse(self, line):
        """
        Parse a single GTF line and return a dict.
        """
        result = {}

        fields = line.rstrip().split('\t')

        for i, col in enumerate(GTF_HEADER):
            result[col] = self._get_value(fields[i])

        # INFO field consists of "key1=value;key2=value;...".
        infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

        for i, info in enumerate(infos, 1):
            try:
                key, _, value = re.split(R_KEYVALUE, info, 1)
            except ValueError:
                key = f"INFO{i}"
                value = info
            if value:
                result[key] = self._get_value(value)

        return result

    @classmethod
    def _get_value(cls, value):
        if not value:
            return None

        # Strip double and single quotes.
        value = value.strip('"\'')

        # Return a list if the value has a comma.
        if ',' in value:
            value = re.split(R_COMMA, value)

        # These values are equivalent to None.
        elif value in ['', '.', 'NA']:
            return None

        return value