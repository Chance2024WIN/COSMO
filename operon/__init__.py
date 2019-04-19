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
from abc import ABC, abstractmethod


def average_coverage(alignment_file, ref, start, end, strand):
    """
    Calculate average coverage of reads in contig region located in
    one strand and defined by the limits start and end.
    """
    average = 0
    if start < end:
        # region = ref + ':' + str(start) + '-' + str(end)
        # reads = alignment_file.fetch(region=region)
        # reads = alignment_file.fetch(start=start, stop=end, region=ref)
        reads = alignment_file.fetch(ref, start=start, stop=end)
        num_bases = 0
        for read in reads:
            if not read.is_secondary and \
               not read.is_supplementary and \
               not read.is_unmapped:
                if (not read.is_reverse and strand == "+") or \
                   (read.is_reverse and strand == "-"):
                    read_len = read.reference_length
                    read_start = read.reference_start
                    read_end = read.reference_end

                    if read_start < start and read_end > end:
                        num_bases += end - start
                    elif read_start < start < read_end:
                        num_bases += read_len - (start - read_start)
                    elif read_end > end > read_start:
                        num_bases += end - read_start
                    else:
                        num_bases += read_len

        average = num_bases / (end - start)

    return round(average, 2)


class Analyze(ABC):
    """
    Classes inheriting this class should implement the 'analyze' method
    """
    @abstractmethod
    def analyze(self):
        pass
