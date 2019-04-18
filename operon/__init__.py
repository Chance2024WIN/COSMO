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
        reads = alignment_file.fetch(start=start, stop=end, region=ref)
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
                    # print(num_bases)

        average = num_bases / (end - start)

    return round(average, 2)


class Analyze(ABC):
    """
    Classes inheriting this class should implement the 'analyze' method
    """
    @abstractmethod
    def analyze(self):
        pass
