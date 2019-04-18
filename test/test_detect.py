import pysam

from operon import average_coverage
from operon.detect import Detect
from argparse import Namespace
import os

from operon.gtf_process import GtfProcess


def test_detect():
    cwd = os.path.dirname(__file__)
    args = Namespace()
    setattr(args, 'ref', 'NC_000962')
    setattr(args, 'length', 4411532)
    setattr(args, 'bam', os.path.join(cwd, 'data/test_file.bam'))
    setattr(args, 'gtf', os.path.join(cwd, 'data/test_file.gtf'))
    detect = Detect(args.ref, args.length, args.bam, args.gtf)
    assert detect.ref_seq == 'NC_000962'
    assert detect.gene_depth == 10
    detect.analyze()
    alignment_file = pysam.AlignmentFile(detect.bam_file, 'rb')
    gtf_process = GtfProcess(detect.gtf_file)
    lines = gtf_process.lines()
    line = lines.__next__()
    cov = average_coverage(alignment_file, detect.ref_seq, int(line['start']), int(line['end']), line['strand'])
    assert round(cov,2) == 107.20
    line = lines.__next__()
    cov = average_coverage(alignment_file, detect.ref_seq, int(line['start']), int(line['end']), line['strand'])
    assert round(cov, 2) == 107.20
    alignment_file.close()
