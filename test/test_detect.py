import pysam

from operon import average_coverage
from operon.detect import Detect
from argparse import Namespace
import os

# Add your bam and bai files to data folder and uncomment this method to test
# detect operon analyze function.

# def test_detect():
#     cwd = os.path.dirname(__file__)
#     args = Namespace()
#     setattr(args, 'ref', 'NC_000962')
#     setattr(args, 'length', 4411532)
#     setattr(args, 'bam', os.path.join(cwd, 'data/test_file.bam'))
#     setattr(args, 'gtf', os.path.join(cwd, 'data/test_file.gtf'))
#     detect = Detect(args.ref, args.length, args.bam, args.gtf)
#     assert detect.ref_seq == 'NC_000962'
#     assert detect.gene_depth == 10
#     # detect.analyze()
#     alignment_file = pysam.AlignmentFile(detect.bam_file, 'rb')
#     cov = average_coverage(alignment_file, detect.ref_seq, 75301, 76212, '+')
#     assert round(cov,2) == 9.5
#     cov = average_coverage(alignment_file, detect.ref_seq, 80624, 81673, '+')
#     assert round(cov,2) == 0.0
