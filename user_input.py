from operon.detect import Detect
import argparse
import sys


def run_command(args: argparse.ArgumentParser):
    """
    Execute user's input command and return detected operons in csv files
    """
    detect = Detect(args.ref, args.length, args.bam, args.gtf)
    detect.analyze()

    
def main():
    parser = argparse.ArgumentParser(description="Detect possible genome operons using RNA expression coverages")
    parser.add_argument("-r", "--ref", help="Name of reference sequence in BAM file", type=str)
    parser.add_argument("-l", "--length", help="Length of reference sequence", type=int)
    parser.add_argument("-b", "--bam", help="Bam input file", type=argparse.FileType())
    parser.add_argument("-g", "--gtf", help="Gtf input file")
    parser.add_argument("-D", "--gdepth", help="Average number of reads per base required to consider a gene expressed", type=int, default=10)
    parser.add_argument("-d", "--idepth", help="Average number of reads per base required to consider a IGR expressed", type=int, default=10)
    parser.add_argument("-F", "--gfactor", help="Allowed difference factor of two gene's coverages to be part of "
                                                "same Operon", type=float, default=4.5)
    parser.add_argument("-f", "--ifactor", help="Allowed difference factor of IGR and adjacent gene's coverages "
                                                "to be part of same operon", type=float, default=5.0)
    parser.add_argument("-p", "--prefix", help="Prefix output files", type=str, default="out-")

    args = parser.parse_args()
    run_command(args)

if __name__ == "__main__":
    main()
