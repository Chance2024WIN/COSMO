# COSMO

Condition-specific Mapping of Operons (COSMO) is a tool for detecting operons in microorganisms

[![CircleCI](https://circleci.com/gh/hocinebendou/tb_operon_detection.svg?style=svg)](https://circleci.com/gh/hocinebendou/tb_operon_detection)

An **Operon** consists of a group of structural genes that codes for enzymes involved in a metabolic pathway.
The genes of an operon are located contiguously on a stretch DNA and are under the control of one promoter.
A single mRNA unit is transcribed and translated into separate proteins.

*COSMO* takes as input, a reference genome aligned BAM file and a GTF file with gene coordinates. 
COSMO predicts operons by calculating average coverages of the genes/CDSs and their intergenic regions (IGRs). 
The user provides four cut-offs which determine whether genes/CDSs form part of the same operon. 
If they fail any one of these cut-offs, then they are not part of the operon.

#### Requirements
1. Python >= 3.7
2. Pysam >=0.15.0

#### Installation (Run from a python script)
```
python setup.py install  
```

#### Usage (Run from command line)
```
user_input.py [-h] [-D GDEPTH] [-d IDEPTH] [-F GFACTOR] [-f IFACTOR]
                     [-p PREFIX]
                     ref length bam gtf

Detect possible genome operons using RNA expression coverages

positional arguments:
  ref                   Name of reference sequence in BAM file
  length                Length of reference genome according to GTF file
  bam                   Bam input file
  gtf                   Gtf input file

optional arguments:
  -h, --help            show this help message and exit
  -D GDEPTH, --gdepth GDEPTH
                        Average number of reads per base required to consider
                        a gene expressed
  -d IDEPTH, --idepth IDEPTH
                        Average number of reads per base required to consider
                        an IGR expressed
  -F GFACTOR, --gfactor GFACTOR
                        Maximum fold difference (FD) allowed between two genes/CDSs of an operon
  -f IFACTOR, --ifactor IFACTOR
                        Maximum FD allowed between an intergenic region (IGR) and its flanking genes/CDSs
  -o OUTPUT, --output OUTPUT
                        CSV output filename (output file is written to "output" folder)
```    
