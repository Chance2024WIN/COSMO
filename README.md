# tb_operon_detection
A tool for detecting operons in Mycobacterium Tuberculosis 

[![CircleCI](https://circleci.com/gh/hocinebendou/tb_operon_detection.svg?style=svg)](https://circleci.com/gh/hocinebendou/tb_operon_detection)

An **Operon** consists of a group of structural genes that codes for enzymes involved in a metabolic pathway.
The genes of an operon are located contiguously on a stretch DNA and are under the control of one promoter.
A single mRNA unit is transcribed and translated into separate proteins.

*tb_operon_detection* takes as input an aligned BAM file to TB reference genome, a GTF file with gene's 
coordinates and predicts operons by calculating average coverages of the genes and their intergenic regions 
in the reference genome. The ratios of the coverages are compared if they don't exceed a threshold value 
given as input. The genes then are part of same operon. If not the genes are in separate operons.

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
  length                Length of reference sequence
  bam                   Bam input file
  gtf                   Gtf input file

optional arguments:
  -h, --help            show this help message and exit
  -D GDEPTH, --gdepth GDEPTH
                        Average number of reads per base required to consider
                        a gene expressed
  -d IDEPTH, --idepth IDEPTH
                        Average number of reads per base required to consider
                        a IGR expressed
  -F GFACTOR, --gfactor GFACTOR
                        Allowed difference factor of two gene's coverages to
                        be part of same Operon
  -f IFACTOR, --ifactor IFACTOR
                        Allowed difference factor of IGR and adjacent gene's
                        coverages to be part of same operon
  -o OUTPUT, --output OUTPUT
                        CSV output filename under output folder
```    
