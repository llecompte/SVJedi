#!/bin/bash

python3 ../../svjedi.py -v test.vcf -r genome.fasta -i simulated-reads.fastq.gz

#Test if output as expected
md5sum -c <<<"e0bcc73775f1b22921ea2ff29991af15 genotype_results.vcf"
