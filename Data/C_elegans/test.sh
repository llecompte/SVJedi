#!/bin/bash

python3 ../../svjedi.py -v test.vcf -r reference.fasta -i simulated-reads.fastq.gz

#Test if output as expected
md5sum -c <<<"0157c1baf07daf604e71395284849834 genotype_results.vcf"
