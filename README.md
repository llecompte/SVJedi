# SVJedi : Genotyping structural variations with long read data

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

SVJedi is a structural variation (SV) genotyper for long read data. 
Based on a representation of the different alleles, it estimates the genotype of each variant from specific alignements obtained.
SVJedi takes as input a *variant file* (VCF), a *reference genome* (fasta) and a *long read file* (fasta/fastq) and 
outputs the initial variant file with an additional column containing genotyping information (VCF).

For the moment, SVJedi process **deletions** and **insertions**.

SVJedi is organized in three main steps:

1. Generate allele reference sequences of a set of deletions given in a vcf file
2. Map reads on previously generated references using Minimap2
3. Genotype deletions and output a vcf

*Jedi comes from the verb jediñ* ['ʒeːdɪ] *in Breton, it means calculate.*


### Requirements

- Python3
- Minimap2
- numpy


### Usage

    python3 svjedi.py -v <set_of_sv.vcf> -r <reference.fasta> -i <long_reads.fastq>
    
Note: Chromosome names in `reference.fasta` and in `set_of_sv.vcf` must be the same. 
Also, the `SVTYPE` tag must be present in the VCF (`SVTYPE=DEL` or as `SVTYPE=INS`).


### Installation

    git clone https://github.com/llecompte/SVJedi.git

### Example

The folder Data includes an example on 5 deletions to genotype with a small synthetic read dataset on *Caenorhabditis elegans*.

Example command line:

    python3 svjedi.py -v Data/5del_C_elegans.vcf -r Data/C_elegans.chr.I.fa -i Data/reads_C_elegans.fastq.gz -o 5del_C_elegans_with_genotypes.vcf
    


### Parameters

SVJedi two different usages from non aligned reads or from aligned reads (PAF format).

```
    python3 svjedi.py -v <set_of_sv.vcf> -r <reference.fasta> -i <long_reads.fastq>
    
    python3 svjedi.py -v <set_of_sv.vcf> -p <alignments.paf>
```

| Option       | Description                               |
| ------------ | ----------------------------------------- | 
| -v/--vcf     | set of SVs in VCF                   |
| -r/--ref     | reference genome in FASTA                 |
| -i/--input   | sequenced long reads in FASTQ or FASTQ.GZ (1 file or multiple files)|
| -p/--paf     | alignments in PAF                         |
| -o/--output  | output file with genotypes in VCF                |
| -ms/--minsupport | minimum number of informative alignments to assign a genotype
| -dover       | Breakpoint distance overlap required (default 100 bp) |
| -dend        | Soft-clipping length allowed to consider a semi-global alignment (default 100 bp) |
| -d/--data    | type of sequencing data, either *ont* or *pb* (default pb)  |
| -t/--threads | number of threads for mapping             |
| -h/--help    | show help                                 |


### Contact

<lolita.lecompte@inria.fr>
