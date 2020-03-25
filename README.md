# SVJedi : Genotyping structural variations with long read data

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

SVJedi is a structural variation (SV) genotyper for long read data. 
Based on a representation of the different alleles, it estimates the genotype of each variant in a given individual sample based on allele-specific alignment counts.
SVJedi takes as input a *variant file* (VCF), a *reference genome* (fasta) and a *long read file* (fasta/fastq) and 
outputs the initial variant file with an additional column containing genotyping information (VCF).

SVJedi processes **deletions**, **insertions**, **inversions** and **translocations**.

SVJedi is organized in three main steps:

1. Generate allele reference sequences of a set of SVs given in a vcf file
2. Map reads on previously generated references using Minimap2
3. Genotype SVs and output a vcf file

*Jedi comes from the verb jediñ* ['ʒeːdɪ] *in Breton, it means calculate.*


### Requirements

- Python3
- Minimap2
- NumPy
- Biopython


### Usage

    python3 svjedi.py -v <set_of_sv.vcf> -r <refgenome.fasta> -i <long_reads.fastq>

Note: Chromosome names in `reference.fasta` and in `set_of_sv.vcf` must be the same. 
Also, the `SVTYPE` tag must be present in the VCF (`SVTYPE=DEL` or as `SVTYPE=INS` or as `SVTYPE=INV` or as `SVTYPE=BND`). 
More details are given in [SV representation in VCF](#SV-representation-in-VCF).


### Installation

    git clone https://github.com/llecompte/SVJedi.git

### Example

The folder Data/HG002_son includes an example of 20 SVs (10 insertions and 10 deletions) to genotype on a subsample of a real human dataset of the Ashkenazim son HG002.

Example command line:

	python3 svjedi.py -v Data/HG002_son/HG002_20SVs_Tier1_v0.6_PASS.vcf -a Data/HG002_son/reference_at_breakpoints.fasta -i Data/HG002_son/PacBio_reads_set.fastq.gz -o Data/HG002_son/genotype_results.vcf


The folder Data/C_elegans includes an example on 5 deletions to genotype with a small synthetic read dataset on *Caenorhabditis elegans*.

Example command line:

    python3 svjedi.py -v Data/C_elegans/5del_C_elegans.vcf -r Data/C_elegans/C_elegans.chr.I.fa -i Data/C_elegans/reads_C_elegans.fastq.gz -o Data/C_elegans/5del_C_elegans_with_genotypes.vcf



### Parameters

SVJedi two different usages from non aligned reads or from aligned reads (PAF format).

```
    python3 svjedi.py -v <set_of_sv.vcf> -r <refgenome.fasta> -i <long_reads.fastq>
    
    python3 svjedi.py -v <set_of_sv.vcf> -a <refallele.fasta> -i <long_reads.fastq>
    
    python3 svjedi.py -v <set_of_sv.vcf> -p <alignments.paf>
```

| Option       | Description                               |
| ------------ | ----------------------------------------- |
| -v/--vcf     | Set of SVs in VCF                         |
| -r/--ref     | Reference genome in FASTA                 |
| -i/--input   | Sequenced long reads in FASTQ or FASTQ.GZ (1 file or multiple files)|
| -a/--allele  | Reference sequences of alleles            |
| -p/--paf     | Alignments in PAF                         |
| -o/--output  | Output file with genotypes in VCF         |
| -ms/--minsupport | Minimum number of informative alignments to assign a genotype
| -dover       | Breakpoint distance overlap required (default 100 bp) |
| -dend        | Soft-clipping length allowed to consider a semi-global alignment (default 100 bp) |
| -d/--data    | Type of sequencing data, either *ont* or *pb* (default pb)  |
| -t/--threads | Number of threads for mapping             |
| -h/--help    | Show help                                 |

### SV representation in VCF

Some information needed for SVJedi to genotype the following variants. All variants must have the ```CHROM``` and ```POS``` fields defined. Then additional information is required according to SV type:

- Deletion
	- Either ```ALT``` field is ```<DEL>``` or ```INFO``` field must contain ```SVTYPE=DEL```
	- ```INFO``` field must contain either ```END=``` or ```SVLEN=``` tag

- Insertion
	- ```INFO``` field must contain ```SVTYPE=INS```
	- ```ALT``` field must contain the sequence of the insertion 
	
- Inversion
	- Either ```ALT``` field is ```<INV>``` or ```INFO``` field must contain ```SVTYPE=INV```
	- ```INFO``` field must contain ```END=``` tag

- Translocation
	- ```INFO``` field must contain ```SVTYPE=BND``` and ```CHR2=``` and ```END=``` tags
	- CHR2 name and sequence must be in reference genome
	- ```ALT``` field ALT must be formated as VCF Specification: ```t[chr:pos[```, ```t]chr:pos]```, ```]chr:pos]t``` or ```[chr:pos[t```

### Contact

SVJedi is a [Genscale](http://team.inria.fr/genscale/) tool developed by Lolita Lecompte <lolita.lecompte@inria.fr>

