# SVJedi : Genotyping structural variations with long read data

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/svjedi/README.html)

SVJedi is a structural variation (SV) genotyper for long read data. 
Based on a representation of the different alleles, it estimates the genotype of each variant in a given individual sample based on allele-specific alignment counts.
SVJedi takes as input a *variant file* (VCF), a *reference genome* (fasta) and a *long read file* (fasta/fastq) and 
outputs the initial variant file with an additional column containing genotyping information (VCF).

SVJedi processes **deletions**, **insertions**, **inversions** and **translocations**.

SVJedi is organized in three main steps:

1. Generate representative allele sequences of a set of SVs given in a vcf file
2. Map reads on previously generated allele sequences using Minimap2
3. Genotype SVs and output a vcf file

*Jedi comes from the verb jediñ* ['ʒeːdɪ] *in Breton, it means calculate.*


### Requirements

- Python3
- Minimap2 (we recommend using version 2.17-r941)
- NumPy
- Biopython


### Usage

    python3 svjedi.py -v <set_of_sv.vcf> -r <refgenome.fasta> -i <long_reads.fastq>

Note: Chromosome names in `reference.fasta` and in `set_of_sv.vcf` must be the same. 
Also, the `SVTYPE` tag must be present in the VCF (`SVTYPE=DEL` or `SVTYPE=INS` or `SVTYPE=INV` or `SVTYPE=BND`). 
More details are given in [SV representation in VCF](#SV-representation-in-VCF).


### Installation

    git clone https://github.com/llecompte/SVJedi.git

SVJedi is also distributed as a [Bioconda package](https://anaconda.org/bioconda/svjedi):

	conda install -c bioconda svjedi	

### Examples

The folder Data/HG002_son includes an example of 20 SVs (10 insertions and 10 deletions) to genotype on a subsample of a real human dataset of the Ashkenazim son HG002.

Example command line:

	python3 svjedi.py -v Data/HG002_son/HG002_20SVs_Tier1_v0.6_PASS.vcf -a Data/HG002_son/reference_at_breakpoints.fasta -i Data/HG002_son/PacBio_reads_set.fastq.gz -o Data/HG002_son/genotype_results.vcf

Note: Genotyping results in `Data/HG002_son/expected_genotype_results.vcf` were obtained using Minimap2 version 2.17-r941.
 
 

The folder Data/C_elegans includes an example on 12 SVs (del, ins, inv, bnd) to genotype with a small synthetic read dataset on a subset of the *Caenorhabditis elegans* genome.

Example command line:

    python3 svjedi.py -v Data/C_elegans/test.vcf -r Data/C_elegans/genome.fasta -i Data/C_elegans/simulated-reads.fastq.gz



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
| -ladj	       | Length of sequences adjacent to each end of breakpoints (default 5,000 bp) |
| -d/--data    | Type of sequencing data, either *ont* or *pb* (default pb)  |
| -t/--threads | Number of threads for mapping             |
| -h/--help    | Show help                                 |

### SV representation in VCF

Here are the information needed for SVJedi to genotype the following SV types. All variants must have the ```CHROM``` and ```POS``` fields defined, with the chromosome names in `reference.fasta` and in `set_of_sv.vcf` that must be the same. Then additional information is required according to SV type:

- Deletion
	- Either ```ALT``` field is ```<DEL>``` or ```INFO``` field must contain ```SVTYPE=DEL```
	- ```INFO``` field must contain either ```END=pos``` (with `pos` being the end position of the deleted segment) or ```SVLEN=len``` (with `len` being the size of the deletion) tags

- Insertion
	- ```INFO``` field must contain ```SVTYPE=INS```
	- ```ALT``` field must contain the sequence of the insertion 
	
- Inversion
	- Either ```ALT``` field is ```<INV>``` or ```INFO``` field must contain ```SVTYPE=INV```
	- ```INFO``` field must contain ```END=pos``` tag, with `pos` being the second breakpoint position

- Translocation
	- ```INFO``` field must contain ```SVTYPE=BND``` and ```CHR2=``` and ```END=``` tags
	- CHR2 name and sequence must be in the reference genome fasta file
	- ```ALT``` field must be formated as: ```t[chr:pos[```, ```t]chr:pos]```, ```]chr:pos]t``` or ```[chr:pos[t```, with `chr`and `pos` indicating the second breakpoint position and brackets directions indicating which parts of the two chromosomes should be joined together 

### Reference

SVJedi: Genotyping structural variations with long reads. Lecompte L, Peterlongo P, Lavenier D, Lemaitre C. *Bioinformatics* 2020, 36(17):4568–4575 [doi:10.1093/bioinformatics/btaa527](https://academic.oup.com/bioinformatics/article/36/17/4568/5841661) ([bioRxiv preprint](https://www.biorxiv.org/content/10.1101/849208v1))


### Contact

SVJedi is a [Genscale](http://team.inria.fr/genscale/) tool developed by Lolita Lecompte. For any bug report or feedback, please use the Issues form of Github. 

