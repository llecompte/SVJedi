#!/usr/bin/env python3

"""*******************************************************************************
    Name: SVjedi 
    Description: SVjedi aims to genotype structural variant with long reads data.
    Author: Lolita Lecompte
    Contact: lolita.lecompte@inria.fr, IRISA/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
    
    Copyright (C) 2019 Inria

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*******************************************************************************"""

import sys
import numpy as np
import argparse
import math
from decimal import *
from modules import entry
from modules import rules


def main(args):
    """ Parsing arguments """
    parser = argparse.ArgumentParser(description="Structural variations genotyping using long reads")
    
    parser.add_argument("-p", "--paf", metavar="<paffile>", nargs=1, help="PAF format", required=True)
    
    parser.add_argument("-v", "--vcf", metavar="<vcffile>", help="vcf format", required=True)
    
    parser.add_argument("-o", "--output", metavar="<output>", nargs=1, help="output file")
    
    parser.add_argument("-dover", metavar="<dist_overlap>", nargs=1, type=int, default=[100], help="breakpoint distance overlap")
    
    parser.add_argument("-dend", metavar="<dist_end>", nargs=1, type=int, default=[100], help="soft clipping length allowed for semi global alingments")

    parser.add_argument("-ladj", metavar="<allele_size>", nargs=1, type=int, default=[5000], help="Sequence allele adjacencies at each side of the SV")

    parser.add_argument(
        "-ms",
        "--minsupport",
        metavar="<minNbAln>",
        type=int,
        default=3,
        help="Minimum number of alignments to genotype a SV (default: 3>=)",
    )

    args = parser.parse_args()

    # check if outputfile
    if args.output is None:
        output = "genotype_results.txt"
    else:
        output = args.output[0]

    vcffile = args.vcf
    paffile = args.paf[0]
    min_support = args.minsupport
    d_over = args.dover[0]
    d_end = args.dend[0]
    l_adj = args.ladj[0]
    genotype(paffile, vcffile, output, min_support, d_over, d_end, l_adj)


def genotype(inputfile, vcf_without_gt, outputfile, min_aln, d_over, d_end, l_adj):
    """ Select alignment if it overlaps a junction and follow specific rules """ 
    dict_of_informative_aln = {}
    dict_for_ambiguous_reads = {}

    with open(inputfile) as file:
        for line in file:
            readId, readLength, readStart, readEnd, _, refId, refLength, refStart, refEnd, match, blockLength, quality, *_ = line.split("\t")
            #refSV = (refId.split("_")[1] + "_" + refId.split("_")[-1]) # ex: ref_2_48990527-636
            refSV = '_'.join(refId.split('_')[1:]) #ref_chrom_1_START_LENGTH > chrom_1_START_LENGTH
            svLength = abs(int(refId.split("-")[-1]))

            readLength, readStart, readEnd = (int(readLength), int(readStart), int(readEnd))
            refLength, refStart, refEnd = int(refLength), int(refStart), int(refEnd)

            aln = entry.alignment(readId, readLength, readStart, readEnd, refId, refLength, refStart, refEnd)
            
            # Quality filter
            if int(quality) < 10:
                continue
            # Overlapping filter : Test if read align on one of the junctions
            if (refStart + d_over) < l_adj < (refEnd - d_over) or (refStart + d_over) < (l_adj + svLength) < (refEnd - d_over):
                # Rules filter
                if rules.all_rules(aln, d_end):
                    fill_sv_dict(aln, dict_of_informative_aln)
                    
                    # Ambiguous read filter : remove redundant reads aligning on complementary ref
                    if readId not in list(dict_for_ambiguous_reads.keys()):
                        dict_for_ambiguous_reads[readId] = {refSV: []}
                        
                        if refId.startswith("r"):
                            dict_for_ambiguous_reads[readId][refSV].append(("r", int(quality)))
                        elif refId.startswith("d"):
                            dict_for_ambiguous_reads[readId][refSV].append(("d", int(quality)))

                        else:
                            if refSV not in list(dict_for_ambiguous_reads[readId].keys()):
                                dict_for_ambiguous_reads[readId] = {refSV: []}

                            if refId.startswith("r"):
                                dict_for_ambiguous_reads[readId][refSV].append(("r", int(quality)))
                            elif refId.startswith("d"):
                                dict_for_ambiguous_reads[readId][refSV].append(("d", int(quality)))

    # ambiguous filter
    for fragment in list(dict_for_ambiguous_reads.keys()):
        for region in list(dict_for_ambiguous_reads[fragment].keys()):
            if len(dict_for_ambiguous_reads[fragment][region]) > 1:
                best_quality = max(dict_for_ambiguous_reads[fragment][region], key=lambda x: x[1])
                fragments_to_remove = dict_for_ambiguous_reads[fragment][region]
                fragments_to_remove.remove(best_quality)

                for element in fragments_to_remove:
                    if element[0] == "d":
                        dict_of_informative_aln[region][1].remove(fragment)
                    elif element[0] == "r":
                        dict_of_informative_aln[region][0].remove(fragment)
    
    decision_vcf(dict_of_informative_aln, vcf_without_gt, outputfile, min_aln, l_adj)


def fill_sv_dict(a, dictReadAtJunction):
    """ If all conditions are true : overlapping and at least one rule is true, save the read """
    read = a.query
    reference = a.target
    #svId = reference.split("_")[1] + "_" + reference.split("_")[-1]
    svId = '_'.join(reference.split('_')[1:])
    
    if svId not in list(dictReadAtJunction.keys()):
        dictReadAtJunction[svId] = [[], []]
        if reference.startswith("r"):
            dictReadAtJunction[svId][0].append(read)
        else:
            dictReadAtJunction[svId][1].append(read)

    else:
        if reference.startswith("r"):
            dictReadAtJunction[svId][0].append(read)
        else:
            dictReadAtJunction[svId][1].append(read)



def encode_genotype(g):
    ''' Encode genotype from 0, 1, 2 to 0/0, 0/1, 1/1 '''
    if g == '0': genotype = "0/0"
    elif g == '1': genotype = "0/1"
    elif g == '2': genotype = "1/1"
    
    return genotype


def allele_normalization(nb_aln_per_allele, svtype, svlength, l_adj):
    ''' Allele length normalization '''
    if svlength > (2*l_adj): svlength = 2*l_adj #for upper bound, if case of sv size > 2XLadj, only 2 sequences of 2Ladj are represented
    
    if svtype == "DEL":
        nb_aln_longest_allele_seq = nb_aln_per_allele[0]
        if nb_aln_longest_allele_seq > 0:
            nb_aln_per_allele[0] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
    elif svtype == "INS":
        nb_aln_longest_allele_seq = nb_aln_per_allele[1]
        if nb_aln_longest_allele_seq > 0:
            nb_aln_per_allele[1] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
    return nb_aln_per_allele


def likelihood(all_count, svtype, svlength, minNbAln, l_adj):
    """ Compute likelihood """
    
    unbalanced_sv = ("DEL", "INS")
    if svtype in unbalanced_sv:
        c1, c2 = allele_normalization(all_count, svtype, svlength, l_adj)  # allelic sequence normalization for unbalanced SV
    else:
        c1, c2 = all_count
    
    rc1 = int(round(c1,0))
    rc2 = int(round(c2,0))
    e = 0.00005 #sequencing err

    lik0 = Decimal(c1*math.log10(1-e)) + Decimal(c2*math.log10(e)) 
    lik1 = Decimal((c1+c2)*math.log10(1/2)) 
    lik2 = Decimal(c2*math.log10(1-e)) + Decimal(c1*math.log10(e))
    
    L = [lik0, lik1, lik2]
    
    index_of_L_max = [i for i, x in enumerate(L) if x == max(L)]
    if len(index_of_L_max) == 1: 
        geno_not_encoded = str(index_of_L_max[0])
        geno = encode_genotype(geno_not_encoded)
        
    else:
        geno = "./."    #no genotype estimation since likelihood are not conclusive
    
    #Check for minimum number of alignment to assign a genotype
    if not sum(all_count) >= minNbAln: # minimum support
        geno = "./."
                
    combination = Decimal(math.factorial(rc1 + rc2)) / Decimal(math.factorial(rc1)) / Decimal(math.factorial(rc2))
    lik0 += combination 
    lik1 += combination
    lik2 += combination

    #phred scaled score
    prob0 = -10*lik0
    prob1 = -10*lik1
    prob2 = -10*lik2                                 
                 
    prob = [str(int(prob0)), str(int(prob1)), str(int(prob2))]
            
    return geno, prob
                    

def decision_vcf(dictReadAtJunction, inputVCF, outputDecision, minNbAln, l_adj):
    """ Output in VCF format and take genotype decision """
    getcontext().prec = 28
    outDecision = open(outputDecision, "w")
    with open(inputVCF) as inputFile:
        for line in inputFile:
            if line.startswith("##"):
                outDecision.write(line)

            elif line.startswith("#C"):
                outDecision.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                outDecision.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
                outDecision.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Depth of each allele by sample">\n')
                outDecision.write(line.rstrip("\n") + "\t" + "\t".join(["FORMAT", "SAMPLE"]) + "\n")

            else:
                in_chrom, in_start, _, __, in_type, ___, ____, in_info, *_ = line.rstrip("\n").split("\t")

                ### get SVTYPE ###
                if 'SVTYPE' in in_info:
                    if in_info.split(';')[-1].startswith('SVTYPE='):
                        svtype = in_info.split('SVTYPE=')[1]
                    else:
                        svtype = in_info.split('SVTYPE=')[1].split(';')[0]
                else:
                    svtype = ''
                
                ### get LENGTH for DELETION ###                 
                if svtype == 'DEL':
                    if "SVLEN=FALSE" in in_info:
                        if in_info.startswith("END="):
                            end = in_info.split("END=")[1].split(";")[0]
                        else:
                            end = in_info.split(";END=")[1].split(";")[0]

                        in_length = int(end) - int(in_start)

                    elif "SVLEN=" in in_info:
                        if in_info.split(';')[-1].startswith('SVLEN='):
                            in_length = abs(int(in_info.split("SVLEN=")[1]))
                        else:
                            in_length = abs(int(in_info.split("SVLEN=")[1].split(";")[0]))
                
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                    
                
                ### get LENGTH for INSERTION ###                    
                elif svtype == 'INS': 
                    in_length = len(in_type) 
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                
                
                ### get LENGTH for INVERSION ###                
                elif svtype == 'INV':
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(';')[0]
                    else:
                        end = in_info.split(";END=")[1].split(';')[0]               
                    in_length = int(end) - int(in_start)
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV


                ### get sv id for TRANSLOCATION ###             
                elif svtype == 'BND': 
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(';')[0]
                    elif ";END=" in in_info:
                        end = in_info.split(";END=")[1].split(';')[0]   
                    elif "[" in  in_type:
                        end = in_type.split(':')[1].split('[')[0]
                    elif "]" in  in_type:
                        end = in_type.split(':')[1].split(']')[0]
                    
                    if 'CHR2=' in in_info: 
                        chr2 = in_info.split('CHR2=')[1].split(';')[0]
                    elif '[' in in_type:
                            chr2 = in_type.split(':')[0].split('[')[1]     #ALT[CHR2:POSTION[
                    elif ']' in in_type:
                            chr2 = in_type.split(':')[0].split(']')[1]     #ALT]CHR2:POSTION]
                            
                    in_sv = in_chrom + "_" + in_start + "-" + chr2 + "-" + end #define sv id for TRANS
                
                
                #######################################################################################
                #Asign genotype 
                
                if svtype in ('DEL', 'INS', 'INV', 'BND') and in_sv in list(dictReadAtJunction.keys()) and abs(in_length) >= 50:
                    nbAln = [len(x) for x in dictReadAtJunction[in_sv]]
                    genotype, proba = likelihood(nbAln, svtype, in_length, minNbAln, l_adj)
            
                else: #if svtype different from DEL, INS, INV, BND or if sv not in supported by alignment
                    nbAln = [0,0]
                    genotype = "./."
                    proba = [".",".","."]
                
                    
                #######################################################################################
                #Output genotype in VCF
                
                numbers = ",".join(str(y) for y in nbAln)
                if len(line.split("\t")) <= 8:
                    new_line = (
                        line.rstrip("\n")
                        + "\t"
                        + "GT:DP:AD:PL"
                        + "\t"
                        + genotype
                        + ":"
                        + str(round(sum(nbAln), 3))
                        + ":"
                        + str(numbers)
                        + ":"
                        + str(','.join(proba))
                    )
                    outDecision.write(new_line + "\n")

                else:
                    line_without_genotype = line.split("\t")[0:8]
                    new_line = (
                        "\t".join(line_without_genotype)
                        + "\t"
                        + "GT:DP:AD:PL"
                        + "\t"
                        + genotype
                        + ":"
                        + str(round(sum(nbAln), 3))
                        + ":"
                        + str(numbers)
                        + ":"
                        + str(','.join(proba))
                    )
                    outDecision.write(new_line + "\n")

    outDecision.close()


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
