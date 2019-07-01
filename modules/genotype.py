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
import argparse
from modules import entry
from modules import rules


def main(args):
    """ Parsing arguments """
    parser = argparse.ArgumentParser(description="Structural variations genotyping using long reads")
    
    parser.add_argument("-p", "--paf", metavar="<paffile>", nargs=1, help="PAF format", required=True)
    
    parser.add_argument("-v", "--vcf", metavar="<vcffile>", help="vcf format", required=True)
    
    parser.add_argument("-o", "--output", metavar="<output>", nargs=1, help="output file")
    
    parser.add_argument("-dover", metavar="<dist_overlap>", nargs=1, type=int, default=100, help="breakpoint distance overlap")
    
    parser.add_argument("-dend", metavar="<dist_end>", nargs=1, type=int, default=100, help="soft clipping length allowed for semi global alingments")

    parser.add_argument(
        "-ms",
        "--minsupport",
        metavar="<minNbAln>",
        type=int,
        default=3,
        help="Minimum number of alignments to consider a deletion (default: 3>=)",
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
    d_over = args.dover
    d_end = args.dend
    genotype(paffile, vcffile, output, min_support, d_over, d_end)


def genotype(inputfile, vcf_without_gt, outputfile, min_aln, d_over, d_end):
    """ Select alignment if it overlaps a junction and follow specific rules """ 
    dict_of_informative_aln = {}
    dict_for_ambiguous_reads = {}

    with open(inputfile) as file:
        for line in file:
            readId, readLength, readStart, readEnd, _, refId, refLength, refStart, refEnd, match, blockLength, quality, *_ = line.split("\t")
            refSV = (refId.split("_")[1] + "_" + refId.split("_")[-1]) # ex: ref_2_48990527-636
            deletionLength = abs(int(refId.split("-")[-1]))

            readLength, readStart, readEnd = (int(readLength), int(readStart), int(readEnd))
            refLength, refStart, refEnd = int(refLength), int(refStart), int(refEnd)

            aln = entry.alignment(readId, readLength, readStart, readEnd, refId, refLength, refStart, refEnd)
            
            # Quality filter
            if int(quality) < 10:
                continue
            # Overlapping filter : Test if read align on one of the junctions
            if (refStart + d_over) < 5000 < (refEnd - d_over) or (refStart + d_over) < (5000 + deletionLength) < (refEnd - d_over):
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
    
    decision_vcf(dict_of_informative_aln, vcf_without_gt, outputfile, min_aln)


def fill_sv_dict(a, dictReadAtJunction):
    """ If all conditions are true : overlapping and at least one rule is true, save the read """
    read = a.query
    reference = a.target
    svId = reference.split("_")[1] + "_" + reference.split("_")[-1]

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


def decision_vcf(dictReadAtJunction, inputVCF, outputDecision, minNbAln):
    """ Output in VCF format and take genotype decision """
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

                if in_type != "<DEL>":
                    continue

                if "SVLEN=FALSE" in in_info:
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(";")[0]
                    else:
                        end = in_info.split(";END=")[1].split(";")[0]

                    in_length = abs(int(end) - int(start))

                elif "SVLEN=" in in_info:
                    in_length = abs(int(in_info.split("SVLEN=")[1].split(";")[0]))

                if abs(in_length) < 50:
                    continue

                in_sv = in_chrom + "_" + in_start + "-" + str(in_length)
                if in_sv not in list(dictReadAtJunction.keys()):
                    continue

                nbAln = [len(x) for x in dictReadAtJunction[in_sv]]

                # normalization
                nbAln_prev = nbAln[0]
                if nbAln[0] != 0:
                    nbAln[0] = round(nbAln_prev * 10000 / (10000 + in_length), 3)

                # genotype estimation
                if (0.2 <= (nbAln[0] / (nbAln[0] + nbAln[1])) <= 0.8):
                    geno = "0/1"
                elif nbAln[0] > nbAln[1]:
                    geno = "0/0"
                elif nbAln[0] < nbAln[1]:
                    geno = "1/1"

                # minimum support
                if sum(nbAln) >= minNbAln:
                    numbers = ",".join(str(y) for y in nbAln)
                    if len(line.split("\t")) <= 8:
                        new_line = (
                            line.rstrip("\n")
                            + "\t"
                            + "GT:DP:AD"
                            + "\t"
                            + geno
                            + ":"
                            + str(round(sum(nbAln), 3))
                            + ":"
                            + str(numbers)
                        )
                        outDecision.write(new_line + "\n")

                    else:
                        line_without_genotype = line.split("\t")[0:8]
                        new_line = (
                            "\t".join(line_without_genotype)
                            + "\t"
                            + "GT:DP:AD"
                            + "\t"
                            + geno
                            + ":"
                            + str(round(sum(nbAln), 3))
                            + ":"
                            + str(numbers)
                        )
                        outDecision.write(new_line + "\n")
    outDecision.close()


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
