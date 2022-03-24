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
from Bio.Seq import Seq


def main(args):
    """ Parsing arguments """
    parser = argparse.ArgumentParser(
        description="Structural variations genotyping using long reads"
    )
    parser.add_argument(
        "-v", "--vcf", metavar="<vcffile>", nargs=1, help="vcf format", required=True
    )
    parser.add_argument(
        "-r", "--ref", metavar="<reffile>", nargs=1, help="fasta format", required=True
    )

    parser.add_argument("-ladj", metavar="<allele_size>", nargs=1, type=int, default=[5000], help="Sequence allele adjacencies at each side of the SV")

    args = parser.parse_args()
    
    genome_file = args.ref[0]
    vcf_file = args.vcf[0]
    l_adj = args.ladj[0]
    create_ref(genome_file, vcf_file, l_adj)


def create_ref(genome, set_of_sv, l_adj):
    """ Generate triplet of reference """
    dict_of_chrom = {}
    list_of_deletions = []
    list_of_insertions = []
    list_of_inversions = []
    list_of_translocations = []

    # original sequence ref required
    with open(genome) as sequenceFile:
        sequence = ""
        for line in sequenceFile:
            if line.startswith(">"):
                if sequence != "":
                    dict_of_chrom[header] = sequence
                sequence = ""
                header = line.rstrip("\n")[1:].split()[0]

            else:
                sequence += line.rstrip("\n")

        # save last one
        dict_of_chrom[header] = sequence

    previous_start = 0
    previous_end = 0

    with open(set_of_sv) as svFile:
        for line in svFile:
            if not line.startswith("#"):
                if len(line.split("\t")) > 6:
                    chrom, start, _, __, alt_field, ___, ____, info, *_ = line.rstrip(
                        "\n"
                    ).split("\t")
                else:
                    chrom, start, _, __, alt_field, ___, ____, info = line.rstrip(
                        "\n"
                    ).split("\t")
                
                if 'SVTYPE' in info:
                    svtype = info.split('SVTYPE=')[1].split(';')[0]

                else:
                    svtype = '' 
                
                #for deletions
                if alt_field == "<DEL>" or svtype == 'DEL':
                
                    start = int(start)
                    if info.split(';')[-1].startswith('SVLEN='):
                        length = abs(int(info.split("SVLEN=")[1]))
                    else:
                        length = abs(int(info.split("SVLEN=")[1].split(";")[0]))
                    end = start + length

                    # focus on >=50bp length deletion
                    if length >= 50:
                        list_of_deletions.append((chrom, start, end, length))
                
                           
                #for insertion
                elif svtype == 'INS':
                    start = int(start)
                    length = len(alt_field)
                                
                    # focus on >=50bp length insertion
                    if length >= 50:
                        list_of_insertions.append((chrom, start, length, alt_field))

                #for inversions
                elif alt_field == '<INV>' or svtype == 'INV':
                    start = int(start)
                    if info.startswith('END='):
                        end = int(info.split('END=')[1].split(';')[0])
                    elif ';END=' in info:
                        end = int(info.split(';END=')[1].split(';')[0])
                    else:
                        continue
                    length = end - start
                     
                    # focus on >=50bp length inversions
                    if length >= 50:
                        list_of_inversions.append((chrom, start, end, length))

                #for translocations
                elif svtype == 'BND':
                    start = int(start)
                    if 'CHR2' in info: 
                        chrom2 = info.split('CHR2=')[1].split(';')[0]
                    elif '[' in alt_field:
                            chrom2 = alt_field.split(':')[0].split('[')[1]     #ALT[CHR2:POSTION[
                    elif ']' in alt_field:
                            chrom2 = alt_field.split(':')[0].split(']')[1]     #ALT]CHR2:POSTION]
                    else:
                        continue
                            
                    if chrom2 not in list(dict_of_chrom.keys()): continue
                    else: 
                        trans_case = translocation_cases(alt_field) #determine the translocation case
                        
                        if ';END=' in info: 
                            end = int(info.split(';END=')[1].split(';')[0])
                        elif '[' in alt_field:
                            end = int(alt_field.split(':')[1].split('[')[0])     #ALT[CHR2:POSTION[
                        elif ']' in alt_field:
                            end = int(alt_field.split(':')[1].split(']')[0])     #ALT]CHR2:POSTION]
                            
                        list_of_translocations.append((chrom, start, chrom2, end, trans_case))
                        
                        
    # output files
    filename_normal = "reference_at_breakpoints.fasta"

    f1 = open(filename_normal, "w")
    for d in list_of_deletions:
        define_references_for_deletions(f1, dict_of_chrom, d, l_adj)
    for ins in list_of_insertions:
        define_references_for_insertions(f1, dict_of_chrom, ins, l_adj)
    for inv in list_of_inversions:
        define_references_for_inversions(f1, dict_of_chrom, inv, l_adj)
    for trans in list_of_translocations:
        define_references_for_translocation(f1, dict_of_chrom, trans, l_adj)
    f1.close()


def define_references_for_deletions(out1, genome, deletion, side_length):
    """ Define the reference duplet """
    local_seq_size = 2 * side_length #size of the generated allelic sequence
    ch, s, e, l = deletion
    
    if abs(l) <= local_seq_size:
        ### breakpoint withOUT deletion
        header = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - (side_length + 1)): e + side_length]
        out1.write(header + seq + "\n")

    else:
        # left breakpoint withOUT deletion
        header = ">refLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - (side_length + 1)) : s + side_length]
        out1.write(header + seq + "\n")

        # right breakpoint withOUT deletion
        header = ">refRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, e - (side_length + 1)) : e + side_length]
        out1.write(header + seq + "\n")

    ### breakpoint WITH the deletion
    header_seq = ">del_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
    seq_del = genome[ch][max(0, s - (side_length + 1)) : s]
    seq_del += genome[ch][e : e + side_length]
    out1.write(header_seq + seq_del + "\n")

def define_references_for_insertions(out1, genome, insertion, side_length):
    """ Define the reference duplet """
    local_seq_size = 2 * side_length #size of the generated allelic sequence
    #side_length = int(local_seq_size / 2)
    ch, s, l, sequence = insertion
    
    #Ref
    header_seq = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
    seq_ins = genome[ch][max(0, s - side_length) : s]
    seq_ins += genome[ch][s : s + side_length]
    out1.write(header_seq + seq_ins + "\n")
    
       
    if abs(l) <= local_seq_size:
        ### breakpoint with insertion
        header = ">ins_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - side_length) : s]
        seq += sequence
        seq += genome[ch][s:s + side_length]
        out1.write(header + seq + "\n")

    else:
        # left breakpoint with insertion
        header = ">insLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - side_length) : s]
        seq += sequence[:side_length]
        out1.write(header + seq + "\n")

        # right breakpoint with insertion
        header = ">insRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = sequence[-side_length:]
        seq += genome[ch][s: s + side_length]
        out1.write(header + seq + "\n")


def define_references_for_inversions(out1, genome, inversion, side_length):
    local_seq_size = 2 * side_length #size of the generated allelic sequence
    #side_length = int(local_seq_size / 2)
    ch, s, e, l = inversion
    
    if abs(l) <= local_seq_size:
        #ref
        header = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - side_length) : e + side_length]
        out1.write(header + seq + "\n")
        
        #alt
        header = ">inv_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
        seq = genome[ch][max(0, s - side_length) : s] #add inv left side
        inversion = Seq(genome[ch][s : e])
        seq += str(inversion.reverse_complement()) #add rev comp seq
        seq += genome[ch][e : e + side_length] #add right side
        out1.write(header + seq + "\n")
        
    else:
        #ref
        header = ">refLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, s - side_length) : s + side_length]
        out1.write(header + seq + "\n")
        
        header = ">refRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
        seq = genome[ch][max(0, e - side_length) : e + side_length]
        out1.write(header + seq + "\n")
   
        #alt
        header = ">invLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
        seq = genome[ch][max(0, s - side_length) : s] #add inv left side
        inversion = Seq(genome[ch][max(0, e - side_length): e])
        seq += str(inversion.reverse_complement()) #add rev comp seq
        out1.write(header + seq + "\n")
       
        header = ">invRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
        inversion = Seq(genome[ch][s : s + side_length])
        seq = str(inversion.reverse_complement()) #add rev comp seq        
        seq += genome[ch][e : e + side_length] #add right side
        out1.write(header + seq + "\n")


def define_references_for_translocation(out1, genome, translocation, side_length):
    ''' '''
    #local_seq_size = 2 * side_length  #size of the generated allelic sequence
    #side_length = int(local_seq_size / 2)
    ch, s, chr2, e, case = translocation

    #ref
    header = ">ref_" + str(ch) + "_" + str(s) + "-" + chr2 + "-" + str(e) + "\n"
    seq = genome[ch][max(0, s - side_length) : s + side_length]
    out1.write(header + seq + "\n")
    
    #alt
    header = ">bnd_" + str(ch) + "_" + str(s) + "-" + chr2 + "-" + str(e) + "\n" 
    seq = ''
    if case == 'after_piece':  # t[p[       piece extending to the right of p is joined after t  (from VCFv4.2.pdf)
        seq = genome[ch][max(0, s - side_length) : s] #BND chr1 left side
        seq += genome[chr2][e : e + side_length] #add BND chr2 right side
    
    elif case == 'after_revcomp': # t]p]    reverse comp piece extending left of p is joined after t
        seq = genome[ch][max(0, s - side_length) : s] #BND chr1 left side
        seq += str(Seq(genome[chr2][e - side_length : e]).reverse_complement()) #add BND chr2 left side (revcomp)
        
    elif case == 'before_piece': # ]p]t     piece extending to the left of p is joined before t
        seq = genome[chr2][max(0, e - side_length) : e] #BND chr2 left side
        seq += genome[ch][s : s + side_length] #add BND chr1 right side
        
    elif case == 'before_revcomp': # [p[t   reverse comp piece extending right of p is joined before t
        seq = str(Seq(genome[chr2][e : e + side_length]).reverse_complement()) #BND chr2 right side (revcomp)
        seq += genome[ch][s : s + side_length]  #add BND chr1 right side
        
    out1.write(header + seq + "\n")



def translocation_cases(alt):
    ''' Determine the translocation cases '''
    if alt.endswith('['):
        #t[p[
        #piece extending to the right of p is joined after t
        case = 'after_piece'
        
    elif alt.endswith(']'):
        #t]p]   
        #reverse comp piece extending left of p is joined after t   
        case = 'after_revcomp'
            
    elif alt.startswith(']'):
        #]p]t
        #piece extending to the left of p is joined before t 
        case = 'before_piece'
        
    elif alt.startswith('['):
        #[p[t[
        #reverse comp piece extending right of p is joined before t 
        case = 'before_revcomp'
        
    return case


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
