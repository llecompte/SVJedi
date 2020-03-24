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
	args = parser.parse_args()

	genome_file = args.ref[0]
	vcf_file = args.vcf[0]

	create_ref(genome_file, vcf_file)


def create_ref(genome, set_of_sv):
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
					chrom, start, _, __, type_sv, ___, ____, info, *_ = line.rstrip(
						"\n"
					).split("\t")
				else:
					chrom, start, _, __, type_sv, ___, ____, info = line.rstrip(
						"\n"
					).split("\t")
				
				if 'SVTYPE' in info:
					svtype = info.split('SVTYPE=')[1].split(';')[0]

				else:
					svtype = '' 
				
				#for deletions
				if type_sv == "<DEL>" or svtype == 'DEL':
				
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
					length = len(type_sv)
								
					# focus on >=50bp length insertion
					if length >= 50:
						list_of_insertions.append((chrom, start, length, type_sv))

				#for inversions
				elif type_sv == '<INV>' or svtype == 'INV':
					start = int(start)
					end = int(info.split(';END=')[1].split(';')[0])
					length = end - start
					 
					# focus on >=50bp length inversions
					if length >= 50:
						list_of_inversions.append((chrom, start, end, length))

				#for translocations
				elif svtype == 'BND':
					start = int(start)
					chrom2 = info.split('CHR2=')[1].split(';')[0]
					if chrom2 not in list(dict_of_chrom.keys()): continue
					else: 
						end = int(info.split(';END=')[1].split(';')[0])
						list_of_translocations.append((chrom, start, chrom2, end))
						
						
	# output files
	filename_normal = "reference_at_breakpoints.fasta"

	f1 = open(filename_normal, "w")
	for d in list_of_deletions:
		define_references_for_deletions(f1, dict_of_chrom, d)
	for ins in list_of_insertions:
		define_references_for_insertions(f1, dict_of_chrom, ins)
	for inv in list_of_inversions:
		define_references_for_inversions(f1, dict_of_chrom, inv)
	for trans in list_of_translocations:
		define_references_for_translocation(f1, dict_of_chrom, trans)
	f1.close()


def define_references_for_deletions(out1, genome, deletion):
	""" Define the reference duplet """
	local_seq_size = 10000 #size of the generated allelic sequence
	side_length = int(local_seq_size / 2)
	ch, s, e, l = deletion
	
	if abs(l) <= local_seq_size:
		### breakpoint withOUT deletion
		header = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : e + side_length]
		out1.write(header + seq + "\n")

	else:
		# left breakpoint withOUT deletion
		header = ">refLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : s + side_length]
		out1.write(header + seq + "\n")

		# right breakpoint withOUT deletion
		header = ">refRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][e - side_length : e + side_length]
		out1.write(header + seq + "\n")

	### breakpoint WITH the deletion
	header_seq = ">del_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
	seq_del = genome[ch][s - side_length : s]
	seq_del += genome[ch][e : e + side_length]
	out1.write(header_seq + seq_del + "\n")


def define_references_for_insertions(out1, genome, insertion):
	""" Define the reference duplet """
	local_seq_size = 10000 #size of the generated allelic sequence
	side_length = int(local_seq_size / 2)
	ch, s, l, sequence = insertion
	
	#Ref
	header_seq = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
	seq_ins = genome[ch][s - side_length : s]
	seq_ins += genome[ch][s : s + side_length]
	out1.write(header_seq + seq_ins + "\n")
	
	   
	if abs(l) <= local_seq_size:
		### breakpoint with insertion
		header = ">ins_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : s]
		seq += sequence
		seq += genome[ch][s:s + side_length]
		out1.write(header + seq + "\n")

	else:
		# left breakpoint with insertion
		header = ">insLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : s]
		seq += sequence[:side_length]
		out1.write(header + seq + "\n")

		# right breakpoint with insertion
		header = ">insRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = sequence[-side_length:]
		seq += genome[ch][s: s + side_length]
		out1.write(header + seq + "\n")


def define_references_for_inversions(out1, genome, inversion):
	local_seq_size = 10000 #size of the generated allelic sequence
	side_length = int(local_seq_size / 2)
	ch, s, e, l = inversion
	
	if abs(l) <= local_seq_size:
		#ref
		header = ">ref_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : e + side_length]
		out1.write(header + seq + "\n")
		
		#alt
		header = ">inv_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
		seq = genome[ch][s - side_length : s] #add inv left side
		inversion = Seq(genome[ch][s : e])
		seq += str(inversion.reverse_complement()) #add rev comp seq
		seq += genome[ch][e : e + side_length] #add right side
		out1.write(header + seq + "\n")
		
	else:
		#ref
		header = ">refLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][s - side_length : s + side_length]
		out1.write(header + seq + "\n")
		
		header = ">refLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n"
		seq = genome[ch][e - side_length : e + side_length]
		out1.write(header + seq + "\n")
   
		#alt
		header = ">invLeft_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
		seq = genome[ch][s - side_length : s] #add inv left side
		inversion = Seq(genome[ch][s : s + side_length])
		seq += str(inversion.reverse_complement()) #add rev comp seq
		out1.write(header + seq + "\n")
	   
		header = ">invRight_" + str(ch) + "_" + str(s) + "-" + str(abs(l)) + "\n" 
		inversion = Seq(genome[ch][e - side_length : e])
		seq += str(inversion.reverse_complement()) #add rev comp seq        
		seq += genome[ch][e : e + side_length] #add right side
		out1.write(header + seq + "\n")


def define_references_for_translocation(out1, genome, translocation):
	''' '''
	local_seq_size = 10000 #size of the generated allelic sequence
	side_length = int(local_seq_size / 2)
	ch, s, chr2, e = translocation

	#ref
	header = ">ref_" + str(ch) + "_" + str(s) + "-" + chr2 + "-" + str(e) + "\n"
	seq = genome[ch][s - side_length : s + side_length]
	out1.write(header + seq + "\n")
	
	#alt
	header = ">bnd_" + str(ch) + "_" + str(s) + "-" + chr2 + "-" + str(e) + "\n" 
	seq = genome[ch][s - side_length : s] #add inv left side
	seq += genome[chr2][e : e + side_length] #add right side
	out1.write(header + seq + "\n")
	
		
if __name__ == "__main__":
	if sys.argv == 1:
		sys.exit("Error: missing arguments")
	else:
		main(sys.argv[1:])
