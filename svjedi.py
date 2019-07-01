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

import os
import sys
import subprocess
import shlex
import argparse
from modules import generateRef
from modules import genotype


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Structural variations genotyping using long reads")
    parser.add_argument("-v", "--vcf", metavar="<vcffile>", help="vcf format", required=True)

    parser.add_argument("-r", "--ref", metavar="<reffile>", nargs=1, help="fasta format")

    parser.add_argument("-i", "--input", metavar="<readfile>", nargs="*", help="reads")

    parser.add_argument("-p", "--paf", metavar="<paffile>", nargs=1, help="PAF format")

    parser.add_argument("-t", "--threads", metavar="<nb_threads>", help="Number of threads")

    parser.add_argument("-o", "--output", metavar="<output>", nargs=1, help="genotype output file")
    
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

    parser.add_argument(
        "-d",
        "--data",
        metavar="<seq data type>",
        nargs="?",
        default="pb",
        choices=["ont", "pb"],
        type=str,

    )

    return parser.parse_args()


def main(args):
    """ Pipeline for genotyping SVs """
    launch_ref, launch_align = True, True
    args = parse_arguments(args)  # parse arguments

    # check if all required arguments are provided
    if not any(
        [
            all([arg is not None for arg in [args.vcf, args.ref, args.input]]),
            all([args.paf is not None, args.vcf is not None]),
        ]
    ):
        sys.exit(
            "User must provide either VCF and REF and READS files OR only PAF file"
        )

    if all([arg is not None for arg in [args.vcf, args.ref, args.input]]):
        vcf_file = args.vcf
        ref_file = args.ref[0]
        reads_file = args.input

    elif args.paf is not None:
        vcf_file = args.vcf
        paf_file = args.paf[0]
        launch_ref, launch_align = False, False  # skip the ref + minimap cmd

    # optional argument
    if args.output is None:
        output_file = "genotype_results.vcf"
    else:
        output_file = args.output[0]

    if args.threads is not None:
        threads = args.threads

    data_type = args.data
    min_support = args.minsupport
    min_support = args.minsupport
    d_over = args.dover
    d_end = args.dend

    # generate ref sequence
    if launch_ref is True:
        generateRef.create_ref(ref_file, vcf_file)

    # map with minimap2
    if launch_align is True:
        paf_file = "minimap_results.paf"
        outPaf = open(paf_file, "w")
        outErr = open("/dev/null", "w")
        if args.threads is not None:
            threads = args.threads
            cmd_minimap = (
                "minimap2 -cx "
                + "map-"
                + str(data_type)
                + " -t "
                + threads
                + " "
                + "reference_at_breakpoints.fasta"
                + " "
                + " ".join(reads_file)
            )
        else:
            cmd_minimap = (
                "minimap2 -cx "
                + "map-"
                + str(data_type)
                + " "
                + "reference_at_breakpoints.fasta"
                + " "
                + " ".join(reads_file)
            )
        
        args = shlex.split(cmd_minimap)
        subprocess.Popen(args, stdout=outPaf, stderr=outErr).communicate()
        outPaf.close()
        outErr.close()

    # compute genotype
    genotype.genotype(paf_file, vcf_file, output_file, min_support, d_over, d_end)


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
