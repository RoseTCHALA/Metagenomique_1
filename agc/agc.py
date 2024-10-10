#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Rose TCHALA SARE"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Rose TCHALA SARE"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Rose TCHALA SARE"
__email__ = "r.tchalasare@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    # Determine the appropriate file opening method based on the file extension
    if str(amplicon_file).endswith(".gz"):
        open_func = gzip.open  # Use gzip.open for compressed files
    else:
        open_func = open  # Use open for regular text files
    
    # Read the file and process the sequences
    with open_func(amplicon_file, "rt") as monfich:
        sequence = ""
        for line in monfich:
            line = line.strip()
            if line.startswith(">"):  # Sequence header
                if len(sequence) >= minseqlen:  # Yield the sequence if it's long enough
                    yield sequence
                sequence = ""  # Reset sequence for the next entry
            else:
                sequence += line  # Append the sequence line

        # After the loop, yield the last sequence if it meets the minimum length
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    stock = {}
    sequences = read_fasta(amplicon_file,minseqlen)

    for sequence in sequences :
        if sequence not in stock.keys() :
            stock[sequence] = 1
        else :
            stock[sequence] += 1


    sorted_stock = sorted(stock.items(), key=lambda x:x[1], reverse=True)
    converted_dict = dict(sorted_stock)

    for key in converted_dict:
        if converted_dict[key] >= mincount :
            yield [key, stock[key]]




def get_identity(alignment_list: List[str]) -> float:
    
    seq = alignment_list[0]
    seq1 = alignment_list[1]

    equal = 0

    if len(seq) != 0 :

        for i in range(0,len(seq)) :
            if seq[i] == seq1[i] :
                equal += 1 

        return (equal/len(seq) * 100)
    else :
        return 0
    

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    
    occurences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))

    total = len(occurences)

    seq = occurences[0]
    
    #La première séquence est potentiellement une OTU
    OTU = [occurences[0]]

    for i in range(1, total):
        seq = occurences[i]

        sequence_test = seq[0]
        abondance_test = seq[1]


        #On regarde si la séquence en question match ou match pas avec une OTU déja présente
        for j in range(0, len(OTU)):
            seq_OTU = OTU[j]
            if abondance_test >= seq_OTU[1] :
                a = nw.global_align(sequence_test, seq_OTU[0], gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
                if get_identity(a) >= 97 :
                    OTU.remove(seq_OTU)
                    OTU.append(seq)
        #Si c'est pas le cas, on la teste avec toutes les autres séquences
        is_OTU = True

        for k in range(i + 1 , total):
            seq1 = occurences[k]
            sequence_test1 = seq1[0]
            abondance_test1 = seq1[1]

            if abondance_test1 > abondance_test:
                a = nw.global_align(seqsequence_testuence1, seq_OTU[0], gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
                
                if get_identity(a) >= 97 :
                    is_OTU = False

        if is_OTU == True :
            OTU.append(occurences[i])
    return OTU


        

    




def write_OTU(OTU_list: List, output_file: Path) -> None:
    with open(output_file, "w",newline="\n") as file:
        for i in range(0, len(OTU_list)):
            seq = OTU_list[i]
            file.write( f">OTU_{i+1} occurence:{seq[1]}\n")
            file.write( f"{textwrap.fill(seq[0], width = 80 )}\n")





#===========================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    sequences = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, 0, 0)
    write_OTU(sequences, args.output_file)
    



if __name__ == '__main__':
    main()
