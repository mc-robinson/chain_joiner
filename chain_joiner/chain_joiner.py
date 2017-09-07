#!/usr/bin/env python
# File name: chop_chain_joiner.py
# Author: Matt Robinson
# Date created: 5/24/2017
# Date last modified: 5/24/2017
# Python Version: 3.6
"""
Description:
This script 'fills in' missing residues in PDB files and creates a fixed PDB file. 
In order to fill in the gaps, Modeller is used to create a homology model with the 
original PDB file serving as a template and the full sequence serving as the target 
sequence to model. 

The process of doing this is as follows:

1. Call make_seq.py to extract the sequence from the original PDB file with missing residues.

2. Call make_alignment.py to create an alignment file, alignment.ali, between the original PDB structure
with missing residues and the full fasta sequence (usually available on PDB website). 

3. Call make_model.py to create the actual homology model.

Please see the headers of each of these scripts for more specific information.

Usage: python chain_joiner.py -p pdbfile.pdb -f fastafile.fasta [options]

For example, to make a loop model of PDB code 1qg8, I would call it as: 
'python chain_joiner.py -p 1qg8.pdb -f 1qg8_full_seq.fasta -a'

Input Arguments:
[optional]
-a, --automodel
    The simplest method for simple comparitive modeling. Will not give 
    great results but suggested when many chain breaks are present. [default: True]

-f, --fixed_automodel
    Builds an automodel and keeps the non-missing residues fixed, 
    whereas they can move in the other methods. [default: False]

-l, --loopmodel 
    Builds a model by refining the loop with the missing residues.
    Suggested when have one small chain break in the PDB. [default: False]

Output: A set of PDB files (number depends on the chosen method)
"""

import argparse
import sys
import os

# import modules
# from chain_joiner import make_seq
# from chain_joiner import make_alignment
# from chain_joiner import make_model

from chain_joiner import make_seq
from chain_joiner import make_alignment
from chain_joiner import make_model

def main():

    parser = argparse.ArgumentParser(
        prog='make_model.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
    """
    This script 'fills in' missing residues in PDB files and creates a fixed PDB file. 
    In order to fill in the gaps, Modeller is used to create a homology model with the 
    original PDB file serving as a template and the full sequence serving as the target 
    sequence to model. 

    The process of doing this is as follows:

    1. Call make_seq.py to extract the sequence from
       the original PDB file with missing residues.

    2. Call make_alignment.py to create an alignment file, alignment.ali, 
       between the original PDB structure with missing residues and the full
       fasta sequence (usually available on PDB website). 

    3. Call make_model.py to create the actual homology model.

    Please see the headers of each of these scripts for more specific information.

    Usage: python chain_joiner.py -p pdbfile.pdb -f fastafile.fasta [options]

    For example, to make a loop model of PDB code 1qg8, I would call it as: 
    'python chain_joiner.py -p 1qg8.pdb -f 1qg8_full_seq.fasta -a'

    @author: Matt Robinson, matthew.robinson@yale.edu
             William L. Jorgensen lab, Yale University
    
    REQUIREMENTS:
    Preferably Anaconda python 3 with following modules:
        argparse
        modeller    
    """
    )

    parser.add_argument(
        "-p", "--pdb", help="path of the pdb file with .pdb file descriptor")
    parser.add_argument(
        "-f", "--fasta", help="path of the fasta file with .fasta file descriptor")
    parser.add_argument(
        "-a", "--automodel", help="the simplest method for simple comparitive modeling", action="store_true")
    parser.add_argument(
        "-fm", "--fixed_automodel", help="builds an automodel and keeps the non-missing residues fixed", action="store_true")
    parser.add_argument(
        "-l", "--loopmodel", help="builds a model by refining the loop with the missing residues", action="store_true")

    args = parser.parse_args()

    # join the chains
    join_chains(args.pdb, args.fasta, args.automodel, args.fixed_automodel, args.loopmodel)

def join_chains(pdb_file, fasta_file, a, fm, l):

    # Get the PDB id from the file
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

    # get the sequence from the PDB file:
    make_seq.get_sequence(pdb_file)

    # make the alignment file:
    make_alignment.align(pdb_file, pdb_id + ".seq", fasta_file)

    # make the model
    make_model.model(pdb_file, a, fm, l)

    # make a folder for the output
    dir_name = './' + pdb_id +'_output/'
    os.makedirs(dir_name)

    # get a list of all output files in the working directory
    output_files = [filename for filename in os.listdir('.') if filename.startswith(pdb_id)]
    # remove the folder name
    if (pdb_id + '_output') in output_files:
        output_files.remove(pdb_id + '_output')

    # mv these files to the output folder
    for file in output_files:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + '_output/')   
        except:
            pass

if __name__ == "__main__":

    main()