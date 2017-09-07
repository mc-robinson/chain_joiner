#!/usr/bin/env python
# File name: make_seq.py
# Author: Matt Robinson
# Date created: 5/23/2017
# Date last modified: 9/05/2017
# Python Version: 3.6
"""
Description:
This script extracts the sequence from a PDB file with missing residues. 
The output is a sequence file of the form [pdb_filename].seq

Usage: python make_pdb_seq.py -p pdbfile.pdb
"""
import sys
import os
import argparse
from modeller import *

def main():

    parser = argparse.ArgumentParser(
        prog='make_seq.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
    """
    This script extracts the sequence from a PDB file with missing residues.
    Output is a sequence file of the form [pdb_filename].seq

    @author: Matt Robinson, matthew.robinson@yale.edu
             William L. Jorgensen lab, Yale University

    Usage: python make_seq.py -p pdbfile.pdb
    
    Example Usage: python make_seq.py -p 1qg8.pdb
    
    REQUIREMENTS:
    Preferably Anaconda python 3 with following modules:
        argparse
        modeller    
    """
    )

    parser.add_argument(
        "-p", "--pdb", help="full path of the pdb file with .pdb file descriptor")

    args = parser.parse_args()

    # call the main function
    get_sequence(args.pdb)

def get_sequence(pdb_file):

    #Get the PDB id from the file (just strips off extension)
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

    # Get the sequence of the PDB file, and write to an alignment file
    code = pdb_file

    e = environ()

    # directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    # Read in HETATM records from template PDBs
    e.io.hetatm = True

    m = model(e, file=code)
    aln = alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=pdb_id+'.seq')

if __name__ == "__main__":

    main()

