#!/usr/bin/env python
# File name: run_chain_joiner.py
# Author: Matt Robinson
# Date created: 06/09/2017
# Date last modified: 06/09/2017
# Python Version: 3.6
"""
This program uses wget to obtain the files needed for chain_joiner from the PDB website.
It then runs chain_joiner in automodel mode. 

usage: python run_chain_joiner.py pdb_code
"""

import os
import sys
import subprocess
import argparse

from chain_joiner import chain_joiner

def main():
    #create parser object
    parser = argparse.ArgumentParser(
        prog='run_chain_joiner_online.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
    """
    SCRIPT TO FIX MISSING RESIDUES IN PDB FILES.
    This version obtains the PDB and FASTA files 
    for the user directly from the PDB website.

    @author: Matthew Robinson, matthew.robinson@yale.edu
    @author: William L. Jorgensen Lab, Yale University 

    Usage: python run_chain_joiner_online.py -p 1qg8 -a
    
    REQUIREMENTS:
    Preferably Anaconda python 3 environment with following modules:    
        Modeller
        pandas 
        biopandas
        argparse
    """
    )

    # defining arguments for parser object
    parser.add_argument(
        "-p", "--pdb_id", help="pdb four letter code for the protein of interest")
    parser.add_argument(
        "-a", "--automodel", help="the simplest method for simple comparitive modeling", action="store_true")
    parser.add_argument(
        "-fm", "--fixed_automodel", help="builds an automodel and keeps the non-missing residues fixed", action="store_true")
    parser.add_argument(
        "-l", "--loopmodel", help="builds a model by refining the loop with the missing residues", action="store_true")

    #parse the arguments from standard input
    args = parser.parse_args()

    # call the model function
    model(args.pdb_id, args.automodel, args.fixed_automodel, args.loopmodel)

def model(pdb_id, a, fm, l):

    # make url's for wget
    pdb_url = 'https://files.rcsb.org/download/' + pdb_id + '.pdb'
    # fasta_url = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=' + pdb_id + '&compressionType=uncompressed'
    fasta_url = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=' + pdb_id

    #make file names
    pdb_fn = pdb_id + '.pdb'
    fasta_fn = pdb_id + '.fasta'

    #get the files from the pdb_database
    os.system('wget --no-check-certificate -O ' + pdb_fn + ' ' + pdb_url)
    os.system('wget -O ' + fasta_fn +'.gz ' + fasta_url + ' && gzip -d '+ fasta_fn +'.gz')
    #os.system('wget -q -O - ' + fasta_url + ' | gzip -d > ' + fasta_fn)
    #os.system('gzip -d ' + fasta_fn + '.gz')

    # run chain joiner
    chain_joiner.join_chains(pdb_fn, fasta_fn, a, fm, l )

if __name__ == "__main__":

    main()