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

def main():

    #create parser object
    parser = argparse.ArgumentParser(
        prog='run_chain_joiner_online.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    SCRIPT TO FIX MISSING RESIDUES IN PDB FILES.
    This version obtains the PDB and FASTA files 
    for the user directly from the web.

    @author: Matthew Robinson matthew.robinson@yale.edu
    @author: William L. Jorgensen Lab 
    
    REQUIREMENTS:
    Preferably Anaconda python 3 environment with following modules:

    Modeller
    pandas 
    biopandas
    argparse
    getopt
    numpy
    """)

    #defining arguments for parser object
    parser.add_argument('pdb_id')
    parser.add_argument(
        "-a", "--automodel", help="The simplest method for simple comparitive modeling.")
    parser.add_argument(
        "-l", "--loopmodel", help="Builds a model by refining the loop with the missing residues.")
    parser.add_argument(
        "-f", "--fixed_automodel", help="Builds an automodel and keeps the non-missing residues fixed.")

    #parse the arguments from standard input
    args = parser.parse_args()

    #pdb_id = sys.argv[1]
    pdb_id = str(args.pdb_id)

    #make url's for wget
    pdb_url = 'https://files.rcsb.org/download/' + pdb_id + '.pdb'
    #fasta_url = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=' + pdb_id + '&compressionType=uncompressed'
    fasta_url = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=' + pdb_id

    #make file names
    pdb_fn = pdb_id + '.pdb'
    fasta_fn = pdb_id + '.fasta'

    #get the files from the pdb_database
    os.system('wget --no-check-certificate -O ' + pdb_fn + ' ' + pdb_url)
    os.system('wget -O ' + fasta_fn +'.gz ' + fasta_url + ' && gzip -d '+ fasta_fn +'.gz')
    #os.system('wget -q -O - ' + fasta_url + ' | gzip -d > ' + fasta_fn)
    #os.system('gzip -d ' + fasta_fn + '.gz')

    #call the functions based on flag
    if args.automodel:
        subprocess.call(['python', 'chain_joiner.py', pdb_id + '.pdb', pdb_id + '.fasta','-a'])

    if args.loopmodel:
        subprocess.call(['python', 'chain_joiner.py', pdb_id + '.pdb', pdb_id + '.fasta','-l'])

    elif args.fixed_automodel:
        subprocess.call(['python', 'chain_joiner.py', pdb_id + '.pdb', pdb_id + '.fasta','-f'])

    else:
        subprocess.call(['python', 'chain_joiner.py', pdb_id + '.pdb', pdb_id + '.fasta','-a'])

    #make a folder for the output
    dir_name = './' + pdb_id +'_output/'
    os.makedirs(dir_name)

    #get a list of all output files in the working directory
    output_files = [filename for filename in os.listdir('.') if filename.startswith(pdb_id)]
    #remove the folder name
    if (pdb_id + '_output') in output_files:
        output_files.remove(pdb_id + '_output')

    #mv these files to the output folder
    for file in output_files:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + '_output/')   
        except:
            pass

if __name__ == "__main__":
    # calling the main function
    main()
