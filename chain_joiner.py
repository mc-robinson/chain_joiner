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

1. Call make_hetatm_seq.py to extract the sequence from the original PDB file with missing residues.

2. Call make_alignment.py to create an alignment file, alignment.ali, between the original PDB structure
with missing residues and the full sequence. 

3. Call make_model.py to create the actual homology model.

Please see the headers of each of these scripts for more specific information.

Usage: python chain_joiner.py pdb_filename.pdb fasta_filename.fasta [options]

For example, to make a loop model of PDB code 1qg8, I would call it as: 'python chain_joiner.py 1qg8.pdb 1qg8_full_seq.fasta -a'
with both 1qg8.pdb and 1qg8_full_seq.fasta in the local directory. 

Note: the PDB file and fasta file of the full sequence must be in this same directory as this code is run in. 

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

import subprocess
import sys, getopt
import os

#Get the PDB id from the file
pdb_id = os.path.splitext(os.path.basename(sys.argv[1]))[0]

def main(argv):

    try:                                
        opts, args = getopt.getopt(argv, "afl", ["automodel", "fixed_automodel", "loopmodel"]) 

    except getopt.GetoptError:           
        usage()                          
        sys.exit(2) 

    # make -a the default
    if (opts == []):
        opts = [('-a', '')]

    for opt, arg in opts:

        if opt in ('-f','--fixed_automodel'):

            #get the sequence from the PDB file:
            subprocess.call(['python', 'make_hetatm_seq.py', str(sys.argv[1])])

            #make the alignment file:
            subprocess.call(['python', 'make_alignment.py', str(sys.argv[1]), pdb_id + ".seq", str(sys.argv[2])])

            #make the model
            subprocess.call(['python', 'make_model.py', str(sys.argv[1]),'-f'])

        elif opt in ('-l','--loopmodel'):
            
            #get the sequence from the PDB file:
            subprocess.call(['python', 'make_hetatm_seq.py', str(sys.argv[1])])

            #make the alignment file:
            subprocess.call(['python', 'make_alignment.py', str(sys.argv[1]), pdb_id + ".seq", str(sys.argv[2])])

            #make the model
            subprocess.call(['python', 'make_model.py', str(sys.argv[1]),'-l'])

        else:
            
            #get the sequence from the PDB file:
            subprocess.call(['python', 'make_hetatm_seq.py', str(sys.argv[1])])

            #make the alignment file:
            subprocess.call(['python', 'make_alignment.py', str(sys.argv[1]), pdb_id + ".seq", str(sys.argv[2])])

            #make the model
            subprocess.call(['python', 'make_model.py', str(sys.argv[1]),'-a'])

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

def usage():
    print("Please see header of python script for details of proper usage")


if __name__ == "__main__":
    main(sys.argv[3:])