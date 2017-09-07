#!/usr/bin/env python
# File name: make_model.py
# Author: Matt Robinson
# Date created: 5/24/2017
# Date last modified: 9/05/2017
# Python Version: 3.6
"""
Description:
This script uses the alignment file, alignment.ali, to make a homology model that includes the previously missing residues. 
In this case, the original PDB structure is the template while the full sequence is the target of the model. 

Usage: python make_model.py pdbfile.pdb [options]

Input Arguments:
[optional]
-a, --automodel
    The simplest method for simple comparitive modeling. 
    Will not give great results but suggested when many chain breaks are present. [default: True]

-f, --fixed_automodel
    Builds an automodel and keeps the non-missing residues fixed,
    whereas they can move in the other methods. [default: False]
    
-l, --loopmodel 
    Builds a model by refining the loop with the missing residues. 
    Suggested when have one small chain break in the PDB. [default: False]

Output: A set of PDB files (number depends on the chosen method)

Note: The alignment file must be in the same directory as this script.
"""
import sys
import argparse
import os
from modeller import *
from modeller.automodel import *    # Load the automodel class

def main():

    parser = argparse.ArgumentParser(
        prog='make_model.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
    """
    This script uses the alignment file, alignment.ali, 
    to make a homology model that includes the previously missing residues. 
    In this case, the original PDB structure is the template
    while the full sequence is the target of the model. 

    @author: Matt Robinson, matthew.robinson@yale.edu
             William L. Jorgensen lab, Yale University

    For simple automodel,
    Usage: python make_model.py -p pdbfile.pdb -a

    For fixed automodel,
    Usage: python make_model.py -p pdbfile.pdb -f

    For loopmodel,
    Usage: python make_model.py -p pdbfile.pdb -l
    
    REQUIREMENTS:
    Preferably Anaconda python 3 with following modules:
        argparse
        modeller    
    """
    )

    parser.add_argument(
        "-p", "--pdb", help="full path of the pdb file with .pdb file descriptor")
    parser.add_argument(
        "-a", "--automodel", help="the simplest method for simple comparitive modeling", action="store_true")
    parser.add_argument(
        "-f", "--fixed_automodel", help="builds an automodel and keeps the non-missing residues fixed", action="store_true")
    parser.add_argument(
        "-l", "--loopmodel", help="builds a model by refining the loop with the missing residues", action="store_true")

    args = parser.parse_args()

    # call the model function
    model(args.pdb, args.automodel, args.fixed_automodel, args.loopmodel)

def model(pdb_file, a, f, l):

    log.verbose()
    env = environ()

    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']
    # Read in HETATM records from template PDBs
    env.io.hetatm = True


    #first need to get PDB data so can find missing residues
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

    with open(pdb_file) as pdb:
        pdb_file_data = pdb.readlines()

    with open('./' + pdb_id + '_alignment.ali') as aln_file:
        aln_data = aln_file.readlines()

    # do modeller stuff
    pdb_seq = get_pdb_seq(aln_data)

    if (f):
        #make the string for selecting the residues to move (the missing residues)
        selection_str = make_sel_str(find_missing_residues(pdb_seq))
        #selection_str = "self.residue_range('647:B', '648:B')"
        print(selection_str)

        #build the model
        class MyModel(automodel):
            def select_atoms(self):
                #select only the missing residues
                return eval('selection(' + selection_str + ')') #need to use eval b/c sel_str is str

        a = MyModel(env, alnfile = (pdb_id + '_alignment.ali'),
                    knowns = pdb_file, sequence = pdb_id + '_fill',
                    assess_methods=(assess.DOPE, assess.GA341))

        #build only one model
        a.starting_model= 1
        a.ending_model  = 1

        a.make()

    elif (l):
        a = loopmodel(env, alnfile = (pdb_id + '_alignment.ali'),
                      knowns = pdb_file, sequence = pdb_id + '_fill',
                      assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model   = 1

        a.loop.starting_model = 1
        a.loop.ending_model   = 2
        a.loop.md_level       = refine.fast

        a.make()

    else:
        a = automodel(env, alnfile = (pdb_id + '_alignment.ali'),
                      knowns = pdb_file, sequence = pdb_id + '_fill',
                      assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model= 1
        a.ending_model  = 1

        a.make()

def find_missing_residues(pdb_seq):
    ##first delete all / from string
    #pdb_seq = pdb_seq.replace("/","")
    missing_res_lists = []
    res_number = 0 

    #create a list for holding chain labels
    chain_labels = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']


    #first split the chains
    chains_l = pdb_seq.split('/')
    # go through and find every dash in string and note index
    chain_number = 0
    res_idx = 0
    for pdb_seq in chains_l:

        missing_res_l = []
        for i in range(len(pdb_seq)):
            if (pdb_seq[i] == '-'):
                missing_res_l.append(res_idx+1) #+1 b/c Modeller starts at 1
            res_idx = res_idx + 1

        #Split missing_res into seperate loops
        chain_missing_res_lists = [] 

        start_idx = 0
        for i in range(len(missing_res_l)-1): 
            #split if not consecutive residues
            if (not(missing_res_l[i] + 1 == missing_res_l[i+1])):
                chain_missing_res_lists.append(missing_res_l[start_idx:i+1])   
                start_idx = i+1   
            #also split if end of list
            if (i==len(missing_res_l)-2):
                chain_missing_res_lists.append(missing_res_l[start_idx:i+2])

        print(chain_missing_res_lists)

        #if only one chain, don't need to add residue numbers. Just make string
        if (len(chains_l)==1):
            for lst in chain_missing_res_lists:
                lst = [str(i) for i in lst]
                missing_res_lists.append(lst)

        #go over all lists and add chain identifier if there is more than one chain 
        else:
            for lst in chain_missing_res_lists:
                lst = [str(i) for i in lst]
                lst_idx = 0
                for atom_str in lst:
                    lst[lst_idx] = atom_str + ':' + chain_labels[chain_number]
                    lst_idx = lst_idx + 1
                missing_res_lists.append(lst)

        chain_number = chain_number + 1


    print(missing_res_l)
    print(missing_res_lists)
    return missing_res_lists            

def make_sel_str(missing_res_ls):
    sel_str = ''
    for l in missing_res_ls:
        # need to add offset since Modeller numbering is different
        first_idx = l[0]
        last_idx = l[-1]
        # make the str to use as an argument
        sel_str = sel_str + "self.residue_range('" + str(first_idx) + "', '" + str(last_idx) + "'),"
    # take off the final comma
    sel_str = sel_str[0:-1]
    return sel_str

def get_pdb_seq(aln_data):
    pdb_seq = ""
    seq_data = aln_data[3:]
    for line in seq_data:
        line = line.rstrip()
        if (line[0]=='>'):
            break
        pdb_seq = pdb_seq + line
    #remove the break character '/'
    #pdb_seq = re.sub('/','',pdb_seq)
    #print(pdb_seq)
    return pdb_seq


if __name__ == "__main__":

    main()

