# !/usr/bin/env python
# File name: run_ligpargen_regression_tests.py
# Author: Matt Robinson, Jorgensen Lab @ Yale
# Email: matthew.robinson@yale.edu
# Python Version: 3.6

# to get the module imports to work, need to add .. to python path
import sys, os
testdir = os.path.dirname(__file__)
srcdir = '..'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

# import the package
import chain_joiner
# import the modules
from chain_joiner import chain_joiner
from chain_joiner import run_chain_joiner_online
# now import all functions
#from LigParGen.Conve import *

molecules_list = ['1qg8','5uiq']

def run_tests():

    FILE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

    # bool that will go to false if one of the tests failed
    passed_all_tests = True

    for mol in molecules_list:

        pdb_path = os.path.join(FILE_DIR, 'test_files', mol + '.pdb')
        fasta_path = os.path.join(FILE_DIR, 'test_files', mol + '.fasta')
        output_path = './' + mol +'_output/'

        if mol == '1qg8':

            test_path = os.path.join(FILE_DIR, 'test_files', mol + '_fixed_automodel.pdb') 
            fixed_pdb_path = './' + mol + '_output/' + mol + '_fill.B99990001.pdb'

            # test fixed_automodel
            try:
                chain_joiner.main(pdb_path, fasta_path, a=False, fm=True, l=False)
                assert(os.path.getsize(test_path)==168437)
                #assert(open(fixed_pdb_path,'r').read() == open(test_path,'r').read())
            except:
                print("FIXED_AUTOMODEL FAILED ON " + mol)
                passed_all_tests = False
                sys.exit()

            # remove folder (be carful with this)
            os.system('rm -rf ' + './' + mol + '_output')

            try:
                chain_joiner.main(pdb_path, fasta_path, a=False, fm=False, l=True)
            except:
                print("LOOPMODEL FAILED ON " + mol)
                passed_all_tests = False

            # remove folder (be carful with this)
            os.system('rm -rf ' + './' + mol + '_output')

        elif mol == '5uiq':

            # test_path = os.path.join(FILE_DIR, 'data', mol + '_automodel.pdb') 

            # test automodel
            try:
                chain_joiner.main(pdb_path, fasta_path, a=True, fm=False, l=False)
                # assert(os.path.getsize(pdb_path) == os.path.getsize(test_path))
                # assert(open(pdb_path,'r').read() == open(test_path,'r').read())
            except:
                print("AUTOMODEL FAILED ON " + mol)
                passed_all_tests = False

            # remove folder (be carful with this)
            os.system('rm -rf ' + './' + mol + '_output')

    if passed_all_tests:
        print('PASSED ALL TESTS')

if __name__ == '__main__':
    run_tests()