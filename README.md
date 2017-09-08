# chain_joiner: a python package to automatically fix chain breaks in PDB files
### author: ###

* [Matt Robinson](https://github.com/mc-robinson) - `<matthew.robinson@yale.edu>`

### Description ###

chain_joiner 'fills in' missing residues in PDB files and creates a new, fixed PDB file. The package mainly uses MODELLER, a homology and comparitive modelling program created by the Sali lab at UCSF. In order to fill in the missing residue gaps, MODELLER is used to create a homology model with the original PDB file serving as a template and the full sequence serving as the target sequence to model. 

The process of doing adding in the missing residues is highly based on this [helpful tutorial on the Sali lab's website](https://salilab.org/modeller/wiki/Missing%20residues):

1. Use Modeller to extract the sequence from the original PDB file with missing residues.

2. Create an alignment file, alignment.ali, between the original PDB structure
with missing residues and the full fasta sequence (usually available on PDB website). 

3. Use Modeller to create the actual homology model.

NOTE: The way the alignment file is constructed is not guarenteed to work for every possible protein. It is often a good idea to check the automatically generated alignment file if the output seems suspect. Also, the results of chain_joiner will not be that great if the loop is of substantial size. For this, more detailed homology modelling is likely needed. 

### Usage ####

Running with locally stored PDB and fasta files:

`join_chains -p pdbfile.pdb -f fastafile.fasta [options]`

Running with PDB and sequence files from the Protein Data Bank's website:

`chain_joiner_online -p PDB_ID [options]`

Input Arguments:

[options]
`-a, --automodel`
    The simplest method for simple comparitive modeling. Will not give 
    great results but suggested when many chain breaks are present. [default: True]

`-f, --fixed_automodel`
    Builds an automodel and keeps the non-missing residues fixed, 
    whereas they can move in the other methods. [default: False]

`-l, --loopmodel` 
    Builds a model by refining the loop with the missing residues.
    Suggested when have one small chain break in the PDB. [default: False]
    
   
### Installation ###

**USING CONDA:**

The simplest way to install chain_joiner is through the [conda package manager](https://conda.io/docs/):
`conda install -c mc-robinson -c conda-forge -c salilab chain_joiner`
Using this method, conda will automatically install both biopandas and modeller, which chain_joiner depends on. 

NOTE: one must first obtain a MODELLER license key (free for academic users) from the Sali lab's [website](https://salilab.org/modeller/registration.html). Then set this license equal to an enviroment variable called `KEY_MODELLER`.

**USING PIP:**

One can also download chain_joiner from the Python Package Index (PyPI). 
`pip install chain_joiner`
However, this will require that both [MODELLER](https://salilab.org/modeller/9.19/release.html#anaconda) and [biopandas](https://github.com/rasbt/biopandas) are installed separetly. MODELLER, will likely involve a conda install anyway. So it is suggested you just install chain_joiner with conda, which will handle the MODELLER install automatically.

### References ###
* Sebastian Raschka. Biopandas: Working with molecular structures in pandas dataframes. The Journal of Open Source Software, 2(14), jun 2017. doi: 10.21105/joss.00279. URL http://dx.doi.org/10.21105/joss.00279.
* B. Webb, A. Sali. Comparative Protein Structure Modeling Using Modeller. Current Protocols in Bioinformatics, John Wiley & Sons, Inc., 5.6.1-5.6.32, 2014.

### Contact ###

If you run into an error/issue with chain_joiner, please let me know at matthew.robinson@yale.edu. Feel free to also contribute and send me a pull request.


