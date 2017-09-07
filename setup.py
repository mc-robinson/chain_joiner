from setuptools import setup, find_packages

setup(
  name = 'chain_joiner',
  packages = ['chain_joiner'], # this must be the same as the name above
  version = '0.1',
  description = 'A package to automatically fix chain breaks in PDB files using MODELLER',
  author = 'Matthew C. Robinson',
  author_email = 'matthew.robinson@yale.edu',
  url = 'https://github.com/mc-robinson/chain_joiner', # use the URL to the github repo,
  keywords = [], # arbitrary keywords
  classifiers = [
	'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
	],
  install_requires=['biopandas','modeller'],
  entry_points={
        'console_scripts': [
            'chain_joiner=chain_joiner.chain_joiner:main',
	    'make_seq=chain_joiner.make_seq:main',
	    'make_alignment=chain_joiner.make_alignment:main',
            'make_model=chain_joiner.make_model:main',
        ],
    },

)
