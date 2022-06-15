#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/retrieve_pdb.py

Usage: python -m src.retrieve_pdb
___
--help | -h Display documentation

Description

This file can also be imported as a module and contains the following functions:
    * retrieve_pdbs: retrieves pdb from Uniprot

TO DO:
    *
"""

from __future__ import absolute_import, division, print_function
import argparse

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"

import os, sys
import warnings
from urllib.error import HTTPError

import pandas as pd
from pypdb import get_pdb_file

warnings.filterwarnings('ignore')

def retrieve_pdbs(df, outdir):
    """https://www.uniprot.org/uploadlists/

    Arguments:
        df {str} -- Filepath from Uniprot Retrieve/ID mapping tool to PDB
        outdir {str} -- Directory to save pdbs
    """
    not_downloaded = []
    for uniprot_id, pdb_id in zip(df['From'], df['To']):
        try:
            pdb_file = get_pdb_file(pdb_id, filetype='pdb', compression=False)
            with open(os.path.join(outdir, f'{pdb_id}-{uniprot_id}.pdb'), 'w') as pdbf:
                pdbf.write(pdb_file)
    #         print(f'{pdb_id}-{uniprot_id}.pdb downloaded')
        except (HTTPError):
            not_downloaded.append(pdb_id)
    #         print(f'{pdb_id}-{uniprot_id}.pdb not downloaded')
            pass
    with open(os.path.join(outdir, f'not_downloaded_pdbs.log'), 'w') as ndwf:
            ndwf.write('\n'.join(not_downloaded))

def main():
    pdbs_dir = os.path.join(os.path.dirname(__file__), '../data/external/pdbs')
    df_map = pd.read_csv(os.path.join(os.path.dirname(__file__), '../data/external/uniprot-to-pdb.tsv'), sep='\t')
    retrieve_pdbs(df_map.iloc[10308:-1], pdbs_dir)

if __name__ == '__main__':
    main()
