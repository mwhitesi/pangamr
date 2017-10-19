#!/usr/bin/env python

"""Functions for working with Phenotype data


"""

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

import argparse
import logging
import pandas as pd
import re
import pickle

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


logger = None

def to_dict(tsvfile, id_col, pheno_col, name_col=None, phenotype_name=None, idformatter=None,
        phenoformatter=None):
    """Convert from tsv file to dictionary structure expected in LocusTree

    Args:
        tsvfile: 

    """

    df = pd.read_table(tsvfile, dtype={id_col: str})
    if name_col:
        df = df[[id_col,name_col,pheno_col]]
    else:
        df = df[[id_col,pheno_col]]

    if idformatter:
        df[id_col] = df[id_col].apply(idformatter)

    if phenoformatter:
        df[pheno_col] = df[pheno_col].apply(phenoformatter)

    if name_col == phenotype_name:
        raise Exception("Missing both argument: name_col and phenotype_name. Must provide column with phenotype name "+
            "(when input file contains multiple phenotypes) "+
            "or the name of phenotype (input file contains single phenotype).")

    if name_col:
        # Multiple phenotypes
        df[name_col] = df[name_col].astype("category")

        phenodict = df.groupby(name_col).apply(
            lambda df: df.groupby(id_col).apply(
                lambda df2: df2[pheno_col].values[0]
            ).to_dict()
        ).to_dict()

    else:

        tmp = df.groupby(id_col).apply(lambda df2: df2[pheno_col].values[0]).to_dict()
        phenodict = {phenotype_name: tmp}

    return phenodict
        


def normalizefn(x):
    return re.sub(r'\.',r'_dot_',str(x))


if __name__ == "__main__":
    """Convert tsv to dict and pickle

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('panseqtrees.phenotype')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('output', help='Picke file location')
    parser.add_argument('input', help='Phenotype .tsv file')
    
    options = parser.parse_args()

    phenodict = to_dict(options.input, 'genome_id', 'resistant_phenotype', name_col='antibiotic', idformatter=normalizefn)

    with open(options.output, 'wb') as fh:
        pickle.dump(phenodict, fh)