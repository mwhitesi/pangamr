#!/usr/bin/env python

"""Utility functions

Commonly-used functions for loading and manipulating inputs

"""

import pandas as pd
import scipy as sp
import numpy as np
import re

def read_panseq(file):
    """Load panseq pangenome file into scipy sparse matrix

    Args:
        file (str): Panseq pan_genome.txt file
     
    Returns:
        scipy sparse matrix with genome (sample) rows and locus (feature) columns


    """

    df = pd.read_table(file, dtype={'Genome': str, 'LocusID': str})
    df = df.dropna(subset=['Contig'])
    df = df[['LocusID','Genome']]
    df = df.sort_values(by='Genome')

    df = df.apply(lambda s: s.astype("category"))

    pg = sp.sparse.coo_matrix((np.ones(df.shape[0]),(df.Genome.cat.codes, df.LocusID.cat.codes)),dtype=np.int8)
    pg = pg.tocsr()

    return (pg,df.Genome.cat.categories,df.LocusID.cat.categories)


def read_amr(file, genomes):
    """Load amr resistant/susceptible file into sparse matrix

    Args:
        file (str): PATRIC amr csv file
     
    Returns:
        sparse matrix with genome (sample) rows and antibiotic (test) columns


    """

    df = pd.read_table(file, dtype={'genome_id': str})
    df = df[['genome_id','antibiotic','resistant_phenotype']]

    df['genome_id'] = df['genome_id'].apply(lambda x: re.sub(r'(\d+)\.(\d+)',r'\1_dot_\2',str(x)))

    i=0
    df['antibiotic'] = df['antibiotic'].astype("category")
    df['resistant_phenotype'] = df['resistant_phenotype'].apply(lambda x: 1 if x == 'Resistant' else 0)
    values = np.empty(shape=(len(genomes),len(df['antibiotic'].unique())))
    values[:] = np.nan 
    for g in genomes:
        subset = df[df['genome_id'] == g]
        values[i, subset.antibiotic.cat.codes] = subset['resistant_phenotype']
        i=i+1

    return (values, subset.antibiotic.cat.categories)


def read_annot(file):
    """Load BLAST hit annotations for pangenome regions

    Args:
        file (str): 3 column tab-delimited file: locus ID, match accession, match description
     
    Returns:
        dataframe matrix


    """

    df = pd.read_table(file, dtype={'LocusID': str})

    return df

    



