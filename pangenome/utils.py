#!/usr/bin/env python

"""Utility functions

Commonly-used functions for loading and manipulating inputs

"""

import pandas as pd
import scipy as sp
import numpy as np
import re
import collections
from sklearn.externals import joblib

import config

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


SalmonellaSet = collections.namedtuple('SalmonellaSet', ['X_train', 'y_train', 'X_test', 'y_test', 'locus_index', 'serovar_index'])

def load_salmonella_data(serovars=None):
    """Load pangenome and amr data for some/all salmonella species

    Args:
        file (str): 3 column tab-delimited file: locus ID, match accession, match description
     
    Returns:
        SalmonellaSet namedtuple with attributes:
            X_train
            X_test
            y_train
            y_test
            locus_index
            serovar_index

    """
    
    amr_df = joblib.load(config.S['amr'])
    test_train_index = joblib.load(config.S['test_train_index'])
    serovar_index = joblib.load(config.S['serovar_index'])
    pg = joblib.load(config.S['pg'])
    locus_index = joblib.load(config.S['locus_index'])

    if serovars and len(serovars):
        filtered_rows = np.in1d(serovar_index, serovars)
        final_amr_df = amr_df.iloc[filtered_rows,:]
        final_test_train_index = test_train_index[filtered_rows]
        final_serovar_index = serovar_index[filtered_rows]
        final_pg = pg[filtered_rows,:]
    else:
        final_amr_df = amr_df
        final_test_train_index = test_train_index
        final_serovar_index = serovar_index        
        final_pg = pg

    training_rows = np.array(final_test_train_index == 'Training')
    X_train = final_pg[training_rows,:]
    X_test = final_pg[~training_rows,:]
    y_train = final_amr_df.iloc[training_rows,:]
    y_test = final_amr_df.iloc[~training_rows,:]

    return SalmonellaSet(X_train,y_train,X_test,y_test,locus_index,final_serovar_index)





