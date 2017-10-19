#!/usr/bin/env python

"""Build tree from input sequences and generate tree image

Script for building phylogenetic tree from alignment

"""


#!/usr/bin/env python

"""Build tree from alignment

Script for building phylogenetic tree from alignment

"""

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

import argparse
import logging

from Bio import SeqIO

from panseqtrees import pangenome, phenotype

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"

logger = None

if __name__ == "__main__":
    """Build tree image

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('panseqtrees.build_locus_tree')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('output', help='.png image file location')
    parser.add_argument('input', help='Fasta file')
    parser.add_argument('--pheno', help='Phenotype .tsv file')
    parser.add_argument('--name', help='Locus name')
   
    options = parser.parse_args()

    # Read sequences from file
    fasta = SeqIO.parse(options.input, 'fasta')
    seqs = dict()
    for record in fasta:
    	seqs[record.id] = record.seq

    # Phenotype provided
    if options.pheno:
    	phenodict = phenotype.to_dict(options.pheno, 'genome_id', 'resistant_phenotype', 
    		name_col='antibiotic', idformatter=phenotype.normalizefn)

    	phenodict = {'ampicillin': phenodict['ampicillin']}

    else:
    	phenodict = None

    l = pangenome.build_tree(seqs, phenodict)
    print(l._tree)
    l.render(options.output)

   
   
        

   
  

