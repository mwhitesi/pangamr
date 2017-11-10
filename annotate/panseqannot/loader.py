#!/usr/bin/env python

"""Load annotation files into DB

"""

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


import argparse
import logging
import os
import pandas
import re

from annot import PanseqAnnot

logger = None


def load_rgi(options):
    """Load RGI annotations into DB

    """
    df = pandas.read_table(options.input)
    pa = PanseqAnnot(options.db)

    for idx, row in df.iterrows():
        contig = row['ORF_ID']
        start = int(row['START'])
        stop = int(row['STOP'])
        annotation_string = row.to_json()

        # Why would you modify the fasta ID
        contig = re.sub(r'_\d+$', '', contig)

        print(pa.get_loci_by_location(contig, start, stop))



if __name__ == "__main__":
    """Run various programs

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('annotate.annot')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('program', choices=['rgi','resfams', 'resfinder'],
        help='Annotation program')
    parser.add_argument('db', help='Panseq directory')
    parser.add_argument('input', help='Annotation file')
    
    
    options = parser.parse_args()

    if options.program == 'rgi':
        load_rgi(options)



