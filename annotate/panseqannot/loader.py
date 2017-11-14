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
        name = row['ORF_ID']
        start = int(row['START'])
        stop = int(row['STOP'])
        annotation_string = row.to_json()

        m = re.search(r'^(?P<contig>.+)_\((?P<begin>\d+)\.\.(?P<end>\d+)\)_\d+$|^(?P<contigonly>.+)_\d+$', 
            str(name))
        if m:
            contig, begin, end, contigonly = m.groups()
            logger.info("Loaded HIT: {}, {}-{}".format(name, start, stop))
            if contigonly:
                fragments = pa.load_annotation(annotation_string, 'rgi', contigonly, start, stop)
                logger.info("\tmapped to {} fragments.".format(len(fragments)))
                for f in fragments:
                    logger.debug("\t\t{}".format(str(f)))
            else:
                relpos = int(begin)
                fragments = pa.load_annotation(annotation_string, 'rgi', contig, start+relpos, stop+relpos)
                logger.info("\tmapped to {} fragments.".format(len(fragments)))
                for f in fragments:
                   logger.debug("\t\t{}".format(str(f)))
                    
        else:
            raise Exception("Invalid ORF ID: {}".format(name))
        

if __name__ == "__main__":
    """Run various programs

    """

    logging.basicConfig(level=logging.DEBUG)
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



