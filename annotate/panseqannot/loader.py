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
from collections import defaultdict

from annot import PanseqAnnot

logger = None

def load_resfams(options):
    """Load Resfams annotations into DB

    """
    df = pandas.read_table(options.input, sep='\s+', comment='#', header=None, 
        names=('target', 'tacc', 'query', 'qacc', 'evalue', 'score', 'bias',
            'domain_evalue', 'domain_score', 'domain_bias',
            'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description'))
    pa = PanseqAnnot(options.db)

    proteins = {}

    for idx, row in df.iterrows():
        target = row['target']

        if target in proteins:
            if row['score'] > proteins[target]['score']:
                proteins[target] = row
        else:
            proteins[target] = row


    for row in proteins.values():
        target = row['target']

        m1 = re.search(r'^(.+)\|strand=(\-?\d);frame=(\d);relstart=(\d+);relend=(\d+)$', target)

        if m1:
            name, strand, frame, rstart, rend = m1.groups()
            start, stop = sorted([int(rstart), int(rend)])

            annotation_string = row.to_json()

            m = re.search(r'^(?P<contig>.+)_\((?P<begin>\d+)\.\.(?P<end>\d+)\)$|^(?P<contigonly>.+)$', 
                str(name))
            if m:
                contig, begin, end, contigonly = m.groups()
                logger.info("Loaded HIT: {}, {}-{}".format(name, start, stop))
                if contigonly:
                    fragments = pa.load_annotation(annotation_string, 'resfams', contigonly, start, stop)
                    logger.info("\tmapped to {} fragments.".format(len(fragments)))
                    for f in fragments:
                        logger.debug("\t\t{}".format(str(f)))
                else:
                    relpos = int(begin)
                    fragments = pa.load_annotation(annotation_string, 'resfams', contig, start+relpos, stop+relpos)
                    logger.info("\tmapped to {} fragments.".format(len(fragments)))
                    for f in fragments:
                       logger.debug("\t\t{}".format(str(f)))
            else:
                raise Exception("Invalid contig name: {}".format(name))
                    
        else:
            raise Exception("Invalid target name: {}".format(target))     


def load_resfinder(options):
    """Load Resfinder annotations into DB

    """
    df = pandas.read_table(options.input, comment='#', header=None, 
        names=('qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 
            'send', 'slen', 'evalue', 'bitscore'))
    pa = PanseqAnnot(options.db)

    loci = {}
    for idx, row in df.iterrows():
        name = row['qseqid']
        start = int(row['qstart'])
        stop = int(row['qend'])
        bitscore = float(row['bitscore'])
        subject = row['sseqid']
        coverage = float(row['length'])/float(row['slen'])

        if coverage > 0.6:
            # over 60% of gene must align with query
            
            # Keep best hit for each position
            addr = locad(name, start, stop)
            if addr in loci:
                if loci[addr]['bitscore'] < bitscore:
                    loci[addr] = row
                elif loci[addr]['bitscore'] == bitscore and loci[addr]['sseqid'] > subject:
                    loci[addr] = row
            else:
                loci[addr] = row

    # Remove overlapping hits, keeping longest
    locations = defaultdict(list)
    for row in loci.values():
        query = row['qseqid']
        locations[query].append(sorted((row['qstart'], row['qend'])))

    hits = non_overlapping(locations)

    for key in hits:
        row = loci[key]

        name = row['qseqid']
        start = int(row['qstart'])
        stop = int(row['qend'])
        annotation_string = row.to_json()

        m = re.search(r'^(?P<contig>.+)_\((?P<begin>\d+)\.\.(?P<end>\d+)\)$|^(?P<contigonly>.+)$', 
            str(name))
        if m:
            contig, begin, end, contigonly = m.groups()
            logger.info("Loaded HIT: {}, {}-{}".format(name, start, stop))
            if contigonly:
                fragments = pa.load_annotation(annotation_string, 'resfinder', contigonly, start, stop)
                logger.info("\tmapped to {} fragments.".format(len(fragments)))
                for f in fragments:
                    logger.debug("\t\t{}".format(str(f)))
            else:
                relpos = int(begin)
                fragments = pa.load_annotation(annotation_string, 'resfinder', contig, start+relpos, stop+relpos)
                logger.info("\tmapped to {} fragments.".format(len(fragments)))
                for f in fragments:
                   logger.debug("\t\t{}".format(str(f)))
                    
        else:
            raise Exception("Invalid contig name: {}".format(name))


def locad(contig, start, stop):
    # Loci address
    s = sorted([start, stop])
    return '{}:{}-{}'.format(contig, *s)


def non_overlapping(locations):
    # Remove hits that are contained within other hits

    def position_sort(x, y):
        if x[0] == y[0]:
            return y[1] - x[1]
        else:
            return x[0] - y[0]

    final = []

    for q in locations.keys():
    
        starts = sorted(locations[q], cmp=position_sort)

        # Save the first one - can't be overlapping
        s = 0
        addr = locad(q, starts[s][0], starts[s][1])
        final.append(addr)

        mx = len(starts)
        while s < (mx-1):
            # Find the next non-overlapping in the list 
            l = [i for i in range(s+1,mx) if starts[i][0] >= starts[s][1]]
            if l:
                s = l[0]
                addr = locad(q, starts[s][0], starts[s][1])
                final.append(addr)
            else:
                # Exit loop, no more non-overlapping
                s = mx
            
    return final


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
    elif options.program == 'resfinder':
        load_resfinder(options)
    elif options.program == 'resfams':
        load_resfams(options)



