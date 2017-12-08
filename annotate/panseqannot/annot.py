#!/usr/bin/env python

"""Annotate panseq pangenomes

Functions for running AMR gene predictors and cross-referencing them against the panseq pangenome

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
import sqlite3
from collections import defaultdict

logger = None

class PanseqAnnot(object):
    """Storage and Retrieval of Panseq Annotations

    Maps annotations on the unfragmented pangenome sequences to the individual 
    pangenome fragments.

    """

    def __init__(self, panseq_directory):
        """Constructor

        SQLite database will be created in panseq_directory

        Args:
            panseq_directory(str): Panseq results directory
           
           
        """

        self._panseqdir = panseq_directory
        self._dbfile = os.path.join(self._panseqdir, 'panseqannot.db')
        self._dbconn = sqlite3.connect(self._dbfile)
        

    def init_db(self):
        """Initialize database

        Initialize SQLite database for lookup of fragments

        Args:
            panseq_directory(str): Panseq results directory
   
        """

        # Create tables
        c = self._dbconn.cursor()
        c.execute('DROP TABLE IF EXISTS fragment')
        c.execute('''
            CREATE TABLE fragment
            (locus_id INT, name TEXT, contig TEXT, st INT, en INT)
            ''')
        c.execute("CREATE UNIQUE INDEX fragment_locus_id_idx ON fragment ('locus_id')")
        c.execute("CREATE UNIQUE INDEX fragment_position_idx ON fragment ('contig','st','en')")
        self._dbconn.commit()
        c.execute('DROP TABLE IF EXISTS annotation')
        c.execute('''
            CREATE TABLE annotation
            (annotation_id INTEGER PRIMARY KEY,
             type TEXT,
             description TEXT,
             contig TEXT,
             st INT,
             en INT,
             ln INT)
            ''')
        c.execute("CREATE INDEX annotation_type_id_idx ON annotation('type')")
        c.execute("CREATE INDEX annotation_position_id_idx ON annotation('contig', 'st', 'en')")
        self._dbconn.commit()
        c.execute('DROP TABLE IF EXISTS fragment_annotation')
        c.execute('''
            CREATE TABLE fragment_annotation
            (fragment_annotation_id INTEGER PRIMARY KEY,
             locus_id REFERENCES fragment('locus_id'),
             annotation_id REFERENCES annotation('annotation_id'),
             st INT, -- relative to the fragment sequence
             en INT, -- relative to the fragment sequence
             ln INT, -- relative to the fragment sequence
             coverage REAL)
            ''')
        c.execute("CREATE INDEX fragment_annotation_coverage_idx ON fragment_annotation('coverage')")
        c.execute("CREATE INDEX fragment_annotation_locus_id_idx ON fragment_annotation('locus_id')")
        c.execute("CREATE INDEX fragment_annotation_annotation_id_idx ON fragment_annotation('annotation_id')")
        self._dbconn.commit()

        # Load data
        pan_genome_file = os.path.join(self._panseqdir, 'pan_genome.txt')

        df = pandas.read_table(pan_genome_file, dtype={'LocusID': str})
        df = df[['LocusID','LocusName']]
        df.drop_duplicates(inplace=True)

        # Insert into DB
        insertion_rows = []
        insertion_cmd = """
            INSERT INTO fragment ('locus_id', 'name', 'contig', 'st', 'en')
            VALUES (?,?,?,?,?)
            """

        for idx, row in df.iterrows():
            locus = row['LocusID']
            name = row['LocusName']
        
            m = re.search(r'^(?P<contig>.+)_\((?P<begin>\d+)\.\.(?P<end>\d+)\)$', str(name))

            if not m:
                raise Exception("Invalid format for contig name: <{}>".format(name))

            if int(m.group('begin')) <= int(m.group('end')):
                start = m.group('begin')
                stop = m.group('end')
            else:
                start = m.group('end')
                stop = m.group('begin')

            insertion_rows.append((locus, name, m.group('contig'), start, stop))

            if insertion_rows > 100:
                c.executemany(insertion_cmd, insertion_rows)
                insertion_rows = []

        if insertion_rows:
            c.executemany(insertion_cmd, insertion_rows)
            insertion_rows = []

        self._dbconn.commit()


    def get_loci_by_location(self, contig, start, end):
        """Retrieve loci that overlap given region on contig

        Args:
            contig(str): Contig ID
            start(int): Start position of region
            end(int): End position of region

        Returns:
            List of lists:
                [0]: locus_id
                [1]: start of annotion on fragment
                [2]: end of annotion on fragment
                [3]: length of annotation overlap
                [4]: proportion of fragment covered by annotation
        
        """

        c = self._dbconn.cursor()

        c.execute("""
            SELECT locus_id, st, en  
            FROM fragment 
            WHERE contig = ? AND 
                st <= ? AND en >= ?
            """, (contig, end, start))

        # Compute overlaps
        results = []
        for row in c:

            if row[1] <= start:
                st = start - row[1] + 1
            else:
                st = 1
            if row[2] >= end:
                en = end - row[1] + 1
            else:
                en = row[2] - row[1] + 1

            ln = en - st + 1.0
            full_ln = row[2] - row[1] + 1
        
            results.append([row[0], st, en, ln, ln/full_ln])

        return results


    def load_annotation(self, annotation, annot_type, contig, start, end):
        """Insert genome region annotation and link to overlapping
        fragments in fragment_annotation table

        Args:
            annotation(str): Annotation assignment in text format (often JSON)
            annot_type(str): resfams|rgi|resfinder
            contig(str): Contig ID
            start(int): Start position of region
            end(int): End position of region

        Returns:
            List of fragment annotations
        
        """

        # Insert annotation
        c = self._dbconn.cursor()
        ln = end-start+1
        c.execute('INSERT INTO annotation (type, description, contig, st, en, ln) VALUES (?, ?, ?, ?, ?, ?)', 
            [annot_type, annotation, contig, start, end, ln])
        annotation_id = c.lastrowid

        # Retrieve fragments that overlap
        overlapping_fragments = self.get_loci_by_location(contig, start, end)
        for row in overlapping_fragments:
            row.append(annotation_id)

        # Link annotation to fragments
        c.executemany("""
            INSERT INTO fragment_annotation ('locus_id','st', 'en', 'ln', 'coverage', 'annotation_id')
            VALUES (?,?,?,?,?,?)
            """, overlapping_fragments)

        self._dbconn.commit()

        return overlapping_fragments


    def get_annotations(self, locus_id, annot_type=None):
        """Retrieve annotations for a locus ID

        Args:
            locus_id(str): Fragment locus_id
            annot_type(str)[OPTIONAL]: resfams|rgi|resfinder
        
        Returns:
            Returns dictionary of lists, one per requested annotation
            type. If annot_type is None, all annotation types are returned

            Empty dictionary is returned for fragments with no fragments
        
        """



        c = self._dbconn.cursor()
        bindings = [int(locus_id)]
        query = """
            SELECT a.type, a.description
            FROM annotation a, fragment_annotation f
            WHERE a.annotation_id = f.annotation_id AND
              f.locus_id = ?
            """
        if annot_type:
            if annot_type != 'rgi' and annot_type != 'resfinder' and annot_type != 'resfams':
                raise Exception('Invalid annot_type parameter: {}'.format(annot_type))
            query = query + ' AND a.type = ?'
            bindings.append(annot_type)

        c.execute(query, bindings)
        
        results = defaultdict(list)
        for r in c:
            results[r[0]].append(r[1])

        return results


    def size(self):
        c = self._dbconn.cursor()
        c.execute('''
            SELECT count(*) FROM fragment
            ''')

        return c.fetchone()[0]


if __name__ == "__main__":
    """Build database

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('annotate.annot')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Panseq directory')
  
    options = parser.parse_args()

    ann = PanseqAnnot(options.input)
    ann.init_db()

    logger.info("Inserted {} rows".format(ann.size()))



