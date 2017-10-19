#!/usr/bin/env python

"""Functions to display tree for a Panseq Loci


"""

import tempfile
import colorsys
import re
from collections import namedtuple
from ete3 import Tree, AttrFace, TextFace, TreeStyle

from panseqtrees.fasttree import FastTreeWrapper
from panseqtrees.seqaligner import SeqAligner


__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class LocusTree(object):
    """Tree representation for pangenome locus


    """

    def __init__(self, treestring, idmap, phenodict=None):
        """Constructor

        Args:
            treestring(str): Newick format tree string for ETE constructor
            idmap(dict): mapping between names in tree string to fasta header ids. Contains
                id2name and name2id dictonaries.
            phenodict(dict)[OPTIONAL]: Two-level dictionary of phenotypes. Each phenotype
                dictionary contains allele names (matching fasta headers) to phenotype values.

        """

        self._tree = Tree(treestring)

        self._id2node = dict()

        # Initialize tree properties
        for l in self._tree.iter_leaves():
            newid = parse_fasta_id(idmap['name2id'][l.name])
            l.add_feature('id', newid)
            self._id2node[str(newid)] = l
        
        # Parse phenotype
        self._parse(phenodict)


    def _parse(self, phenodict):
        """Load phenotype and other attributes into memory

        Args:
            phenodict(dict): two-level dictionary. Each phenotype 'class' dictionary has allele
                to value mapping. Allele keys must match fasta headers in sequence file.

        """

        self._pheno_na = 'NA'
        
        if phenodict:

            self._pheno_classes = dict()
            
            for p in phenodict.keys():

                self._pheno_classes[p] = set()

                for l in self._tree.iter_leaves():
                    gn = l.id.genome
                    if gn in phenodict[p]:
                        l.add_feature(p, phenodict[p][gn])

                        if not phenodict[p][gn] in self._pheno_classes[p]:
                            self._pheno_classes[p].add(phenodict[p][gn])
                    else:
                        l.add_feature(p, self._pheno_na)

                # Assign order to phenotype variants
                ordered_variants = dict()
                i = 0
                for v in sorted(self._pheno_classes[p]):
                    ordered_variants[v] = i
                    i+=1
                self._pheno_classes[p] = ordered_variants

            # Initialize visualization defaults
            self._pheno_colors = ncolors([len(self._pheno_classes[p]) for p in self._pheno_classes])
            self._pheno_na_color = 'black'

            # Assign colors to nodes
            for l in self._tree.iter_leaves():
                self._init_faces(l)

        else:
            self._pheno_classes = None


    def _init_faces(self, node):

        if node.is_leaf():
            c = 0
            for p in self._pheno_classes:
                f = AttrFace(p)
                pheno = getattr(node, p, self._pheno_na)
                if pheno == self._pheno_na:
                    f.background.color = self._pheno_na_color
                else:
                    f.background.color = self._pheno_colors[c][self._pheno_classes[p][pheno]]
                node.add_face(f, column=c, position="aligned")
                c+=1

            node.add_face(TextFace(str(node.id)[0:20]), column=0, position="branch-right")


    def render(self, file):
        """Generate tree image and save to file

        Args:
            file(str): File location

        """

        ts = TreeStyle()
        ts.show_leaf_name = False
        self._tree.render( file, w=200, units="mm", tree_style=ts)


def parse_fasta_id(id):
    """Extract Genome and gene identifier fasta header

    From the ID component of a fasta header, determine
    for a particular genome copy, whether there is multiple 
    alleles of a gene in the genome.  Format is:

    >[lcl|]genome|allele

    Args:
        id (str): Fasta ID

    Returns:
        AlleleID namedtuple

    """

    # Remove database id
    id = re.sub(r'^(?:lcl\|)|(?:gi\|)', '', id)

    # Extract genome and gene Ids if they exist
    parts = re.split(r'\|', id, maxsplit=1)
    
    return AlleleID(*parts)


# Factory for AlleleID class
# Stores genome and allele ID components
class AlleleID(namedtuple('AlleleID', ['genome', 'allele'])):
    __slots__ = ()
    def __new__(cls, genome, allele=None):
        # add default values
        return super(AlleleID, cls).__new__(cls, genome, allele)
    def __str__(self):
        # String
        allele = '|'+str(self.allele) if self.allele else ''
        return '{}{}'.format(self.genome, allele)


def ncolors(levels):
    """Generate len(levels) colors with darker-to-lighter variants
    for each color

    Args:
        levels(list): List with one element per required color. Element
            indicates number of variants for that color

    Returns:
        List of lists with variants for each rgb color hue

    """

    n = len(levels)

    def rgb2hex(rgb):
        return '#%02x%02x%02x' % rgb
    def hls2hex(h, l, s):
        return rgb2hex( tuple([int(x*255) for x in colorsys.hsv_to_rgb(h, l, s)]))

    colors = []

    for x in range(n):
        colors.append([])
        nlevels = levels[x]
        saturation_step = 1/(nlevels + 1)

        for l in range(nlevels):
            hsv_triple = (x*1.0/n, saturation_step*(l+1), 0.5)
        
            colors[x].append(hls2hex(*hsv_triple))
    
    return colors


def build_tree(seqdict, phenodict=None):
    """Create LocusTree object from sequence dictionary

    Args:
        seqdict(dict): Parsed fasta input into headers (keys) and sequences (values)
        phenodict(dict)[OPTIONAL]: Two-level dictionary of phenotypes. Each phenotype is
            dictionary of allele names (matching fasta headers) to phenotype values.  

    Returns:
        LocusTree

    """

    # Rename sequences
    id2name = { 'name2id': {},'id2name': {} }

    ft = FastTreeWrapper()
    sa = SeqAligner()

    with tempfile.NamedTemporaryFile(mode='w') as fastafh, tempfile.NamedTemporaryFile() as alnfh, tempfile.NamedTemporaryFile(mode='r') as treefh:

        # Write sequences to temporary files with safe genome names
        a = 0
        for name,seq in seqdict.items():

            allele_name = 'allele'+str(a)
            a += 1
            id2name['name2id'][allele_name] = name
            id2name['id2name'][name] = allele_name

            fastafh.write('>{}\n{}\n'.format(allele_name, seq))

        fastafh.flush()

        # Align and build tree
        sa.align(fastafh.name, alnfh.name)
        ft.build(alnfh.name, treefh.name, nt=True, fast=False)

        # Load tree string into ETE tree object
        return( LocusTree(str(treefh.read()).strip(), id2name, phenodict) )




    





