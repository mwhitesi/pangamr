#!/usr/bin/env python

"""Multiple Sequence Alignment utilities.

Tools for building and manipulating alignments

"""

from subprocess import check_output, CalledProcessError, STDOUT

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class SeqAligner(object):
    """Run MSA program.

    Wrapper around the multiple sequence alignment program

    """

    def __init__(self):
        """Constructor

        Args:
            NA

        """

        self._mafft = 'mafft'
        self._mafft_args = {
            'auto': '--auto',
            'add': '--add'
        }

        self._trimal = 'trimal'
        self._trimal_args = {
            'auto': '-gappyout -in',
            'html': '-htmlout',
            'matrix': '-sident',
            'output': '-out'
        }


    @property
    def aligner(self):
        """str: Path to executable."""

        return self._mafft

    @property
    def trimmer(self):
        """str: Path to executable."""

        return self._trimal


    @property
    def aligner_args(self):
        """str: Command-line args."""

        return self._mafft_args['auto']

    @property
    def aligner_add_args(self):
        """str: Command-line args."""
        
        return self._mafft_args['add']

    @property
    def trim_args(self):
        """str: Command-line args."""
        
        return self._trimal_args['auto']

    @property
    def trim_output_args(self):
        """str: Command-line args."""
        
        return self._trimal_args['output']

   
    def align(self, fasta_file, alignment_file):
        """Run aligner program

        Args:
            fasta_file (str): Filepath to input fasta file
            alignment_file (str): Filepath for output from alignment


        Returns:
            None

        """

        cmd_args = self.aligner_args
        cmd = "{} {} {} > {}".format(self.aligner, cmd_args, fasta_file, alignment_file)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "Aligner failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None


    def add(self, fasta_file, existing_alignment_file, output_alignment_file):
        """Add new sequences to existing alignment

        Args:
            fasta_file (str): Filepath to input fasta file
            existing_alignment_file (str): Filepath to aligned fasta file
            output_alignment_file (str): Filepath for output from alignment

        Returns:
            None

        """
        cmd_args = self.aligner_add_args
        cmd = "{} {} {} {} > {}".format(self.aligner, cmd_args, fasta_file, existing_alignment_file, 
            output_alignment_file)

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "Aligner --add failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None


    def trim(self, alignment_file, trimmed_file, 
        trimming_summary_file=False, ident_matrix_file=False):
        """Trim alignment using external program trimal

        Args:
            alignment_files (str): Filepath to aligned fasta file
            trimmed_files (str): Filepath for output with trimmed alignment
            trimming_summary_file (str)[OPTIONAL]: Filepath for HTML output showing trimmed 
                alignment columns
            ident_matrix_file (str)[OPTIONAL]: Filepath for output figure 
                showing pairwise sequence identities

        Returns:
            None

        """
        
        cmd_args = self.trim_args
        output_args = "{} {}".format(self.trim_output_args, trimmed_file)
        cmd = "{} {} {} {}".format(self.trimmer, cmd_args, alignment_file, 
            output_args)

        if trimming_summary_file:
            cmd += " {} {}".format(self._trimal_args['html'], trimming_summary_file)

        if ident_matrix_file:
            cmd += " {} > {}".format(self._trimal_args['matrix'], ident_matrix_file)

        #print cmd

        try:
            check_output(cmd, stderr=STDOUT, shell=True, universal_newlines=True)                         
        except CalledProcessError as e:
            msg = "Alignment trimming failed: {} (return code: {}).".format(e.output, e.returncode)                                                                                                   
            raise Exception(msg)

        None


