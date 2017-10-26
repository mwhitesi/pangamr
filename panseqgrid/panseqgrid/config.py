#!/usr/bin/env python

"""Module to manipulate Panseq conf files


"""

import pandas as pd
from itertools import product

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "mwhiteside@canada.ca"


class Config(object):
    """Generates Panseq config files

    Creates new Panseq config files for every combination of
    permuted value inserted into the base config file


    """

    def __init__(self, base_config_file, permutations):
        """Constructor

        Args:
            base_config_file(str): Panseq config file with default values
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values
           
        """

        self._load(base_config_file)

        #self._permute()


    def _load(self, base_config_file):
        """Load default config values from file

        Args:
            base_config_file(str): Panseq config file with default values

        """

        self.conf = pd.read_table(base_config_file, header=None, names=['param','value'], index_col=0, dtype={'param': str, 'value': str})

    def _permute(self, permutations):
        """Create config combinations 

        Cartesean product of everything in dictionary

        Args:
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values
           
        """

        self._combinations = [dict(zip(permutations, v)) for v in product(*permutations.values())]

        # Add unique file names

        # Add coordinated configs

    def _permutable(self, permutations):
        """Check requested permute keywords

        Args:
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values

        Returns
            bool: True if everything checks out
           
        """

    
        
class Runner(object):
    """Generates Panseq config files

    Creates new Panseq config files for every combination of
    permuted value inserted into the base config file


    """

    def __init__(self, base_config_file, permutations):
        """Constructor

        Args:
            base_config_file(str): Panseq config file with default values
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values
           
        """

        self._load(base_config_file)

        self._permute()


    def _load(self, base_config_file):
        """Load default config values from file

        Args:
            base_config_file(str): Panseq config file with default values

        """

        pd.read

    