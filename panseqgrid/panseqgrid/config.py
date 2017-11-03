#!/usr/bin/env python

"""Module to manipulate Panseq conf files


"""

import os
import pandas as pd
import tempfile
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

        self._valid = ['percentIdentityCutoff', 'fragmentationSize', 'coreGenomeThreshold']
        self._tracking_parameters = { 'fragmentationSize': 'minimumNovelRegionSize' }

        self._load(base_config_file)

        # Create root directory if it does not exist
        if not self._conf['baseDirectory']:
            raise ValueError('Missing parameter: baseDirectory')
        else:
            if not os.path.exists(self._conf['baseDirectory']):
                os.makedirs(self._conf['baseDirectory'])

            self._base_directory = self._conf['baseDirectory']

        self._permute(permutations)


    def _load(self, base_config_file):
        """Load default config values from file

        Args:
            base_config_file(str): Panseq config file with default values

        """

        tmp = pd.read_table(base_config_file, header=None, names=['param','value'], 
            index_col=0, dtype={'param': str, 'value': str}).to_dict()

        self._conf = tmp['value']


    def _permute(self, permutations):
        """Create config combinations 

        Cartesean product of everything in dictionary

        Args:
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values
           
        """

        ok, bad_param = self._permutable(permutations)
        if ok:
            self._combinations = [dict(zip(permutations, v)) for v in product(*permutations.values())]

            for c in self._combinations:
                # Add unique base name
                param_string = ''
                param_string = '__'.join(['{}{}'.format(k,v) for k,v in c.items() ])
                filename = os.path.join(self._base_directory, param_string)

                c['baseDirectory'] = filename

                # Some parameters need to be set to same value
                for p in self._tracking_parameters:
                    if p in c:
                        c[self._tracking_parameters[p]] = c[p]

        else:
            raise ValueError('{} parameter cannot be permuted'.format(bad_param))



    def _permutable(self, permutations):
        """Check requested permute keywords

        Args:
            permutations(dict): dictionary of Panseq config keywords pointing to lists of config values

        Returns
            bool: True if everything checks out
           
        """

        for p in permutations:
            if not p in self._valid:
                return (False, p)

        return (True, None)


    def range(self):
        """Config generator

        Iterates over all config dictionaries in self._combinations

        Returns:
            dictionary
           
        """

        newconfig = { k: v for k,v in self._conf.items() }
        for c in self._combinations:
           
            for k,v in c.items():
                newconfig[k] = v

            yield newconfig


    def write(self, conf, filename=None):
        """Save config dictionary to file

        User is responsible for deleting file

        Args:
            conf(dict): Panseq config dictionary
            filename[OPTIONAL](str): File location to write to. If not provided a temporary file will
                be generated and returned.

        Returns:
            filename(str)
           
        """

        if not filename:
            fh = tempfile.NamedTemporaryFile(delete=False)
            filename = fh.name
            fh.close()

        with open(filename, mode='w') as outfh:
            for k, v in conf.items():
                outfh.write("{}\t{}\n".format(k,v))

        return filename
    
        
