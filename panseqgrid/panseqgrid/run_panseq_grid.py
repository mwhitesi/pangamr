#!/usr/bin/env python

"""Run panseq with all combinations of parameters


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
import subprocess

from panseqgrid.config import Config

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"

logger = None

def runner(configfile):
    """Run panseq

    Args:
        configfile(str): Location of panseq config file

    """

    subprocess.call(["panseq", configfile])


if __name__ == "__main__":
    """Build tree image

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('panseqgrid.run_panseq_grid')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Panseq config file')
  
    options = parser.parse_args()

    params = {
        'percentIdentityCutoff': [90,95],
        'fragmentationSize': [1000,500]
    }
    conf = Config(options.input, params)

    for c in conf.range():
        f = conf.write(c)
        logger.info("\n\n***Running panseq with config:\n{}\nFile: {}\n".format(c, f))
        runner(f)
        logger.info("\n\t...complete.\n")
        os.remove(f)


   
   
        

   
  

