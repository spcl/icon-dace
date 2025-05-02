#!/usr/bin/env python
"""Wrapper for FORD fortran documentation generator"""

import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__),'3rdparty/ford'))

from ford import run

if __name__ == '__main__' :
    sys.exit(run())
