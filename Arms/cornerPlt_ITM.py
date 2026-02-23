"""Thin wrapper: ITM corner plot via cornerPlt.make_corner().

Usage::

    python cornerPlt_ITM.py Data/ITM/ITM_MC.hdf5
"""

import sys
from cornerPlt import make_corner

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python cornerPlt_ITM.py <MC_output.hdf5>")
        sys.exit(1)
    make_corner(sys.argv[1], mirror_type='ITM')
