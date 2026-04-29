#!/usr/bin/env python3
# meant to be run with access to a slices/ directory
import sys
import glob
import numpy as np
from pyLCIO import IOIMPL
from ROOT import TH1F, TFile
from helpers.main import process_set

import subprocess

COLLECTIONS = [
    "MCParticle",
]

NBINS = 10
THETA_RANGE = (-np.pi / 2, np.pi / 2)
PT_RANGE = (0, 1000)

def main():

    MAX_EVENTS = int(sys.argv[1]) if len(sys.argv) > 2 else -1
    pattern = sys.argv[2]
    output = sys.argv[3]
    print("pattern: ", pattern)
    hists = process_set(pattern, MAX_EVENTS) 

    output = TFile(output, 'RECREATE')
    for key, hist in hists.items():
        hist.Write()
    output.Close()


if __name__ == "__main__":
    main()    