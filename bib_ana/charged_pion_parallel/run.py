#!/usr/bin/env python3

import numpy as np
import ROOT
from ROOT import TH1F, TFile, TCanvas, gPad
import subprocess
from argparse import ArgumentParser

from helpers.main import process_set

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('-i', '--inputFileDir', type=str, default='/data/gpenn/v8_noFragRem/pion_1GeV_1TeV_noBIB/')
    parser.add_argument('-o', '--outputFile', type=str, default='data/pi_bib_ana.root')
    args = parser.parse_args()

    MAX_EVENTS = 500

    # running over partial slices to minimize memory usage (still unadvised)
    indices = list(range(0,10))

    # process individual slices
    # for index in indices:
    #     print("Processing slice: ", index)
    #     # processing only a tenth of the files at a time for memory management purposes
    #     subprocess.run(["python", "helpers/process_reco.py", str(index), str(MAX_EVENTS)])

    # print("Processing gen files")
    hists = process_set(f"{args.inputFileDir}*.*.slcio", MAX_EVENTS)

    output = TFile(args.outputFile, 'RECREATE')
    for key, hist in hists.items():
        hist.Write()
    output.Close()

    print("Done.")