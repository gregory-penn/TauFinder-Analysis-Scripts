#!/usr/bin/env python3

import numpy as np
import ROOT
from ROOT import TH1F, TFile
import subprocess
from argparse import ArgumentParser

from helpers.main import process_set

if __name__ == "__main__":

    parser = ArgumentParser()
    # BIB file 1
    # parser.add_argument('-i', '--inputFileDir', type=str, default='/data/gpenn/v8_noFragRem/pion_1GeV_1TeV/reco_output_sim_piguns_v8')
    # no BIB file 1
    # parser.add_argument('-i', '--inputFileDir', type=str, default='/data/gpenn/v8_noFragRem/pion_1GeV_1TeV_noBIB/reco_noFragRem_noBIB_output_sim_piguns_v8')
    # BIB file 2 
    #parser.add_argument('-i', '--inputFileDir', type=str, default='/data/gpenn/v8_noFragRem/round2_pion_1GeV_1TeV/reco_output_sim_piguns_1G_1T')
    # no BIB file 2
    parser.add_argument('-i', '--inputFileDir', type=str, default='/data/gpenn/v8_noFragRem/round2_pion_1GeV_1TeV_noBIB/reco_noFragRem_noBIB_output_sim_piguns_1G_1T')
    parser.add_argument('-o', '--outputFile', type=str, default='data/pi_no_bib_file2')
    args = parser.parse_args()

    MAX_EVENTS = 500

    # running over partial slices to minimize memory usage (still unadvised)
    index1 = list(range(0,10))
    index2 = list(range(0,10))

    # process individual slices
    for ind1 in index1:
        for ind2 in index2:
            print("Processing slice: ", ind1, ind2)
            pattern = f"{args.inputFileDir}.{ind1}{ind2}*.slcio"
            outFile = f"{args.outputFile}.{ind1}{ind2}.root"
            # processing only a hundredth of the files at a time for memory management purposes
            subprocess.run(["python", "helpers/process_reco.py", str(MAX_EVENTS), str(pattern), str(outFile)])

    # # print("Processing gen files")
    # hists = process_set(f"{args.inputFileDir}*.*.slcio", MAX_EVENTS)

    # output = TFile(args.outputFile, 'RECREATE')
    # for key, hist in hists.items():
    #     hist.Write()
    # output.Close()

    print("Done.")