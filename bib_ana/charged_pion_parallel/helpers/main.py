#!/usr/bin/env python3

import glob
import math
import numpy as np
import ROOT
from ROOT import TH1F, TFile, TCanvas, gPad
from pyLCIO import IOIMPL
from .geometry import eta, theta_region, delta_phi


# main processing function
# TODO: Clean some of the below parameters 
# they should be input args in some sense

allowed_pdgs = {
    'plus': {211},
    'minus': {-211},
    'both': {211, -211},
    'none': {211, -211}
}

# hard-coding a positively charged pion
charge = "none"
ignore_charge = True

# Hist setup
PT_MIN, PT_MAX, PT_BINS = 0.0, 1000.0, 160

ETA_MIN, ETA_MAX, ETA_BINS = 0.0, 2.5, 20

THETA_BINS = 20

PHI_BINS = 24

EFF_MIN, EFF_MAX = 0.0, 1.2

def book(h):
    h.SetDirectory(0)
    return h

# Regional pT histograms
regions = ['barrel', 'centbarrel', 'transition', 'endcap']
regional_maxs = {
    'plus': [800, 600, 400, 200],
    'minus': [800, 600, 400, 200],
    'both': [1600, 1200, 800, 400],
    'none': [800, 600, 400, 200]
}

regional_max = regional_maxs[charge]

def process_set(pattern, max_events):

    hists = {}

    # booking a ton of TH1Fs
    for region in regions:
        # initializing ROOT TH1Fs for filling    
        # Histograms for counting number of MCPs 
        hists[f"fMCPt_{region}"] = book(TH1F(f'mc_pion_pt_{region}', f'MC Charged Pion p_{{T}} ({region})', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fMCTheta_{region}"] = book(TH1F(f'mc_pion_theta_{region}', f'MC Charged Pion #theta ({region});#theta_reco [rad];Entries', THETA_BINS, 0, np.pi))

        # # Histograms for counting tracking efficiency
        # hists[f"fTrackPt_{region}"] = book(TH1F(f'matched_track_pt_{region}', f'Matched Track p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        # hists[f"fTrackTheta_{region}"] = book(TH1F(f'matched_track_theta_{region}', f'Matched Track #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

        # # Histograms for counting track-cluster matching efficiency
        # hists[f"fTrkClsPt_{region}"] = book(TH1F(f'trk_cls_match_pt_{region}', f'Matched Track p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        # hists[f"fTrkClsTheta_{region}"] = book(TH1F(f'trk_cls_match_theta_{region}', f'Matched Track p_{{T}} #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

        # Histograms for counting reco charged pions
        hists[f"fMatchedPt_{region}"] = book(TH1F(f'mc_matched_pt_{region}', f'Matched Best Reco Charged Pion MC p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fMatchedTheta_{region}"] = book(TH1F(f'mc_matched_theta_{region}', f'Matched Best Charged Reco Pion MC #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

    files = sorted(glob.glob(pattern))
    selected_pdgs = allowed_pdgs[charge]
    event_count = 0

    (f"Found {len(files)} files")

    for fname in files:

        if event_count >= max_events:
            break

        reader = IOIMPL.LCFactory.getInstance().createLCReader()
        print("fname: ", fname)
        reader.setReadCollectionNames(['MCParticle', 'PandoraPFOs'])
        reader.open(fname)

        # try:
        #     while True:

        if event_count >= max_events:
            break

        evt = reader.readNextEvent()
        event_count += 1
        print("Event count: ", event_count)

        # removing loop over MCPs - the first one should always be the true one
        # I suspect others having (hopefully slightly) higher pT implies a material interaction resulting in a transverse boost
        # strange, but is very rare
        # for mc in mcs:
        #     if mc.getPDG() not in selected_pdgs: continue

        #     mcMom = mc.getMomentum()
        #     mcPt = math.hypot(mcMom[0], mcMom[1])

        #     if mcPt > best_mc_pt:
        #         best_mc_pt = mcPt
        #         best_mc = mc

        mcs = evt.getCollection('MCParticle')
        # Best MC charged pion
        best_mc = mcs[0]
        if best_mc is None:
            continue
        mcPDG = best_mc.getPDG()
        mcMom = best_mc.getMomentum()
        mcPt = math.hypot(mcMom[0], mcMom[1])
        mcTheta = math.acos(mcMom[2] / math.sqrt(mcPt**2 + mcMom[2]**2))
        mcEta = eta(mcTheta)
        mcPhi = math.atan2(mcMom[1], mcMom[0])
        mcE = best_mc.getEnergy()

        # find what region the MCP is in
        regs = theta_region(mcTheta)

        # fill histograms according to region
        if regs:
            for reg in regs:
                hists[f"fMCPt_{reg}"].Fill(mcPt)
                hists[f"fMCTheta_{reg}"].Fill(mcTheta)

        # initialize reco pis, to be found
        best_reco_pi = None
        best_reco_pi_pt = -1.0
        best_reco_pi_theta = -1.0
        
        pfos = evt.getCollection('PandoraPFOs')

        for pfo in pfos:
            if ignore_charge and abs(pfo.getType()) != abs(mcPDG): continue # no charge matching case
            elif pfo.getType() != mcPDG: continue

            # Pion kinematics
            recoPiMomDefault = pfo.getMomentum()
            recoPiPtDefault = math.hypot(recoPiMomDefault[0], recoPiMomDefault[1])
            pfoE = pfo.getEnergy()
            recoPiTheta = math.acos(recoPiMomDefault[2] / math.sqrt(recoPiPtDefault ** 2 + recoPiMomDefault[2] ** 2))
            recoPiEta = eta(recoPiTheta)
            recoPiPhi = math.atan2(recoPiMomDefault[1], recoPiMomDefault[0])

            # dR matching
            dphi = delta_phi(mcPhi, recoPiPhi)
            dR = math.sqrt(dphi*dphi + (mcEta - recoPiEta)**2)

            if dR < 0.1:
                if recoPiPtDefault > best_reco_pi_pt:
                    best_reco_pi_pt = recoPiPtDefault
                    best_reco_pi_theta = recoPiTheta
                    best_reco_pi = pfo

        if best_reco_pi is None: # check if no default match was found, if so, skip this event
            continue

        # fill histograms according to region
        if regs:
            for reg in regs:
                hists[f"fMatchedPt_{reg}"].Fill(mcPt)
                hists[f"fMatchedTheta_{reg}"].Fill(mcTheta)

        # except Exception:
        #     # EOF reached
        #     pass

        del evt
        del mcs
        del best_mc
        del pfos
        del best_reco_pi

    return hists
