#!/usr/bin/env python3

import glob
import math
import numpy as np
import ROOT
from ROOT import TH1F, TFile
import pyLCIO
from pyLCIO import IOIMPL
from .geometry import eta, theta_region, delta_phi
from .track_truth_match import build_rel_nav, system_to_relname


# main processing function
# TODO: Clean some of the below parameters 
# they should be input args in some sense

hit_collection_names = [
    "VBTrackerHitsConed",
    "VETrackerHitsConed",
    "IBTrackerHitsConed",
    "IETrackerHitsConed",
    "OBTrackerHitsConed",
    "OETrackerHitsConed",
    "VBTrackerHitsRelationsConed",
    "VETrackerHitsRelationsConed",
    "IBTrackerHitsRelationsConed",
    "IETrackerHitsRelationsConed",
    "OBTrackerHitsRelationsConed",
    "OETrackerHitsRelationsConed",
    "VertexBarrelCollectionConed",
    "VertexEndcapCollectionConed",
    "InnerTrackerBarrelCollectionConed",
    "InnerTrackerEndcapCollectionConed",
    "OuterTrackerBarrelCollectionConed",
    "OuterTrackerEndcapCollectionConed"
]

# hit_collection_names = [
#     "VBTrackerHits",
#     "VETrackerHits",
#     "IBTrackerHits",
#     "IETrackerHits",
#     "OBTrackerHits",
#     "OETrackerHits",
#     "VBTrackerHitsRelations",
#     "VETrackerHitsRelations",
#     "IBTrackerHitsRelations",
#     "IETrackerHitsRelations",
#     "OBTrackerHitsRelations",
#     "OETrackerHitsRelations",
#     "VertexBarrelCollection",
#     "VertexEndcapCollection",
#     "InnerTrackerBarrelCollection",
#     "InnerTrackerEndcapCollection",
#     "OuterTrackerBarrelCollection",
#     "OuterTrackerEndcapCollection"
# ]

# copying this in from the tracking script although I don't think it works
hit_collection_mask = {key:True for key in hit_collection_names}

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
        hists[f"fTrackPt_{region}"] = book(TH1F(f'matched_track_pt_{region}', f'Matched Track p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fTrackTheta_{region}"] = book(TH1F(f'matched_track_theta_{region}', f'Matched Track #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

        # Histograms for counting track-cluster matching efficiency
        hists[f"fTrkClsPt_{region}"] = book(TH1F(f'trk_cls_match_pt_{region}', f'Matched Track p_{{T}} ({region});p_{{T}}_true;Entries', PT_BINS, PT_MIN, PT_MAX))
        hists[f"fTrkClsTheta_{region}"] = book(TH1F(f'trk_cls_match_theta_{region}', f'Matched Track p_{{T}} #theta ({region});#theta_true [rad];Entries', THETA_BINS, 0, np.pi))

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
        collection_names = ['MCParticle', 'PandoraPFOs', 'SelectedTracks']#, 'MCParticle_SelectedTracks']# + hit_collection_names
        #reader.setReadCollectionNames(collection_names)
        reader.open(fname)

        evt = reader.readNextEvent()
        event_count += 1
        print("Event count: ", event_count)

        if "MCParticle" not in evt.getCollectionNames():
            # there seems to be an issue with some non-BIB files. The total number is small, so skipping them should be fine.
            print("Event seems bugged! Skipping.")
            continue

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

        tracks = evt.getCollection('SelectedTracks')
        relationsContainer = evt.getCollection('MCParticle_SelectedTracks')
        relation = pyLCIO.UTIL.LCRelationNavigator(relationsContainer)
        related_tracks = relation.getRelatedToObjects(best_mc)
        print("number of relation tracks: ", len(related_tracks))
        # auto-continue if there are no tracks or truth-matched tracks in the event
        if len(tracks) == 0 or len(related_tracks) == 0:
            continue

        #### Below is for proper track truth matching. Currently commented out. See comments for why. ####
        
        # # build relation between hit collections and sub-detector
        # rel_nav = build_rel_nav(evt)

        # hit_collections = []
        # for hname in hit_collection_names:
        #     if(not hit_collection_mask[hname]):
        #         print("I should never hit this, right??")
        #         continue
        #     try:
        #         hit_collections.append(evt.getCollection(hname))
        #     except:
        #         hit_collection_mask[hname] = False
        #         print('\tDid not find hit collection: {}. Disabling...'.format(hname))
        #         pass        

        # apparently SelectedTracks don't have any hits associated
        # this is a known bug
        # the tracking results assume that these tracks have high hit purity
        # therefore they will pass the truth-matching requirement
        # therefore take any event with a SelectedTrack as passing tracking requirements
        # for track in tracks:
        #     print("Looping over track")
        #     print("Track has number of hits:", len(track.getTrackerHits()))
        #     print("track omega: ", track.getOmega())
        #     truth_matched_hits = 0            
        #     for hit in track.getTrackerHits():
        #         position = hit.getPosition()
        #         print("hit pos:", position)
        #         encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
        #         decoder = pyLCIO.UTIL.BitField64(encoding)
        #         cellID = int(hit.getCellID0())
        #         decoder.setValue(cellID)
        #         detector = decoder["system"].value()
        #         layer = decoder['layer'].value()
        #         if detector == 1 or detector == 2:
        #             LC_pixel_nhit += 1
        #         if detector == 3 or detector == 4:
        #             LC_inner_nhit += 1
        #         if detector == 5 or detector == 6: 
        #             LC_outer_nhit += 1
        #         print("detector:", detector)
        #         print("len(getrel objects):", len(rel_nav[system_to_relname[detector]].getRelatedToObjects(hit)))
        #         for sim_hit in rel_nav[system_to_relname[detector]].getRelatedToObjects(hit):
        #             print(type(sim_hit))
        #             print(sim_hit)
        #             mcp_true = sim_hit.getMCParticle()
        #             if mcp_true and abs(mcp_true.getPDG()) == 211:
        #                 truth_matched_hits += 1
        #         truth_hit_ratio = truth_matched_hits / len(track.getTrackerHits())
        #         print("truth_hit_ratio: ", truth_hit_ratio)
            
        #### end track truth matching section ####

        # this satisfies our tracking efficiency requirements
        # fill tracking efficiency plots
        if regs:
            for reg in regs:
                hists[f"fTrackPt_{reg}"].Fill(mcPt)
                hists[f"fTrackTheta_{reg}"].Fill(mcTheta)

        # initialize reco pis, to be found
        best_reco_charged = None
        best_reco_charged_pt = -1.0
        
        pfos = evt.getCollection('PandoraPFOs')

        for pfo in pfos:
            # allowing all charged particles, to filter by charged pions later
            if abs(pfo.getType()) != abs(mcPDG) and abs(pfo.getType()) != 11 and abs(pfo.getType()) != 13: continue # no charge matching case
            # Pion kinematics
            recoChargedMomDefault = pfo.getMomentum()
            recoChargedPtDefault = math.hypot(recoChargedMomDefault[0], recoChargedMomDefault[1])
            recoChargedTheta = math.acos(recoChargedMomDefault[2] / math.sqrt(recoChargedPtDefault ** 2 + recoChargedMomDefault[2] ** 2))
            recoChargedEta = eta(recoChargedTheta)
            recoChargedPhi = math.atan2(recoChargedMomDefault[1], recoChargedMomDefault[0])

            # dR matching
            dphi = delta_phi(mcPhi, recoChargedPhi)
            dR = math.sqrt(dphi*dphi + (mcEta - recoChargedEta)**2)

            if dR < 0.1:
                if recoChargedPtDefault > best_reco_charged_pt:
                    best_reco_charged_pt = recoChargedPtDefault
                    best_reco_charged = pfo

        if best_reco_charged is None: # check if no default match was found, if so, skip this event
            continue

        # fill histograms according to region
        if regs:
            for reg in regs:
                hists[f"fTrkClsPt_{reg}"].Fill(mcPt)
                hists[f"fTrkClsTheta_{reg}"].Fill(mcTheta)
                # now add to the charged pion ID histogram
                if abs(best_reco_charged.getType()) == abs(mcPDG):
                    hists[f"fMatchedPt_{reg}"].Fill(mcPt)
                    hists[f"fMatchedTheta_{reg}"].Fill(mcTheta)

        del evt
        del mcs
        del best_mc
        del tracks
        del relationsContainer
        del related_tracks
        del pfos
        del best_reco_charged

    return hists
