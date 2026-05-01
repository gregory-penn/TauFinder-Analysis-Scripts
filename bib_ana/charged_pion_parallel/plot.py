import ROOT
from argparse import ArgumentParser
import os
from array import array
import mplhep as hep
from numpy import linspace

ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()
parser.add_argument('--bibFile', type=str, default='data/total_bib/total_bib.root')
parser.add_argument('--noBIBFile', type=str, default='data/noBIB/pi_total_nobib.root')
parser.add_argument('--label1', type=str, default='BIB')
parser.add_argument('--label2', type=str, default='Non BIB')
parser.add_argument('--outdir', type=str, default='plots')
parser.add_argument('--charge', type=str, default='none', choices=['plus', 'minus', 'both', 'none'],
                    help='Charge state of the pions ("plus", "minus", "both", or "none" for no charge matching)')
parser.add_argument('--plotTypes', type=str, default='strict', choices=['default', 'strict', 'both'], help='Plot types to overlay')

args = parser.parse_args()

##################
# CHANGE BINS HERE
# rebins = [0, 50, 100, 200, 300, 400, 700, 1000]
rebins = [0, 50, 100, 300, 500, 700, 1000]
# rebins = [0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
#rebins = [0, 25, 50, 75, 100, 200, 300, 400, 500, 525, 550, 570, 600, 625, 650, 675, 700, 725, 750, 775, 800, 900, 1000]
#rebins = linspace(0, 1000, 81)
n_rebins = len(rebins) - 1
##################

# plots to overlay
to_include = str(args.plotTypes).lower()
if to_include not in ['default', 'strict', 'both']: raise ValueError("Invalid plot type. Must be 'default', 'strict', or 'both'.")
colors = {
    'both': [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow],
    'default': [ROOT.kRed, ROOT.kGreen],
    'strict': [ROOT.kP6Yellow, ROOT.kP6Red, ROOT.kP6Blue]
}
markers = {
    'both': [ROOT.kFullCircle, ROOT.kOpenCircle, ROOT.kFullSquare, ROOT.kOpenSquare],
    'default': [ROOT.kFullCircle, ROOT.kOpenCircle],
    'strict': [ROOT.kOpenCircle, ROOT.kOpenCircle, ROOT.kFullCircle]
}
labels = {
    'both': [f'{args.label1} Default (0.25)', f'{args.label1} Strict (0.1)', f'{args.label2} Default (0.25)', f'{args.label2} Strict (0.1)'],
    'default': [f'{args.label1} Default (0.25)', f'{args.label2} Default (0.25)'],
    'strict': ["Tracking", "Track-Cluster", "Identification"]
}
label = labels[to_include]
marker = markers[to_include]
color = colors[to_include]

# get charge pion we are interested in
charge = str(args.charge).lower()
regional_latex = {
    'barrel': "0.70 < #it{#theta} < 2.44",
    'centbarrel': "0.99 < #it{#theta} < 2.15",
    'endcap': "0.175 < #it{#theta} < 0.7 or 2.44 < #it{#theta} < 2.96",
    'transition': "0.7 < #it{#theta} < 0.99 or 2.15 < #it{#theta} < 2.44"
}
plot_titles = {
    'plus': 'Pi+',
    'minus': 'Pi-',
    'both': 'Charged Pion',
    'none': 'Charged Pion (no charge matching)'
}
plot_title = plot_titles[charge]

os.makedirs(args.outdir, exist_ok=True)

files = []
try:
    f1 = ROOT.TFile(args.bibFile, 'READ')
    files.append(f1)
    f2 = ROOT.TFile(args.noBIBFile, 'READ')
    files.append(f2)
except IOError:
    print('Could not open one or more ROOT files!')
    exit()

# Global efficiency histogram names
eff_names = ['pt','theta']
variables = ['p_{T}', '#theta']
units = ['GeV', 'rad']

regions = ['barrel', 'centbarrel', 'transition', 'endcap']

def style_hist(h, color, marker, linestyle=1):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetLineWidth(2)
    h.SetMarkerStyle(marker)
    h.SetMarkerSize(1.3)
    h.SetLineStyle(linestyle)
    h.SetStats(0)

    h.GetYaxis().SetRangeUser(0.0, 1.3)

    # h.GetXaxis().SetTitleFont(62)
    # h.GetXaxis().SetLabelFont(62)
    # h.GetYaxis().SetTitleFont(62)
    # h.GetYaxis().SetLabelFont(62)

# Eff plots
def make_eff(num, den, name, xtitle, rebin=False):
    if rebin:
        rebin_arr = array('d', rebins) # make rebin array into cpp array

        num_rebin = num.Rebin(n_rebins, num.GetName()+"_rebin", rebin_arr)
        den_rebin = den.Rebin(n_rebins, den.GetName()+"_rebin", rebin_arr)
        eff = num_rebin.Clone(name)
        eff.Divide(num_rebin, den_rebin, 1, 1, 'B')
    else:
        eff = num.Clone(name)
        eff.Divide(num, den, 1, 1, 'B')

    eff.SetLineWidth(2)
    # eff.SetTitle(title)
    eff.GetXaxis().SetTitle(xtitle)
    eff.GetYaxis().SetTitle('Efficiency')
    eff.GetYaxis().SetRangeUser(0.0, 1.2)
    eff.SetStats(0)

    return eff


def draw_overlay(hists, name, xtitle, regional=None, individual="", biblabel = ""):

    c = ROOT.TCanvas('c', 'overlay', 800, 600)
    # Legend

    if individual == "":
        legend = ROOT.TLegend(0.71, 0.75, 0.89, 0.88)        
        for i in range(len(hists)):
            hist = hists[i]
            if i == 0:
                hist.SetTitle("")
                hist.GetXaxis().SetTitle(xtitle)
                hist.GetYaxis().SetTitle('Efficiency')
                style_hist(hist, color[i], marker[i], 1)
                hist.Draw('E1')
                legend.AddEntry(hist, label[i], 'lp')
            else:
                style_hist(hist, color[i], marker[i], 1)
                hist.Draw('E1 SAME')
                legend.AddEntry(hist, label[i], 'lp')
    else: 
        legend = ROOT.TLegend(0.78, 0.75, 0.89, 0.88)        
        # this assumes the order
        histBIB = hists[0]
        histNoBIB = hists[1]

        # want to draw noBIB first
        histNoBIB.SetTitle("")
        histNoBIB.GetXaxis().SetTitle(xtitle)
        histNoBIB.GetYaxis().SetTitle(f'{individual} Efficiency')
        style_hist(histBIB, ROOT.kBlack, ROOT.kFullCircle, 1)
        style_hist(histNoBIB, ROOT.kGray, ROOT.kFullCircle, 1)
        legend.AddEntry(histNoBIB, "No BIB", 'lp')        
        legend.AddEntry(histBIB, "BIB", 'lp')        
        histNoBIB.Draw('E1')
        histBIB.Draw('E1 SAME')

    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetFillStyle(0)
    legend.Draw()
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextAlign(13)
    text.SetTextSize(0.04)
    text.DrawLatex(0.14, 0.87, "#bf{#it{Muon Collider}}")
    text.SetTextSize(0.03)
    text.DrawLatex(0.687, 0.935, "#it{MAIA} Detector Concept")
    text.DrawLatex(0.14, 0.80, "#sqrt{s} = 10 TeV")
    if regional:
        text.DrawLatex(0.14, 0.77, regional_latex[regional])
        text.DrawLatex(0.14, 0.83, f"Simulation {biblabel}(EU24 Lattice)")
    else: text.DrawLatex(0.14, 0.83, f"Simulation {biblabel}(EU24 Lattice)")

    c.Update()

    individual = individual.replace(" ", "")
    biblabel = biblabel.replace(" ", "_")

    c.SaveAs(os.path.join(args.outdir, f'{name}_{individual}_{biblabel}_overlay.png'))
    del c


# ---- Global plots ----
for i in range(len(eff_names)):
    var = eff_names[i]
    print("var: ", var)
    # title = f'{plot_title} Efficiency vs {variables[i]}'
    xtitle = f'True Pion {variables[i]}'
    if units[i]:
        xtitle += f' [{units[i]}]'
    
    for region in regions:
        # expects the first file to be BIB, the second to be without BIB
        # space after BIB is needed
        biblabel = ["with BIB ", "without BIB "]
        bibindex = -1
        trk_cls_match = []
        tracking_alone = []
        for file in files:
            to_overlay = []
            bibindex += 1
            num = file.Get(f'mc_matched_{var}_{region}')
            num_trkcls = file.Get(f'trk_cls_match_{var}_{region}')
            num_trk = file.Get(f'matched_track_{var}_{region}')
            den = file.Get(f'mc_pion_{var}_{region}')

            # if var == "pt" and region == "barrel":
            #     for bin in range(1, num.GetNbinsX() + 1):
            #         print("Bin number: ", bin)
            #         print("Numerator content: ", num_trk.GetBinContent(bin))
            #         print("Denominator content: ", den.GetBinContent(bin))

            #     exit()


            name = f"{var}_{region}"
            rebin = False
            if var == 'pt':
                rebin = True
            trackingE = make_eff(num_trk, den, f'{var}_eff_tracking', xtitle, rebin=rebin)
            to_overlay.append(trackingE)
            tracking_alone.append(trackingE)
            to_overlay.append(make_eff(num_trkcls, den, f'{var}_eff_trk_cls', xtitle, rebin=rebin))
            to_overlay.append(make_eff(num, den, f'{var}_eff', xtitle, rebin=rebin))
            trk_cls_match.append(make_eff(num_trkcls, num_trk, f'{var}_eff_trk_cls_alone', xtitle, rebin=rebin))

            draw_overlay(to_overlay, name, xtitle, regional = region, biblabel=biblabel[bibindex])

        print("Drawing direct comparison plots...")
        draw_overlay(trk_cls_match, name, xtitle, regional = region, individual = "Track-Cluster Matching")
        draw_overlay(tracking_alone, name, xtitle, regional = region, individual = "Tracking")

print("Done. Overlays saved in:", args.outdir)
