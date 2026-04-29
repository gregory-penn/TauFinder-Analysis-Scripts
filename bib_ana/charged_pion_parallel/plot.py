import ROOT
from argparse import ArgumentParser
import os
from array import array
import mplhep as hep
from numpy import linspace

ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()
parser.add_argument('--bibFile', type=str, required=False)
parser.add_argument('--noBIBFile', type=str, default='data/bibFile1/pi_test.root')
parser.add_argument('--label1', type=str, default='BIB')
parser.add_argument('--label2', type=str, default='Non BIB')
parser.add_argument('--outdir', type=str, default='plots')
parser.add_argument('--charge', type=str, default='none', choices=['plus', 'minus', 'both', 'none'],
                    help='Charge state of the pions ("plus", "minus", "both", or "none" for no charge matching)')
parser.add_argument('--plotTypes', type=str, default='strict', choices=['default', 'strict', 'both'], help='Plot types to overlay')

args = parser.parse_args()

##################
# CHANGE BINS HERE
rebins = [0, 25, 50, 75, 100, 200, 300, 400, 700, 1000]
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
    'strict': [ROOT.kRed, ROOT.kGreen]
}
markers = {
    'both': [ROOT.kFullCircle, ROOT.kOpenCircle, ROOT.kFullSquare, ROOT.kOpenSquare],
    'default': [ROOT.kFullCircle, ROOT.kOpenCircle],
    'strict': [ROOT.kFullCircle, ROOT.kOpenCircle]
}
labels = {
    'both': [f'{args.label1} Default (0.25)', f'{args.label1} Strict (0.1)', f'{args.label2} Default (0.25)', f'{args.label2} Strict (0.1)'],
    'default': [f'{args.label1} Default (0.25)', f'{args.label2} Default (0.25)'],
    'strict': [f'{args.label1} Strict (0.1)', f'{args.label2} Strict (0.1)']
}
label = labels[to_include]
marker = markers[to_include]
color = colors[to_include]

# get charge pion we are interested in
charge = str(args.charge).lower()
latexs = {
    'plus': "Simulated #pi^{+} Gun",
    'minus': "Simulated #pi^{-} Gun",
    'both': "Simulated #pi^{+} and #pi^{-} Guns",
    'none': "Simulated Charged #pi Gun (no charge matching)"
}
regional_latex = {
    'barrel': "Barrel Region (0.70 < #theta < 2.44)",
    'centbarrel': "Central Barrel region (0.99 < #theta < 2.15)",
    'endcap': "Endcap Region (0.175 < #theta < 0.7 or 2.44 < #theta < 2.96)",
    'transition': "Transition Region (0.7 < #theta < 0.99 or 2.15 < #theta < 2.44)"
}
if charge not in latexs: raise ValueError("Invalid charge state. Must be 'plus', 'minus', or 'both'.")
plot_titles = {
    'plus': 'Pi+',
    'minus': 'Pi-',
    'both': 'Charged Pion',
    'none': 'Charged Pion (no charge matching)'
}
latex = latexs[charge]
plot_title = plot_titles[charge]

os.makedirs(args.outdir, exist_ok=True)

files = []
try:
    # f1 = ROOT.TFile(args.bibFile, 'READ')
    # files.append(f1)
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

    h.GetYaxis().SetRangeUser(0.0, 1.2)

    h.GetXaxis().SetTitleFont(62)
    h.GetXaxis().SetLabelFont(62)
    h.GetYaxis().SetTitleFont(62)
    h.GetYaxis().SetLabelFont(62)

# Eff plots
def make_eff(num, den, name, title, xtitle, rebin=False):
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
    eff.SetTitle(title)
    eff.GetXaxis().SetTitle(xtitle)
    eff.GetYaxis().SetTitle('#epsilon')
    eff.GetYaxis().SetRangeUser(0.0, 1.2)
    eff.SetStats(0)

    return eff


def draw_overlay(hists, name, title, xtitle, regional=None):

    c = ROOT.TCanvas('c', 'overlay', 800, 600)
    # Legend
    legend = ROOT.TLegend(0.63, 0.77, 0.9, 0.90)

    print("Length of hists:" , len(hists))
    for i in range(len(hists)):
        print("i:", i)
        hist = hists[i]
        if i == 0:
            hist.SetTitle(title)
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle('#epsilon')
            style_hist(hist, color[i], marker[i], 1)
            hist.Draw('E1')
            legend.AddEntry(hist, label[i], 'lp')
        else:
            style_hist(hist, color[i], marker[i], 1)
            hist.Draw('E1 SAME')
            legend.AddEntry(hist, label[i], 'lp')

    legend.SetBorderSize(1)
    legend.SetTextSize(0.03)
    legend.SetFillStyle(0)
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextAlign(13)
    text.SetTextSize(0.04)
    text.DrawLatex(0.12, 0.89, "#bf{#it{Muon Collider}}")
    text.SetTextSize(0.03)
    text.DrawLatex(0.12, 0.93, "#it{MAIA Detector Concept}")
    text.DrawLatex(0.12, 0.85, "#sqrt{s} = 10 TeV")
    if regional:
        text.DrawLatex(0.12, 0.82, regional_latex[regional])
        text.DrawLatex(0.12, 0.79, latex)
    else: text.DrawLatex(0.12, 0.82, latex)

    c.Update()

    c.SaveAs(os.path.join(args.outdir, name + '_overlay.png'))
    del c


# ---- Global plots ----
for i in range(len(eff_names)):
    var = eff_names[i]
    print("var: ", var)
    title = f'{plot_title} Efficiency vs {variables[i]}'
    xtitle = f'True {variables[i]}'
    if units[i]:
        xtitle += f' [{units[i]}]'
    
    for region in regions:
        to_overlay = []
        for file in files:
            num = file.Get(f'mc_matched_{var}_{region}')
            num_trkcls = file.Get(f'trk_cls_match_{var}_{region}')
            den = file.Get(f'mc_pion_{var}_{region}')

            # if var == "pt" and region == "barrel":
            #     for bin in range(1, num.GetNbinsX() + 1):
            #         print("Bin number: ", bin)
            #         print("Numerator content: ", num.GetBinContent(bin))
            #         print("Denominator content: ", den.GetBinContent(bin))

                #exit()

    
            print('type of num: ', type(num_trkcls), num_trkcls)
            print('type of num_trkcls: ', type(num_trkcls), num_trkcls)

            name = f"{var}_{region}"

            if var == 'pt':
                to_overlay.append(make_eff(num, den, f'{var}_eff', title, xtitle, rebin=True))
                print("length after appending once: ", len(to_overlay))
                to_overlay.append(make_eff(num_trkcls, den, f'{var}_eff_trk_cls', title, xtitle, rebin=True))
                print("length after appending: ", len(to_overlay))
            else:
                to_overlay.append(make_eff(num, den, f'{var}_eff', title, xtitle))
                # to_overlay.append(make_eff(num_trkcls, den, f'{var}_eff_trk_cls', title, xtitle))

        print(len(to_overlay), to_overlay)
        

        draw_overlay(to_overlay, name, title, xtitle)
        exit()


# # ---- Regional pT plots ----
# for r in regions:
#     name = f'pt_{r}'
#     title = f'{plot_title} Efficiency vs p_{{T}}'
#     xtitle = 'True p_{T} [GeV/c]'

#     to_overlay = []
#     for file in files:
#         den = file.Get(f'mc_pion_pt_{r}')

#         makeD, makeS = False, False
#         if to_include == 'both' or to_include == 'default':
#             numDefault, makeD = file.Get(f'mc_matched_pt_{r}'), True
#         if to_include == 'both' or to_include == 'strict':
#             numStrict, makeS = file.Get(f'mc_matched_pt_{r}_strict'), True

#         if makeD: to_overlay.append(make_eff(numDefault, den, f'{name}_eff', title, xtitle, rebin=True))
#         if makeS: to_overlay.append(make_eff(numStrict, den, f'{name}_eff_strict', title, xtitle, rebin=True))

#     draw_overlay(to_overlay, name, title, xtitle, regional=r)


print("Done. Overlays saved in:", args.outdir)
