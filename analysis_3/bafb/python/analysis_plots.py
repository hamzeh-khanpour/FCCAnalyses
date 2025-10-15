# analysis_plots.py
# Visualize AFB^b Stage-2 outputs for an exclusive Z->bb sample (no backgrounds).
# Produces:
#   (1) Signed cosθ histogram overlaid with best-fit template
#   (2) ω(|signed cosθ|) TGraphErrors if truth hists exist
#   (3) Response matrix heatmap if present
#   (4) Lightweight monitors if present (Q_primary, Q0, Q1, cos0, cos1, sprime)
#   (5) Numbers panel from JSON summary (A_FB^true, A_FB^meas, ω, N, nbins, provenance highlights)
#
# Defaults:
#   ROOT: outputs/bafb/stage2/afbb_stage2.root
#   JSON: outputs/bafb/stage2/results/AFBb_fit.json
#   OUT : outputs/bafb/stage2/plots
#
# Usage:
#   python analysis_plots.py --root outputs/bafb/stage2/afbb_stage2.root \
#                            --json outputs/bafb/stage2/results/AFBb_fit.json \
#                            --out  outputs/bafb/stage2/plots
#
# Dependencies: ROOT + Python stdlib only.

import os, json, math, argparse
import ROOT

def ensure_dir(p):
    if not os.path.isdir(p):
        os.makedirs(p, exist_ok=True)

def save_canvas(c, outdir, base):
    png = os.path.join(outdir, base + ".png")
    pdf = os.path.join(outdir, base + ".pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)
    print(f"[plots] saved: {png} , {pdf}")

def get_hist(f, key):
    h = f.Get(key)
    return h if (h and not isinstance(h, ROOT.TObject) and not h.IsZombie()) else f.Get(key)  # fallback
    # (ROOT sometimes returns a proxy; above ensures we try twice)

def make_numbers_pave(summary):
    pave = ROOT.TPaveText(0.13,0.62,0.52,0.88,"NDC")
    pave.SetFillColor(0)
    pave.SetBorderSize(1)
    Atrue  = summary.get("A_FB_true", None)
    sAtrue = summary.get("A_FB_true_err", None)
    Ameas  = summary.get("A_FB_meas", None)
    omega  = summary.get("omega_used", None)
    N      = summary.get("N", None)
    nbins  = summary.get("nbins", None)
    acc    = summary.get("acceptance_folded", False)

    if Atrue is not None and sAtrue is not None:
        pave.AddText(f"A_FB^true = {Atrue:.6f} \u00B1 {sAtrue:.6f}")
    if Ameas is not None:
        pave.AddText(f"A_FB^meas = {Ameas:.6f}")
    if omega is not None:
        pave.AddText(f"omega (used) = {omega:.4f}")
    if N is not None and nbins is not None:
        pave.AddText(f"N = {N}   nbins = {nbins}")
    pave.AddText(f"acceptance folded: {'yes' if acc else 'no'}")

    # provenance highlights if present
    prov = summary.get("provenance", {})
    kappa = prov.get("cfg_kappa", None)
    pmin  = prov.get("cfg_track_p_min", None)
    cmax  = prov.get("cfg_abs_costheta_max", None)
    pjmin = prov.get("cfg_jet_p_min", None)
    opp   = prov.get("cfg_opposite_sign", None)
    sps   = prov.get("cfg_sprime_cut", None)
    line = []
    if kappa is not None: line.append(f"kappa={float(kappa):.2f}")
    if pmin  is not None: line.append(f"p_min={float(pmin):.2f} GeV")
    if cmax  is not None: line.append(f"|cosθ|<{float(cmax):.2f}")
    if pjmin is not None: line.append(f"|p|>{float(pjmin):.0f} GeV")
    if opp   is not None: line.append(f"opp={int(opp)}")
    if sps   is not None: line.append(f"s'/s≥{float(sps):.2f}")
    if line:
        pave.AddText(",  ".join(line))
    return pave

def omega_graph_from_pair(h_wrong, h_tot):
    """Build TGraphErrors for ω(|x|) = wrong/total with binomial errors."""
    if not h_wrong or not h_tot:
        return None
    nb = h_tot.GetNbinsX()
    gr = ROOT.TGraphErrors()
    for i in range(1, nb+1):
        tot = float(h_tot.GetBinContent(i))
        wrg = float(h_wrong.GetBinContent(i))
        x   = h_tot.GetXaxis().GetBinCenter(i)
        if tot <= 0:
            continue
        w = wrg / tot
        err = math.sqrt(max(0.0, w*(1.0-w)/tot))
        ip = gr.GetN()
        gr.SetPoint(ip, x, w)
        gr.SetPointError(ip, 0.0, err)
    gr.SetTitle("Wrong-sign fraction vs |signed cos#theta|;|signed cos#theta|;#omega(|x|)")
    return gr

def draw_overlay(root_path, json_path, outdir):
    f = ROOT.TFile.Open(root_path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")
    h_data = get_hist(f, "AFBb/h_signed_cosTheta")
    h_fit  = get_hist(f, "AFBb/h_fit_overlay")
    if not h_data:
        raise RuntimeError("Missing required histogram: AFBb/h_signed_cosTheta")
    # load JSON summary
    with open(json_path, "r") as js:
        s = json.load(js)

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c_signed","A_FB^b — signed cos#theta", 900, 650)
    h_data.SetLineWidth(2)
    h_data.SetMarkerStyle(20)
    h_data.SetTitle("A_{FB}^{b} — signed cos#theta;cos#theta;Events")
    h_data.Draw("E1")
    if h_fit:
        h_fit.SetLineColor(ROOT.kRed+1)
        h_fit.SetLineWidth(3)
        h_fit.SetFillStyle(0)
        h_fit.Draw("HIST SAME")
        leg = ROOT.TLegend(0.62,0.72,0.88,0.88)
        leg.AddEntry(h_data, "Data", "lep")
        leg.AddEntry(h_fit,  "Best-fit template", "l")
        leg.Draw()
    pave = make_numbers_pave(s)
    pave.Draw()
    save_canvas(c, outdir, "AFBb_signed_fit")

    # optional monitors
    monitors = [
        ("AFBb/h_Q_primary", "AFBb_Q_primary"),
        ("AFBb/h_Q0",        "AFBb_Q0"),
        ("AFBb/h_Q1",        "AFBb_Q1"),
        ("AFBb/h_cos0",      "AFBb_cos0"),
        ("AFBb/h_cos1",      "AFBb_cos1"),
        ("AFBb/h_sprime",    "AFBb_sprime"),
    ]
    for key, base in monitors:
        h = get_hist(f, key)
        if not h: 
            continue
        cM = ROOT.TCanvas(f"c_{base}", base, 800, 600)
        h.SetLineWidth(2)
        h.Draw("HIST")
        save_canvas(cM, outdir, base)

    # ω(|x|) if truth hists exist
    h_wrong = get_hist(f, "AFBb/h_abs_wrong_truth")
    h_tot   = get_hist(f, "AFBb/h_abs_tot_truth")
    if h_wrong and h_tot:
        gr = omega_graph_from_pair(h_wrong, h_tot)
        if gr and gr.GetN()>0:
            cW = ROOT.TCanvas("c_omega_abs", "omega_vs_abs", 800, 600)
            cW.SetGrid(1,1)
            gr.SetMarkerStyle(21); gr.SetMarkerSize(1.0); gr.SetLineWidth(2)
            gr.Draw("AP")
            save_canvas(cW, outdir, "AFBb_omega_vs_abs")

    # response matrix if present
    h_resp = get_hist(f, "AFBb/resp_truth_vs_signed")
    if h_resp:
        cR = ROOT.TCanvas("c_resp", "response_matrix", 800, 650)
        ROOT.gStyle.SetNumberContours(100)
        h_resp.GetXaxis().SetTitle("cos#theta_{truth}")
        h_resp.GetYaxis().SetTitle("cos#theta_{signed}")
        h_resp.Draw("COLZ")
        save_canvas(cR, outdir, "AFBb_response_matrix")

    # numbers panel only (compact)
    cN = ROOT.TCanvas("c_numbers","AFBb numbers", 700, 400)
    frame = ROOT.TH1F("frame","",1,0,1); frame.SetStats(0)
    frame.GetXaxis().SetLabelSize(0); frame.GetYaxis().SetLabelSize(0); frame.Draw()
    pave2 = make_numbers_pave(s)
    pave2.Draw()
    save_canvas(cN, outdir, "AFBb_numbers")

    f.Close()

def main():
    parser = argparse.ArgumentParser(description="AFB^b Stage-2 plotting (exclusive Z->bb)")
    parser.add_argument("--root", default="outputs/bafb/stage2/afbb_stage2.root")
    parser.add_argument("--json", default="outputs/bafb/stage2/results/AFBb_fit.json")
    parser.add_argument("--out",  default="outputs/bafb/stage2/plots")
    args = parser.parse_args()

    ensure_dir(args.out)
    draw_overlay(args.root, args.json, args.out)
    print(f"[plots] All figures saved under: {args.out}")

if __name__ == "__main__":
    main()
