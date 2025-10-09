# analysis_plots.py
# Reads stage-2 outputs and produces:
#  - cosTheta and |cosTheta| histograms (PNG/PDF)
#  - A_FB numbers panel (raw + corrected, plus P and ω)
#  - A_FB(|cosθ|) and ω(|cosθ|) graphs from JSON arrays
#  - response matrix heatmap if available

import os, json, math, ROOT

input_dir = "outputs/bafb/stage2"   # where stage2 wrote ROOT + JSON
outdir    = "outputs/bafb/plots"

def _ensure_dir(p):
    if not os.path.isdir(p):
        os.makedirs(p, exist_ok=True)

def _save(c, base):
    png = os.path.join(outdir, base + ".png")
    pdf = os.path.join(outdir, base + ".pdf")
    c.SaveAs(png)
    c.SaveAs(pdf)

def _graph_from_hist(h):
    """Make a TGraphErrors from a TH1: x at bin center, ex=bin width/2, y=content, ey=sqrt(N)."""
    nb = h.GetNbinsX()
    g = ROOT.TGraphErrors(nb)
    for i in range(1, nb+1):
        x  = h.GetXaxis().GetBinCenter(i)
        ex = 0.5 * h.GetXaxis().GetBinWidth(i)
        y  = h.GetBinContent(i)
        ey = math.sqrt(max(0.0, y))
        g.SetPoint(i-1, x, y)
        g.SetPointError(i-1, ex, ey)
    return g

def _graph_from_binedges(edges, yvals, yerrs=None, title=";|cos#theta|;value"):
    """
    Build TGraphErrors from bin edges + y (+/- err).
    edges: length nb+1
    yvals: length nb
    """
    nb = len(edges) - 1
    g = ROOT.TGraphErrors(nb)
    for i in range(nb):
        xL, xR = edges[i], edges[i+1]
        x  = 0.5*(xL + xR)
        ex = 0.5*(xR - xL)
        y  = float(yvals[i]) if i < len(yvals) else 0.0
        ey = float(yerrs[i]) if (yerrs is not None and i < len(yerrs)) else 0.0
        g.SetPoint(i, x, y)
        g.SetPointError(i, ex, ey)
    g.SetTitle(title)
    return g

class Plot:
    def __init__(self, args=None):
        _ensure_dir(outdir)
        # a clean, readable style
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetNumberContours(100)

    def plots(self):
        # --- open ROOT outputs ---
        root_path = os.path.join(input_dir, "bafb_outputs.root")
        f = ROOT.TFile.Open(root_path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Missing or invalid {root_path}")

        # Base histograms
        h_cos = f.Get("h_cos")
        h_abs = f.Get("h_abs")
        if not h_cos or not h_abs:
            raise RuntimeError("Expected h_cos and h_abs in stage-2 ROOT output")

        # Optional extras
        h_absF      = f.Get("h_absF")
        h_absB      = f.Get("h_absB")
        h_resp      = f.Get("h_resp")       # 2D response (truth vs reco)
        h_cth_truth = f.Get("h_cth_truth")
        h_wrong     = f.Get("h_abs_wrong")
        h_tot_t     = f.Get("h_abs_tot_t")

        # --- hist plots ---
        c1 = ROOT.TCanvas("c1","cos",800,600)
        h_cos.SetLineWidth(2)
        h_cos.Draw("hist")
        _save(c1, "cosTheta")

        c2 = ROOT.TCanvas("c2","abscos",800,600)
        h_abs.SetLineWidth(2)
        h_abs.Draw("hist")
        _save(c2, "abs_cosTheta")

        if h_resp:
            cR = ROOT.TCanvas("cR","response",720,600)
            h_resp.GetXaxis().SetTitle("cos#theta_{truth}")
            h_resp.GetYaxis().SetTitle("cos#theta_{reco}")
            h_resp.Draw("COLZ")
            _save(cR, "response_matrix")

        # --- read JSON numbers ---
        jpath = os.path.join(input_dir, "summary.json")
        with open(jpath) as j:
            s = json.load(j)

        # Extract core numbers (with fallbacks)
        Nf = int(s.get("N_F", 0))
        Nb = int(s.get("N_B", 0))
        N  = int(s.get("N_tot", Nf + Nb))

        Araw   = float(s.get("A_FB_raw", 0.0))
        sAraw  = float(s.get("A_FB_raw_err", 0.0))

        P      = float(s.get("purity_assumed", 1.0))
        omega  = float(s.get("omega_measured", s.get("omega_assumed", 0.0)))

        Acorr  = s.get("A_FB_corr", None)
        sAcorr = s.get("A_FB_corr_err", None)

        Afit   = s.get("A_FB_fit", None)
        sAfit  = s.get("A_FB_fit_err", None)

        # --- numbers panel ---
        c3 = ROOT.TCanvas("c3","Afb",900,380)
        frame = ROOT.TH1F("frame","",1,0,1); frame.SetStats(0)
        frame.GetXaxis().SetLabelSize(0); frame.GetYaxis().SetLabelSize(0); frame.Draw()
        txt = ROOT.TLatex(); txt.SetNDC(True); txt.SetTextSize(0.05)

        y = 0.82
        txt.DrawLatex(0.08, y,      f"N_F = {Nf:,d},  N_B = {Nb:,d},  N = {N:,d}"); y -= 0.12
        txt.DrawLatex(0.08, y,      f"A_{{FB}}^{{raw}} = {Araw:.4f} #pm {sAraw:.4f}"); y -= 0.12

        if Acorr is not None:
            txt.DrawLatex(0.08, y,  f"A_{{FB}}^{{corr}} = {Acorr:.4f}" + (f" #pm {float(sAcorr):.4f}" if sAcorr is not None else ""))
            y -= 0.12

        txt.DrawLatex(0.08, y,      f"P = {P:.3f},   #omega = {omega:.3f}"); y -= 0.12

        if Afit is not None:
            txt.DrawLatex(0.08, y,  f"Fit: A_{{FB}} = {float(Afit):.4f}" + (f" #pm {float(sAfit):.4f}" if sAfit is not None else "")); y -= 0.12

        if "input" in s and s["input"]:
            txt.SetTextSize(0.04)
            txt.DrawLatex(0.08, 0.12, f"input: {os.path.basename(s['input'])}  (threads={s.get('n_threads','?')})")

        _save(c3, "afb_numbers")

        # --- A_FB(|cosθ|) from JSON arrays (if present) ---
        afb_vs_abs = s.get("afb_vs_abs", None)
        if afb_vs_abs and afb_vs_abs.get("edges") and afb_vs_abs.get("A"):
            edges = afb_vs_abs["edges"]
            Avals = afb_vs_abs["A"]
            Aerrs = afb_vs_abs.get("sigma", None)
            gA = _graph_from_binedges(edges, Avals, Aerrs, title=";|cos#theta|;A_{FB}(|cos#theta|)")
            cA = ROOT.TCanvas("cA","afb_vs_abs",800,600)
            cA.SetGrid(1,1)
            gA.SetMarkerStyle(20); gA.SetMarkerSize(1.0); gA.SetLineWidth(2)
            gA.Draw("AP")
            _save(cA, "afb_vs_abs")

        # --- ω(|cosθ|) from JSON arrays (if present) ---
        om_vs_abs = s.get("omega_vs_abs", None)
        if om_vs_abs and om_vs_abs.get("edges") and om_vs_abs.get("omega"):
            edges = om_vs_abs["edges"]
            ovals = om_vs_abs["omega"]
            oerrs = om_vs_abs.get("sigma", None)
            gW = _graph_from_binedges(edges, ovals, oerrs, title=";|cos#theta|;#omega(|cos#theta|)")
            cW = ROOT.TCanvas("cW","omega_vs_abs",800,600)
            cW.SetGrid(1,1)
            gW.SetMarkerStyle(21); gW.SetMarkerSize(1.0); gW.SetLineWidth(2)
            gW.Draw("AP")
            _save(cW, "omega_vs_abs")

        # Optional: draw F/B per |cosθ| if present in ROOT (useful for quick QA)
        if h_absF and h_absB:
            cFB = ROOT.TCanvas("cFB","FvsB",800,600)
            h_absF.SetLineColor(ROOT.kBlue+1); h_absF.SetLineWidth(2)
            h_absB.SetLineColor(ROOT.kRed+1);  h_absB.SetLineWidth(2)
            h_absF.SetTitle("|cos#theta|;|cos#theta|;Events")
            h_absF.Draw("hist")
            h_absB.Draw("hist same")
            leg = ROOT.TLegend(0.65,0.75,0.88,0.88)
            leg.AddEntry(h_absF,"Forward","l")
            leg.AddEntry(h_absB,"Backward","l")
            leg.Draw()
            _save(cFB, "forward_backward_abs")

        # done
        print(f"[plots] Wrote figures to: {outdir}")

if __name__ == "__main__":
    Plot().plots()
