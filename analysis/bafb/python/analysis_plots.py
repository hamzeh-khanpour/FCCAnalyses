# analysis_plots.py
# Reads stage-2 outputs and produces:
#  - cosTheta and |cosTheta| histograms (PNG)
#  - A_FB numbers panel, including raw/corrected and (P, ω)

import os, json, ROOT

input_dir = "outputs/bafb/stage2"   # stage2 wrote here
outdir    = "outputs/bafb/plots"

class Plot:
    def __init__(self, args=None):
        os.makedirs(outdir, exist_ok=True)

    def plots(self):
        # ROOT output
        root_path = os.path.join(input_dir, "bafb_outputs.root")
        f = ROOT.TFile.Open(root_path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Missing or invalid {root_path}")

        h_cos = f.Get("h_cos")
        h_abs = f.Get("h_abs")
        if not h_cos or not h_abs:
            raise RuntimeError("Expected histograms h_cos and h_abs not found in stage-2 ROOT output")

        c1 = ROOT.TCanvas("c1","cos",800,600); h_cos.Draw("hist"); c1.SaveAs(os.path.join(outdir,"cosTheta.png"))
        c2 = ROOT.TCanvas("c2","abscos",800,600); h_abs.Draw("hist"); c2.SaveAs(os.path.join(outdir,"abs_cosTheta.png"))

        # Numbers panel from JSON
        jpath = os.path.join(input_dir,"summary.json")
        with open(jpath) as j:
            s = json.load(j)

        txt = ROOT.TLatex(); txt.SetNDC(True)
        c3 = ROOT.TCanvas("c3","Afb",800,320)
        frame = ROOT.TH1F("frame","",1,0,1); frame.SetStats(0)
        frame.GetXaxis().SetLabelSize(0); frame.GetYaxis().SetLabelSize(0); frame.Draw()

        txt.DrawLatex(0.10,0.80, f"N_F = {s['N_F']:,d}, N_B = {s['N_B']:,d}, N = {s.get('N_tot', s['N_F']+s['N_B']):,d}")
        txt.DrawLatex(0.10,0.65, f"A_FB^{{raw}} = {s['A_FB_raw']:.4f} #pm {s['A_FB_raw_err']:.4f}")

        if s.get("A_FB_corr") is not None:
            txt.DrawLatex(0.10,0.50, f"A_FB^{{corr}} = {s['A_FB_corr']:.4f} #pm {s.get('A_FB_corr_err',0.0):.4f}")

        if "purity_assumed" in s and "omega_assumed" in s:
            txt.DrawLatex(0.10,0.35, f"P = {s['purity_assumed']:.3f},  #omega = {s['omega_assumed']:.3f}")

        if "input" in s and s["input"]:
            txt.DrawLatex(0.10,0.20, f"input: {os.path.basename(s['input'])}  (threads={s.get('n_threads','?')})")

        c3.SaveAs(os.path.join(outdir,"afb_numbers.png"))

if __name__ == "__main__":
    Plot().plots()
