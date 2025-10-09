# analysis_stage2.py
# Stage-2: build histograms and A_FB from the Stage-1 slim tree.
# - Filters sentinels, makes signed/|cosθ| histograms
# - Computes N_F / N_B, inclusive A_FB (raw) with binomial error
# - If truth exists, measures wrong-sign fraction ω and corrects A_FB
# - Fills a truth↔reco response matrix
# - Optional linear-fit cross-check to f(x) ∝ (1 + x^2) + 2 A x
# - Builds per-|cosθ| A_FB bins and ω(|cosθ|)
# - Writes ROOT (histos) + JSON (numbers, arrays, provenance)

import os, json, math, ROOT
from utils_bafb import afb_corrected

output_dir = "outputs/bafb/stage2"

def _boolify(x, default=True):
    """Robust bool coercion for CLI args."""
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    s = str(x).strip().lower()
    if s in ("1", "true", "yes", "on"):
        return True
    if s in ("0", "false", "no", "off"):
        return False
    return default

class Analysis:
    def __init__(self, args):
        os.makedirs(output_dir, exist_ok=True)

        # knobs (can be passed via --n-threads, --do-fit, --nCosBins)
        self.n_threads = int(getattr(args, "n_threads", 4))
        self.input_path = getattr(args, "input", None)
        self.do_fit    = _boolify(getattr(args, "do_fit", True), True)
        self.nbins     = int(getattr(args, "nCosBins", 20))

        # simple provenance (extend later with key4hep/git, etc.)
        self.cfg = {
            "n_threads": self.n_threads,
            "input": self.input_path,
            "binning": {"nCosBins": self.nbins, "range": [0.0, 1.0]},
            "note": "omega measured from truth if available; purity P kept at 1.0 (no b-tag yet).",
        }
        print(f"[stage2] (framework-managed threads), input={self.input_path}, do_fit={self.do_fit}, nbins={self.nbins}")

    def analyzers(self, df):
        # base df (drop sentinels; cache |cosθ| and F/B flags)
        d = (
            df
            .Filter("cosTheta>-1.5f && cosTheta<1.5f")
            .Define("abscth","fabs(cosTheta)")
            .Define("isF","cosTheta>0.0")
            .Define("isB","cosTheta<0.0")
        )

        # signed and absolute histograms
        self.h_cos = d.Histo1D(("h_cos","signed cos#theta;cos#theta;Events", self.nbins, -1, 1), "cosTheta")
        self.h_abs = d.Histo1D(("h_abs","|cos#theta|;|cos#theta|;Events", self.nbins, 0, 1), "abscth")

        # forward/backward counts (inclusive)
        self.nF = d.Sum("isF")
        self.nB = d.Sum("isB")

        # per-|cosθ| forward/backward spectra (for binned A_FB)
        self.h_absF = d.Filter("isF").Histo1D(("h_absF","F: |cos#theta|;|cos#theta|;Events", self.nbins, 0, 1), "abscth")
        self.h_absB = d.Filter("isB").Histo1D(("h_absB","B: |cos#theta|;|cos#theta|;Events", self.nbins, 0, 1), "abscth")

        # truth-based diagnostics (if available)
        self.h_cth_truth = None
        self.h_resp      = None
        self.omega_mean  = None
        self.h_abs_tot_t = None   # total per |cosθ| (with truth)
        self.h_abs_wrong = None   # wrong-sign per |cosθ| (with truth)
        try:
            dt = d.Filter("cth_truth>-1.5f && cth_truth<1.5f")
            self.h_cth_truth = dt.Histo1D(("h_cth_truth","truth cos#theta_{b};cos#theta_{b};Events", self.nbins, -1, 1), "cth_truth")
            self.h_resp      = dt.Histo2D(("h_resp",";cos#theta_{truth};cos#theta_{reco}", self.nbins, -1, 1, self.nbins, -1, 1),
                                          "cth_truth", "cosTheta")
            # inclusive ω (mean of 0/1 wrong flag)
            self.omega_mean  = dt.Define("wrong","(cth_truth>0)!=(cosTheta>0)").Mean("wrong")
            # ω(|cosθ|): need wrong-count and total-count per |cosθ| bin
            dtw = dt.Define("w_wrong","((cth_truth>0)!=(cosTheta>0)) ? 1.f : 0.f")
            self.h_abs_tot_t = dt.Histo1D(("h_abs_tot_t","truth-matched total;|cos#theta|;Events", self.nbins, 0, 1), "abscth")
            self.h_abs_wrong = dtw.Histo1D(("h_abs_wrong","truth-matched wrong;|cos#theta|;Wrong", self.nbins, 0, 1), "abscth", "w_wrong")
        except Exception:
            pass

        return d

    def _binned_afb(self, hF, hB):
        """Compute per-bin A_FB and errors from TH1F F/B spectra."""
        F = hF.GetArray()
        B = hB.GetArray()
        nb = hF.GetNbinsX()
        A, dA = [], []
        for i in range(1, nb+1):
            nf = hF.GetBinContent(i)
            nbk = hB.GetBinContent(i)
            n = nf + nbk
            if n <= 0:
                A.append(0.0); dA.append(0.0); continue
            a = (nf - nbk) / n
            err = math.sqrt(max(0.0, (1.0 - a*a)/n))
            A.append(a); dA.append(err)
        # bin edges
        edges = [hF.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)]
        edges.append(hF.GetXaxis().GetBinUpEdge(nb))
        return edges, A, dA

    def _omega_vs_abs(self, h_wrong, h_tot):
        """Per-|cosθ| omega = wrong / total, with simple binomial error."""
        nb = h_tot.GetNbinsX()
        om, dom = [], []
        for i in range(1, nb+1):
            tot = h_tot.GetBinContent(i)
            wrg = h_wrong.GetBinContent(i)
            if tot <= 0:
                om.append(0.0); dom.append(0.0); continue
            w = wrg / tot
            err = math.sqrt(max(0.0, w*(1.0 - w)/tot))
            om.append(w); dom.append(err)
        edges = [h_tot.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)]
        edges.append(h_tot.GetXaxis().GetBinUpEdge(nb))
        return edges, om, dom

    def output(self):
        # inclusive counts
        nf, nb = int(self.nF.GetValue()), int(self.nB.GetValue())
        N = nf + nb
        if N <= 0:
            A_raw, sig_raw = 0.0, 0.0
        else:
            A_raw = (nf - nb) / N
            sig_raw = math.sqrt(max(0.0, (1.0 - A_raw*A_raw) / N))

        # corrections (purity placeholder, ω from truth if available)
        P = 1.0
        omega = float(self.omega_mean.GetValue()) if self.omega_mean else 0.0
        A_corr, sig_corr = afb_corrected(A_raw, sig_raw, P, omega)

        # linear-fit cross-check (optional)
        A_fit, A_fit_err = None, None
        if self.do_fit:
            f = ROOT.TF1("afb","[0]*(1+x*x+2*[1]*x)", -1, 1)  # [1] = A_FB
            self.h_cos.GetValue().Fit(f, "Q")
            A_fit, A_fit_err = float(f.GetParameter(1)), float(f.GetParError(1))

        # per-|cosθ| A_FB bins
        hF = self.h_absF.GetValue()
        hB = self.h_absB.GetValue()
        edges_absc, Afb_bins, dAfb_bins = self._binned_afb(hF, hB)

        # ω(|cosθ|) bins (requires truth)
        om_edges, om_bins, dom_bins = None, None, None
        if self.h_abs_wrong and self.h_abs_tot_t:
            om_edges, om_bins, dom_bins = self._omega_vs_abs(self.h_abs_wrong.GetValue(),
                                                             self.h_abs_tot_t.GetValue())

        # ROOT outputs
        fout = ROOT.TFile(os.path.join(output_dir, "bafb_outputs.root"), "RECREATE")
        self.h_cos.Write(); self.h_abs.Write()
        self.h_absF.Write(); self.h_absB.Write()
        if self.h_cth_truth: self.h_cth_truth.Write()
        if self.h_resp: self.h_resp.Write()
        if self.h_abs_wrong: self.h_abs_wrong.Write()
        if self.h_abs_tot_t: self.h_abs_tot_t.Write()
        fout.Close()

        # JSON summary (numbers + arrays + provenance)
        summary = {
            "N_F": nf, "N_B": nb, "N_tot": N,
            "A_FB_raw": A_raw, "A_FB_raw_err": sig_raw,
            "purity_assumed": P,
            "omega_measured": omega,
            "A_FB_corr": A_corr, "A_FB_corr_err": sig_corr,
            "A_FB_fit": A_fit, "A_FB_fit_err": A_fit_err,
            "afb_vs_abs": {
                "edges": edges_absc,
                "A": Afb_bins,
                "sigma": dAfb_bins
            },
            "omega_vs_abs": {
                "edges": om_edges,
                "omega": om_bins,
                "sigma": dom_bins
            }
        }
        summary.update(self.cfg)

        with open(os.path.join(output_dir, "summary.json"), "w") as jf:
            json.dump(summary, jf, indent=2)

        # return list is ignored in managed mode; keep minimal
        return ["cosTheta", "abscth"]
