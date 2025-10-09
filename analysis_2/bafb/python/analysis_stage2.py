# analysis_stage2.py
# Stage-2: build histograms and A_FB from the Stage-1 slim tree, with
#   - optional quality cuts on |DeltaQ| and jet charged multiplicity (n_ch of b-jet)
#   - inclusive & binned A_FB (raw), ω (wrong-sign), and guarded corrected A_FB
#   - truth↔reco response matrix
#   - optional linear-fit cross-check
#   - rich JSON summary with provenance and guard flags
#
# CLI knobs you can pass via `fccanalysis run ... --args ...`:
#   --n-threads 16
#   --nCosBins 20
#   --do-fit true|false
#   --minDeltaQ 0.0         # apply |DeltaQ| >= value if > 0
#   --minJetNch 0           # require nch_b >= value if > 0
#   --minScaleGuard 0.10    # if D = 1-2ω < minScaleGuard, skip corrected A_FB
#
# Notes:
#  - This stage expects Stage-1 to have stored: cosTheta, abscth (we recompute),
#    DeltaQ, jet_nch (RVec<int>), bjet_idx, cth_truth (optional).
#  - If any of those are missing, quality cuts are simply disabled and we continue.

import os, json, math, ROOT
from utils_bafb import afb_corrected

output_dir = "outputs/bafb/stage2"


def _boolify(x, default=True):
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    s = str(x).strip().lower()
    if s in ("1", "true", "yes", "on"):  return True
    if s in ("0", "false", "no", "off"): return False
    return default


class Analysis:
    def __init__(self, args):
        os.makedirs(output_dir, exist_ok=True)

        # knobs (read from CLI; safe defaults)
        self.n_threads     = int(getattr(args, "n_threads", 4))
        self.input_path    = getattr(args, "input", None)
        self.do_fit        = _boolify(getattr(args, "do_fit", True), True)
        self.nbins         = int(getattr(args, "nCosBins", 20))
        self.minDeltaQ     = float(getattr(args, "minDeltaQ", 0.0))
        self.minJetNch     = int(getattr(args, "minJetNch", 0))
        self.minScaleGuard = float(getattr(args, "minScaleGuard", 0.10))  # guard on D = 1-2ω

        self.cfg = {
            "n_threads": self.n_threads,
            "input": self.input_path,
            "binning": {"nCosBins": self.nbins, "range": [-1.0, 1.0]},
            "quality": {
                "minDeltaQ": self.minDeltaQ,
                "minJetNch": self.minJetNch,
                "minScaleGuard": self.minScaleGuard,
            },
            "note": "ω measured from truth if present; purity P kept at 1.0 (no b-tag yet). Guard skips A_FB^corr if D<minScaleGuard.",
        }
        print(f"[stage2] (framework-managed threads) input={self.input_path}  "
              f"nbins={self.nbins}  do_fit={self.do_fit}  "
              f"minDeltaQ={self.minDeltaQ}  minJetNch={self.minJetNch}  "
              f"minScaleGuard={self.minScaleGuard}")

    # -------------------------- analysis graph --------------------------- #
    def analyzers(self, df):
        # keep only valid signed cosθ; cache useful columns
        d0 = (
            df
            .Filter("cosTheta>-1.5f && cosTheta<1.5f", "remove sentinels")
            .Define("abscth", "fabs(cosTheta)")
            .Define("isF", "cosTheta>0.0")
            .Define("isB", "cosTheta<0.0")
        )

        # Optional quality cuts if Stage-1 columns exist: DeltaQ, jet_nch (vector), bjet_idx
        cols = set(map(str, df.GetColumnNames()))
        have_quality = {"DeltaQ", "jet_nch", "bjet_idx"}.issubset(cols)

        if have_quality and (self.minDeltaQ > 0.0 or self.minJetNch > 0):
            d = (
                d0
                .Define("nch_b", "(bjet_idx>=0 && bjet_idx<(int)jet_nch.size()) ? jet_nch[bjet_idx] : 0")
                .Define("passCuts",
                        f"(fabs(DeltaQ) >= {self.minDeltaQ}) && (nch_b >= {self.minJetNch})")
                .Filter("passCuts", "quality cuts: |DeltaQ| and jet n_ch")
            )
            self.nAll = d0.Count()
            self.nSel = d.Count()
        else:
            d = d0
            self.nAll = d.Count()
            self.nSel = self.nAll
            if not have_quality and (self.minDeltaQ > 0.0 or self.minJetNch > 0):
                print("[stage2][warn] Requested quality cuts but DeltaQ/jet_nch/bjet_idx not in tree; cuts disabled.")

        # Signed and |cosθ| distributions
        self.h_cos = d.Histo1D(("h_cos", "signed cos#theta;cos#theta;Events", self.nbins, -1, 1), "cosTheta")
        self.h_abs = d.Histo1D(("h_abs", "|cos#theta|;|cos#theta|;Events", self.nbins, 0, 1), "abscth")

        # Inclusive forward/backward counts
        self.nF = d.Sum("isF")
        self.nB = d.Sum("isB")

        # F/B vs |cosθ| for binned A_FB
        self.h_absF = d.Filter("isF").Histo1D(("h_absF", "Forward;|cos#theta|;Events", self.nbins, 0, 1), "abscth")
        self.h_absB = d.Filter("isB").Histo1D(("h_absB", "Backward;|cos#theta|;Events", self.nbins, 0, 1), "abscth")

        # Truth-based diagnostics (if available)
        self.h_cth_truth = None
        self.h_resp = None
        self.omega_mean = None
        self.h_abs_tot_t = None   # total truth-matched per |cosθ|
        self.h_abs_wrong = None   # wrong-sign (weight=1) per |cosθ|
        try:
            dt = d.Filter("cth_truth>-1.5f && cth_truth<1.5f", "require truth cosθ")
            # Truth spectrum
            self.h_cth_truth = dt.Histo1D(("h_cth_truth", "truth cos#theta_{b};cos#theta_{b};Events", self.nbins, -1, 1), "cth_truth")
            # Response matrix (truth on X, reco on Y)
            self.h_resp = dt.Histo2D(("h_resp", ";cos#theta_{truth};cos#theta_{reco}", self.nbins, -1, 1, self.nbins, -1, 1),
                                     "cth_truth", "cosTheta")
            # Inclusive ω (mean wrong flag)
            self.omega_mean = dt.Define("wrong", "(cth_truth>0)!=(cosTheta>0)").Mean("wrong")
            # ω(|cosθ|): build wrong (weight=1) and total per abscth bin
            dtw = dt.Define("w_wrong", "((cth_truth>0)!=(cosTheta>0)) ? 1.f : 0.f")
            self.h_abs_tot_t = dt.Histo1D(("h_abs_tot_t", "truth-matched total;|cos#theta|;Events", self.nbins, 0, 1), "abscth")
            self.h_abs_wrong = dtw.Histo1D(("h_abs_wrong", "truth-matched wrong;|cos#theta|;Wrong", self.nbins, 0, 1),
                                           "abscth", "w_wrong")
        except Exception:
            # No truth available → ω, response, and ω(|cosθ|) will be absent
            pass

        return d

    # ------------------------------ helpers ------------------------------ #
    def _binned_afb(self, hF, hB):
        """Compute per-|cosθ| A_FB and binomial errors from TH1 histos."""
        nb = hF.GetNbinsX()
        A, dA = [], []
        for i in range(1, nb + 1):
            nf  = hF.GetBinContent(i)
            nbk = hB.GetBinContent(i)
            n   = nf + nbk
            if n <= 0:
                A.append(0.0); dA.append(0.0); continue
            a = (nf - nbk) / n
            err = math.sqrt(max(0.0, (1.0 - a*a) / n))
            A.append(a); dA.append(err)
        edges = [hF.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)]
        edges.append(hF.GetXaxis().GetBinUpEdge(nb))
        return edges, A, dA

    def _omega_vs_abs(self, h_wrong, h_tot):
        """Per-|cosθ| ω = wrong / total with binomial errors."""
        nb = h_tot.GetNbinsX()
        om, dom = [], []
        for i in range(1, nb + 1):
            tot = h_tot.GetBinContent(i)
            wrg = h_wrong.GetBinContent(i)
            if tot <= 0:
                om.append(0.0); dom.append(0.0); continue
            w = wrg / tot
            err = math.sqrt(max(0.0, w * (1.0 - w) / tot))
            om.append(w); dom.append(err)
        edges = [h_tot.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)]
        edges.append(h_tot.GetXaxis().GetBinUpEdge(nb))
        return edges, om, dom

    # ------------------------------ outputs ------------------------------ #
    def output(self):
        # Event accounting
        n_all = int(self.nAll.GetValue())
        n_sel = int(self.nSel.GetValue())
        frac  = (float(n_sel) / n_all) if n_all > 0 else 0.0

        # Inclusive A_FB (raw)
        nf = int(self.nF.GetValue())
        nb = int(self.nB.GetValue())
        N  = nf + nb
        if N <= 0:
            A_raw, sig_raw = 0.0, 0.0
        else:
            A_raw  = (nf - nb) / N
            sig_raw = math.sqrt(max(0.0, (1.0 - A_raw * A_raw) / N))

        # ω (inclusive) and scale factor D = 1-2ω
        P = 1.0  # no b-tag purity yet
        omega = float(self.omega_mean.GetValue()) if self.omega_mean else None
        D = (1.0 - 2.0 * omega) if omega is not None else None

        # Guard: if D is too small (near-random tag), skip corrected A_FB to avoid exploding values
        corr_applied = False
        A_corr, sig_corr, corr_reason = None, None, None
        if D is None:
            corr_reason = "no_truth"
        else:
            if abs(D) < self.minScaleGuard:
                corr_reason = f"low_scale(D={D:.4f}<minScaleGuard={self.minScaleGuard:.3f})"
            else:
                A_corr, sig_corr = afb_corrected(A_raw, sig_raw, P, omega)
                corr_applied = True

        # Linear-fit cross-check (on signed distribution)
        A_fit, A_fit_err = None, None
        if self.do_fit:
            f = ROOT.TF1("afb", "[0]*(1+x*x+2*[1]*x)", -1, 1)  # [1] is A_FB
            self.h_cos.GetValue().Fit(f, "Q")                  # quiet fit
            A_fit, A_fit_err = float(f.GetParameter(1)), float(f.GetParError(1))

        # Binned A_FB(|cosθ|) (raw)
        hF, hB = self.h_absF.GetValue(), self.h_absB.GetValue()
        edges_absc, Afb_bins, dAfb_bins = self._binned_afb(hF, hB)

        # ω(|cosθ|) (if truth is present) and per-bin corrected A_FB with guard
        om_edges = om_bins = dom_bins = None
        Afb_bins_corr = dAfb_bins_corr = None
        if self.h_abs_wrong and self.h_abs_tot_t:
            om_edges, om_bins, dom_bins = self._omega_vs_abs(self.h_abs_wrong.GetValue(),
                                                             self.h_abs_tot_t.GetValue())
            # Per-bin corrected A_FB using D_i = 1-2ω_i
            Afb_bins_corr, dAfb_bins_corr = [], []
            for a, da, w in zip(Afb_bins, dAfb_bins, om_bins):
                D_i = 1.0 - 2.0 * w
                if abs(D_i) < self.minScaleGuard:
                    Afb_bins_corr.append(None)
                    dAfb_bins_corr.append(None)
                else:
                    a_c, da_c = afb_corrected(a, da, 1.0, w)
                    Afb_bins_corr.append(a_c)
                    dAfb_bins_corr.append(da_c)

        # ROOT outputs
        fout = ROOT.TFile(os.path.join(output_dir, "bafb_outputs.root"), "RECREATE")
        self.h_cos.Write(); self.h_abs.Write()
        self.h_absF.Write(); self.h_absB.Write()
        if self.h_cth_truth: self.h_cth_truth.Write()
        if self.h_resp:      self.h_resp.Write()
        if self.h_abs_wrong: self.h_abs_wrong.Write()
        if self.h_abs_tot_t: self.h_abs_tot_t.Write()
        fout.Close()

        # JSON summary
        summary = {
            "n_all": n_all, "n_sel": n_sel, "sel_fraction": frac,
            "N_F": nf, "N_B": nb, "N_tot": N,
            "A_FB_raw": A_raw, "A_FB_raw_err": sig_raw,
            "purity_assumed": P,
            "omega_measured": omega,
            "dilution_D": D,
            "corr_applied": corr_applied,
            "corr_reason": corr_reason,
            "A_FB_corr": A_corr, "A_FB_corr_err": sig_corr,
            "A_FB_fit": A_fit, "A_FB_fit_err": A_fit_err,
            "afb_vs_abs": {
                "edges": edges_absc,
                "A_raw": Afb_bins,
                "sigma_raw": dAfb_bins,
                "A_corr": Afb_bins_corr,
                "sigma_corr": dAfb_bins_corr,
            },
            "omega_vs_abs": {
                "edges": om_edges,
                "omega": om_bins,
                "sigma": dom_bins,
            },
        }
        summary.update(self.cfg)

        with open(os.path.join(output_dir, "summary.json"), "w") as jf:
            json.dump(summary, jf, indent=2)

        # (Return list unused in managed mode)
        return ["cosTheta", "abscth"]
