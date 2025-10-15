# analysis_stage2.py
# Stage-2 for A_FB^b with an exclusive Z->bb sample (no backgrounds).
# Consumes Stage-1 slim tree columns (notably: signed_cos, Q_primary, cth_truth, sprime_over_s, cfg_* echoes)
# Builds a signed cos(theta) histogram, measures (optionally) omega from truth, and fits the normalized
# template with a binned Poisson likelihood:
#
#   f(x | A_true, ω) ∝ (1 + x^2) + (8/3) * (1 - 2ω) * A_true * x,  x ∈ [-1,1]
#
# Options:
#   - fixed ω (default from truth if available else CLI prior; else 0.20)
#   - profile ω with a Gaussian prior (omega_sigma > 0)
#   - acceptance folding ε(x) from JSON (optional)
#
# Outputs:
#   - ROOT file with data histos and a TH1 overlay for the best-fit template (under AFBb/ TDirectory)
#   - PNG/PDF plot with the overlay
#   - JSON + TXT summary with A_FB_true, A_FB_meas, ω, errors, and provenance

import os
import json
import math
import ROOT

# Headless plotting (safer on lxplus / batch)
ROOT.gROOT.SetBatch(True)

# I/O locations
output_dir = "outputs/bafb/stage2"
results_dir = os.path.join(output_dir, "results")
plots_dir = os.path.join(output_dir, "plots")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)


def _boolify(x, default=False):
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


def _edges(nbins):
    return [-1.0 + i * (2.0 / nbins) for i in range(nbins + 1)]


def _I(a, b):
    # ∫(1 + x^2) dx
    return (b - a) + (b ** 3 - a ** 3) / 3.0


def _J(a, b):
    # ∫ x dx
    return (b ** 2 - a ** 2) / 2.0


def _normalized_bin_fracs(edges, Atrue, omega, eff=None):
    """Return normalized expected fractions per bin with optional acceptance (eff per bin)."""
    C = (8.0 / 3.0) * (1.0 - 2.0 * omega) * Atrue
    # total integral (odd term vanishes over symmetric [-1,1])
    tot = _I(-1.0, 1.0) + C * _J(-1.0, 1.0)

    vals = []
    for i in range(len(edges) - 1):
        a, b = edges[i], edges[i + 1]
        v = _I(a, b) + C * _J(a, b)
        if eff is not None:
            v *= max(0.0, eff[i])
        vals.append(v)

    if eff is not None:
        tot = sum(vals)

    # normalize
    if tot <= 0:
        return [1.0 / (len(edges) - 1)] * (len(edges) - 1)
    return [v / tot for v in vals]


def _nll_poisson(counts, mu):
    # -2 ln L up to a constant
    s = 0.0
    for n, m in zip(counts, mu):
        if m <= 0.0:
            return float("inf")
        if n > 0.0:
            s += 2.0 * (m - n + n * math.log(n / m))
        else:
            s += 2.0 * m
    return s


def _fit_binned_poisson(counts, edges, omega0, omega_sigma=0.0, eff=None):
    """Grid-scan fit over A_true (and ω if profiled) returning best-fit and error on A."""
    N = float(sum(counts))
    if N <= 0:
        return {"Atrue": 0.0, "Aerr": float("nan"), "omega": omega0, "nll": float("nan")}

    # scan ranges
    A_scan = [-0.30 + i * 0.0005 for i in range(int((0.30 - (-0.30)) / 0.0005) + 1)]
    if omega_sigma > 0.0:
        om_lo, om_hi = max(0.0, omega0 - 0.25), min(0.49, omega0 + 0.25)
        O_scan = [om_lo + i * 0.002 for i in range(int((om_hi - om_lo) / 0.002) + 1)]
    else:
        O_scan = [omega0]

    best = {"Atrue": 0.0, "omega": omega0, "nll": float("inf")}
    for A in A_scan:
        for om in O_scan:
            fr = _normalized_bin_fracs(edges, A, om, eff=eff)
            mu = [N * f for f in fr]
            nll = _nll_poisson(counts, mu)
            if omega_sigma > 0.0:
                nll += ((om - omega0) / omega_sigma) ** 2
            if nll < best["nll"]:
                best = {"Atrue": A, "omega": om, "nll": nll}

    # curvature on A at best (fix ω at best)
    h = 1e-3

    def nllA(Av):
        fr = _normalized_bin_fracs(edges, Av, best["omega"], eff=eff)
        mu = [N * f for f in fr]
        return _nll_poisson(counts, mu)

    n0 = nllA(best["Atrue"])
    np = nllA(best["Atrue"] + h)
    nm = nllA(best["Atrue"] - h)
    curv = (np - 2.0 * n0 + nm) / (h * h)
    Aerr = math.sqrt(2.0 / curv) if curv > 0 else float("nan")
    return {"Atrue": best["Atrue"], "Aerr": Aerr, "omega": best["omega"], "nll": best["nll"]}


def _load_acceptance_json(path, nbins):
    try:
        with open(path, "r") as f:
            data = json.load(f)
        edges = data.get("edges", None)
        eff = data.get("eff", None)
        if not isinstance(edges, list) or not isinstance(eff, list):
            print(f"[stage2] acceptance JSON missing 'edges'/'eff' lists -> ignore")
            return None, None
        if len(edges) != nbins + 1 or len(eff) != nbins:
            print(f"[stage2] acceptance JSON bins mismatch (nbins={nbins}) -> ignore")
            return None, None
        eff = [max(0.0, float(x)) for x in eff]
        return edges, eff
    except Exception as e:
        print(f"[stage2] could not load acceptance JSON: {e}")
        return None, None


class Analysis:
    def __init__(self, args):
        # CLI knobs
        self.n_threads = int(getattr(args, "n_threads", 4))
        self.tree_name = getattr(args, "tree", "events")
        self.nbins = int(getattr(args, "nbins", 20))
        self.omega_prior = getattr(args, "omega", None)  # None -> from truth or default 0.20
        self.omega_sigma = float(getattr(args, "omega_sigma", 0.0))
        self.fit_unbinned = _boolify(getattr(args, "fit_unbinned", False), False)
        self.fold_acceptance = getattr(args, "fold_acceptance", "")
        self.sprime_min = float(getattr(args, "sprime_min", -1.0))

        # provenance bag; will be extended in analyzers() once cfg_* are read
        self.prov = {
            "nbins": self.nbins,
            "omega_prior": self.omega_prior,
            "omega_sigma": self.omega_sigma,
            "fit_unbinned": self.fit_unbinned,
            "fold_acceptance": self.fold_acceptance if self.fold_acceptance else None,
            "sprime_min": self.sprime_min if self.sprime_min > 0 else None,
        }
        print(
            f"[stage2] cfg: nbins={self.nbins}  omega_prior={self.omega_prior}  "
            f"omega_sigma={self.omega_sigma}  unbinned={self.fit_unbinned}  "
            f"fold_acceptance={bool(self.fold_acceptance)}  sprime_min={self.sprime_min}"
        )

        # placeholders for RResultPtr histos and counters
        self.h_signed = None
        self.h_Qpri = self.h_Q0 = self.h_Q1 = None
        self.h_cos0 = self.h_cos1 = None
        self.h_sprime = None
        self.h_cth_truth = None
        self.h_resp = None
        self.h_abs_wrong = self.h_abs_tot = None
        self.nF = self.nB = None
        self.omega_mean = None
        self.cfg_echo_means = {}

    def _maybe_mean(self, df, col):
        try:
            m = df.Mean(col)
            self.cfg_echo_means[col] = m
            return True
        except Exception:
            return False

    def analyzers(self, df):
        # Filter out sentinels; define abs and F/B; apply sprime cut if requested
        d = (
            df.Filter("signed_cos>-1.5f && signed_cos<1.5f", "valid signed_cos")
            .Define("absx", "fabs(signed_cos)")
            .Define("isF", "signed_cos>0.0")
            .Define("isB", "signed_cos<0.0")
        )
        if self.sprime_min > 0:
            d = d.Filter(f"sprime_over_s >= {self.sprime_min}", f"s'/s >= {self.sprime_min}")

        # Histograms (IMPORTANT: no '/' in object names)
        self.h_signed = d.Histo1D(("h_signed_cosTheta", "signed cos#theta;cos#theta;Events", self.nbins, -1, 1), "signed_cos")
        self.h_Qpri = d.Histo1D(("h_Q_primary", "Q_{primary};Q;Events", 60, -1, 1), "Q_primary")
        try:
            self.h_Q0 = d.Histo1D(("h_Q0", "Q_{0};Q;Events", 60, -1, 1), "Q0")
            self.h_Q1 = d.Histo1D(("h_Q1", "Q_{1};Q;Events", 60, -1, 1), "Q1")
            self.h_cos0 = d.Histo1D(("h_cos0", "cos#theta_{0};cos#theta;Events", 40, -1, 1), "cos0")
            self.h_cos1 = d.Histo1D(("h_cos1", "cos#theta_{1};cos#theta;Events", 40, -1, 1), "cos1")
        except Exception:
            pass
        try:
            self.h_sprime = d.Histo1D(("h_sprime", "s'/s; s'/s; Events", 50, 0, 1.2), "sprime_over_s")
        except Exception:
            pass

        # counts
        self.nF = d.Sum("isF")
        self.nB = d.Sum("isB")

        # Truth-based ω (if truth exists)
        try:
            dt = d.Filter("cth_truth>-1.5f && cth_truth<1.5f", "valid truth cth")
            self.h_cth_truth = dt.Histo1D(("h_cth_truth", "truth cos#theta_b;cos#theta_b;Events", self.nbins, -1, 1), "cth_truth")
            self.h_resp = dt.Histo2D(("resp_truth_vs_signed", ";cos#theta_{truth};cos#theta_{signed}", self.nbins, -1, 1, self.nbins, -1, 1), "cth_truth", "signed_cos")
            # inclusive ω
            self.omega_mean = dt.Define("wrong", "(cth_truth>0)!=(signed_cos>0)").Mean("wrong")
            # ω(|x|)
            dtw = dt.Define("w_wrong", "((cth_truth>0)!=(signed_cos>0)) ? 1.f : 0.f")
            self.h_abs_tot = dt.Histo1D(("h_abs_tot_truth", "|x| total (truth-matched);|x|;Events", self.nbins, 0, 1), "absx")
            self.h_abs_wrong = dtw.Histo1D(("h_abs_wrong_truth", "|x| wrong (truth-matched);|x|;Wrong", self.nbins, 0, 1), "absx", "w_wrong")
        except Exception:
            pass

        # Capture cfg_* echoes if present
        for c in ("cfg_kappa", "cfg_track_p_min", "cfg_abs_costheta_max", "cfg_jet_p_min", "cfg_opposite_sign", "cfg_sprime_cut"):
            self._maybe_mean(d, c)

        return d

    def _afb_counting(self, nf, nb):
        N = nf + nb
        if N <= 0:
            return 0.0, 0.0
        A = (nf - nb) / N
        sig = math.sqrt(max(0.0, (1.0 - A * A) / N))
        return A, sig

    def _omega_vs_abs(self, h_wrong, h_tot):
        nb = h_tot.GetNbinsX()
        E, W, dW = [], [], []
        for i in range(1, nb + 1):
            tot = h_tot.GetBinContent(i)
            wrg = h_wrong.GetBinContent(i)
            eL = h_tot.GetXaxis().GetBinLowEdge(i)
            E.append(eL)  # we'll append the last upper after loop
            if tot <= 0:
                W.append(0.0)
                dW.append(0.0)
            else:
                w = wrg / tot
                err = math.sqrt(max(0.0, w * (1.0 - w) / tot))
                W.append(w)
                dW.append(err)
        E.append(h_tot.GetXaxis().GetXmax())
        return E, W, dW

    def _build_fit_overlay(self, h_data, edges, Atrue, omega, eff=None):
        """Return a TH1 with the best-fit template (scaled to N)."""
        N = h_data.Integral()
        nb = len(edges) - 1
        # IMPORTANT: never use '/' in object names
        h_fit = ROOT.TH1F("h_fit_overlay_tmp", "Best-fit template;cos#theta;Events", nb, -1, 1)
        fr = _normalized_bin_fracs(edges, Atrue, omega, eff=eff)
        for i in range(nb):
            h_fit.SetBinContent(i + 1, N * fr[i])
        h_fit.SetLineColor(ROOT.kRed + 1)
        h_fit.SetLineWidth(3)
        h_fit.SetFillStyle(0)
        return h_fit

    def output(self):
        # fetch histos and counts
        h_signed = self.h_signed.GetValue()
        nf = int(self.nF.GetValue())
        nb = int(self.nB.GetValue())
        N = nf + nb

        # counting A_FB (cross-check)
        A_count, A_count_err = self._afb_counting(nf, nb)

        # omega from truth if available
        omega_meas = self.omega_mean.GetValue() if self.omega_mean else None
        omega_prior = self.omega_prior if self.omega_prior is not None else (omega_meas if omega_meas is not None else 0.20)
        omega_sigma = self.omega_sigma

        # acceptance folding (optional)
        eff_edges = eff_vals = None
        if self.fold_acceptance:
            eff_edges, eff_vals = _load_acceptance_json(self.fold_acceptance, self.nbins)

        # Build counts and edges for fit
        edges = _edges(self.nbins)
        counts = [h_signed.GetBinContent(i + 1) for i in range(self.nbins)]

        # Fit (binned Poisson) - profile ω if sigma>0
        fitres = _fit_binned_poisson(counts, edges, omega_prior, omega_sigma, eff=eff_vals)
        Atrue = fitres["Atrue"]
        Aerr = fitres["Aerr"]
        omega_used = fitres["omega"]

        A_meas = (1.0 - 2.0 * omega_used) * Atrue

        # ---- Save ROOT artifacts (data + overlay + monitors) ----
        fout = ROOT.TFile(os.path.join(output_dir, "afbb_stage2.root"), "RECREATE")
        dir_afb = fout.mkdir("AFBb")
        dir_afb.cd()

        # Data histogram: clone to detach from RDF ownership; use clean name
        h_data = h_signed.Clone("h_signed_cosTheta")
        h_data.SetTitle("signed cos#theta;cos#theta;Events")
        h_data.SetDirectory(0)
        h_data.Write()
        h_data.SetDirectory(0)


        # Fit overlay (build -> rename -> write)
        h_fit = self._build_fit_overlay(h_signed, edges, Atrue, omega_used, eff=eff_vals)
        if not h_fit:
            raise RuntimeError("Internal error: fit overlay histogram is null")
        h_fit.SetName("h_fit_overlay")
        h_fit.SetTitle("Best-fit template;cos#theta;Events")
        h_fit.SetDirectory(0)
        h_fit.Write()
        h_fit.SetDirectory(0)

        # Optional monitors (write under same dir; ensure no '/' in names)
        for h in (
            self.h_Qpri,
            self.h_Q0,
            self.h_Q1,
            self.h_cos0,
            self.h_cos1,
            self.h_sprime,
            self.h_cth_truth,
            self.h_resp,
            self.h_abs_wrong,
            self.h_abs_tot,
        ):
            if h is None:
                continue
            try:
                obj = h.GetValue()
                if obj:
                    if "/" in obj.GetName():
                        obj.SetName(obj.GetName().replace("/", "_"))
                    obj.Write()
            except Exception:
                pass

        fout.Write()
        fout.Close()

        # ---- Make overlay plot (safe guards) ----
        if not h_data or not h_fit:
            raise RuntimeError("Null histogram(s) before plotting; aborting canvas draw")

        # Optional debug prints (comment out if noisy)
        print("[stage2/dbg] h_data ok:", bool(h_data), "entries:", (h_data.GetEntries() if h_data else -1))
        print("[stage2/dbg] h_fit  ok:", bool(h_fit), "integral:", (h_fit.Integral() if h_fit else -1))

        c = ROOT.TCanvas("c_afbb", "A_FB^b fit", 800, 600)
        h_data.SetLineWidth(2)
        h_data.SetMarkerStyle(20)
        h_data.Draw("E1")
        h_fit.SetLineColor(ROOT.kRed + 1)
        h_fit.SetLineWidth(3)
        h_fit.SetFillStyle(0)
        h_fit.Draw("HIST SAME")

        leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
        leg.AddEntry(h_data, "Data", "lep")
        leg.AddEntry(h_fit, "Best-fit template", "l")
        leg.Draw()

        pave = ROOT.TPaveText(0.15, 0.70, 0.52, 0.88, "NDC")
        pave.SetFillColor(0)
        pave.SetBorderSize(1)
        # Use Unicode ± (U+00B1) to avoid mojibake
        pave.AddText(f"A_FB^true = {Atrue:.5f} \u00B1 {Aerr:.5f}")
        pave.AddText(f"omega = {omega_used:.4f}")
        pave.AddText(f"N = {N}    bins = {self.nbins}")
        pave.Draw()

        c.SaveAs(os.path.join(plots_dir, "AFBb_signed_fit.png"))
        c.SaveAs(os.path.join(plots_dir, "AFBb_signed_fit.pdf"))
        del c

        # Compose provenance (include cfg_* echoes if present)
        prov = dict(self.prov)
        for k, m in self.cfg_echo_means.items():
            try:
                prov[k] = float(m.GetValue())
            except Exception:
                pass

        # JSON + TXT summary
        summary = {
            "N": N,
            "nbins": self.nbins,
            "A_FB_true": Atrue,
            "A_FB_true_err": Aerr,
            "A_FB_meas": A_meas,
            "omega_used": omega_used,
            "omega_prior": omega_prior,
            "omega_prior_sigma": omega_sigma,
            "counting_crosscheck": {"A_FB": A_count, "sigma": A_count_err},
            "acceptance_folded": bool(eff_vals is not None),
            "provenance": prov,
        }
        with open(os.path.join(results_dir, "AFBb_fit.json"), "w") as fjs:
            json.dump(summary, fjs, indent=2)
        with open(os.path.join(results_dir, "AFBb_fit.txt"), "w") as ftx:
            ftx.write(
                f"N = {N}\n"
                f"A_FB^true = {Atrue:.6f} +/- {Aerr:.6f}\n"
                f"omega_used = {omega_used:.4f}  (prior={omega_prior}, sigma={omega_sigma})\n"
                f"A_FB^meas  = {A_meas:.6f}\n"
            )

        print("====================================================")
        print(" AFBb Stage-2 fit (signed cos(theta))")
        print(f" N                : {N}")
        print(f" omega (used)     : {omega_used:.4f}  (prior={omega_prior}, sigma={omega_sigma})")
        print(f" A_FB^true        : {Atrue:.6f} +/- {Aerr:.6f}")
        print(f" A_FB^meas        : {A_meas:.6f}")
        print(f" ROOT out         : {os.path.join(output_dir, 'afbb_stage2.root')}")
        print(f" Results JSON/TXT : {os.path.join(results_dir, 'AFBb_fit.json')}")
        print(f" Plots            : {os.path.join(plots_dir, 'AFBb_signed_fit.[png|pdf]')}")
        print("====================================================")

        # managed-mode expects a list; it's not used downstream
        return ["signed_cos"]
