# analysis_stage1.py
# Stage-1 of the b-quark A_FB analysis using FCCAnalyses.
# - Reclusters events into exactly two e+e- jets (Durham/ee-kT).
# - Computes per-jet charge with |p|^kappa weighting and chooses the b-jet
#   as the jet with the more negative charge (conventional rule).
# - Defines the A_FB observable as cos(theta) of the chosen b-jet w.r.t. +z (e- beam).
# - Stores an optional truth-level cos(theta_b) with a proper guard (sentinel = -2.f).
# - Adds diagnostics: DeltaQ between the two jets’ charges, jet charged multiplicities,
#   and a simple visible s'/s proxy with an optional cut to emulate “pole” conditions.
#
# Requires utils_bafb.py (jet/MC helpers) importable via:
#   export PYTHONPATH="$(pwd)/analysis/bafb/python:${PYTHONPATH:-}"

import os
import ROOT
from utils_bafb import *  # injects C++ helpers: bafb_jet_charge, bafb_jet_nch, bafb_choose_bjet, bafb_is_b, bafb_sum_mask, bafb_argmax_masked


output_dir = "outputs/bafb/stage1"


class Analysis:
    def __init__(self, args):
        # -------------------- configurable knobs --------------------
        self.kappa = float(getattr(args, "kappa", 0.5))        # jet-charge exponent
        self.sprime_cut = float(getattr(args, "sprime_cut", 0.0))  # visible s'/s cut; 0 disables
        self.nThreads = int(getattr(args, "n_threads", 50))     # respected by FCCAnalyses managed mode

        # Managed-mode boilerplate
        self.outputDir = output_dir
        self.doTree = True
        self.doSkim = False
        self.doPlot = False

        # In managed mode we rely on --input/--input-file-list or --test; keep these empty.
        self.processList = {}
        self.prodTag = ""

        # Tiny test file to validate chain when --test is used (framework handles this).
        self.testFile = "https://fccsw.web.cern.ch/fccsw/analysis/test-samples/edm4hep099/p8_ee_WW_ecm240_edm4hep.root"

        os.makedirs(self.outputDir, exist_ok=True)

        # Let FCCAnalyses/ROOT manage threads; just print what we see so logs are clear.
        try:
            act = ROOT.ROOT.GetThreadPoolSize() if hasattr(ROOT.ROOT, "GetThreadPoolSize") else -1
        except Exception:
            act = -1
        print(f"[stage1] (framework-managed threads)  kappa={self.kappa}  sPrimeCut={self.sprime_cut}  threads~{act}")

    def analyzers(self, df):
        """
        Build the RDataFrame pipeline:
          - ReconstructedParticles -> pseudojets -> ee_kt exclusive 2-jet clustering
          - Per-jet charge from constituents (|p|^kappa)
          - Pick b-jet index; define cosTheta as cos(theta_b)
          - Truth helpers via robust C++ helpers (no inline JIT lambdas)
          - Diagnostics: DeltaQ, jet charged multiplicity, s'/s proxy (optional cut)
        """
        d = (
            df
            # ---- EDM4hep aliases ----
            .Alias("RP", "ReconstructedParticles")
            .Alias("MC", "Particle")

            # ---- particle kinematics / charges ----
            .Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(RP)")
            .Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(RP)")
            .Define("RP_pz", "FCCAnalyses::ReconstructedParticle::get_pz(RP)")
            .Define("RP_e",  "FCCAnalyses::ReconstructedParticle::get_e(RP)")
            .Define("RP_p",  "FCCAnalyses::ReconstructedParticle::get_p(RP)")
            .Define("RP_q",  "FCCAnalyses::ReconstructedParticle::get_charge(RP)")

            # ---- recluster exclusive 2 jets in e+e- mode (Durham/ee-kT) ----
            # NOTE: This mirrors your current usage; keep as-is for compatibility.
            .Define("fj_pseudojets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
            .Define("fj_clusters",    "JetClustering::clustering_ee_kt(2, 2, 1, 0)(fj_pseudojets)")
            .Define("jets",           "FCCAnalyses::JetClusteringUtils::get_pseudoJets(fj_clusters)")
            .Define("jet_theta",      "FCCAnalyses::JetClusteringUtils::get_theta(jets)")
            .Define("jet_constits",   "FCCAnalyses::JetClusteringUtils::get_constituents(fj_clusters)")

            # ---- per-jet charge with |p|^kappa ----
            .Define("jet_Q",     f"bafb_jet_charge(RP_q, RP_p, jet_constits, {self.kappa})")
            .Define("bjet_idx",  "bafb_choose_bjet(jet_Q)")
            .Define("other_idx", "(bjet_idx==0)?1:0")

            # ΔQ = Q_b - Q_{bbar} (guard indices)
            .Define(
                "DeltaQ",
                " (bjet_idx>=0 && other_idx>=0 && bjet_idx<(int)jet_Q.size() && other_idx<(int)jet_Q.size()) "
                "   ? (jet_Q[bjet_idx] - jet_Q[other_idx]) : 0.f"
            )

            # jet charged multiplicity per jet (counts RP with non-zero charge among constituents)
            .Define("jet_nch", "bafb_jet_nch(RP_q, jet_constits)")

            # ---- AFB observable: cos(theta) of the chosen b-jet (beam +z is forward) ----
            .Define("cosTheta",  "(bjet_idx>=0 && bjet_idx<(int)jet_theta.size()) ? std::cos(jet_theta[bjet_idx]) : -2.f")
            .Define("isForward", "cosTheta>0.f")

            # ---- (optional) truth helpers for dilution/purity studies (robust C++ helpers) ----
            .Define("MC_pdg",   "FCCAnalyses::MCParticle::get_pdg(MC)")
            .Define("MC_theta", "FCCAnalyses::MCParticle::get_theta(MC)")
            .Define("MC_p",     "FCCAnalyses::MCParticle::get_p(MC)")
            .Define("is_b",     "bafb_is_b(MC_pdg)")                      # RVec<int> 0/1 mask
            .Define("n_b",      "bafb_sum_mask(is_b)")                    # int
            .Define("b_idx",    "bafb_argmax_masked(MC_p, is_b)")         # int (or -1)
            .Define("cth_truth","(n_b>0 && b_idx >= 0 && b_idx < (int)MC_theta.size()) ? std::cos(MC_theta[b_idx]) : -2.f")

            # ---- simple visible s'/s proxy from reconstructed four-momentum (optional ISR control) ----
            # s' ≈ m_vis^2; we form (m_vis^2 / s) with sqrt(s)=91.1876 GeV unless you later expose it as a knob.
            .Define("RP_px_sum","ROOT::VecOps::Sum(RP_px)")
            .Define("RP_py_sum","ROOT::VecOps::Sum(RP_py)")
            .Define("RP_pz_sum","ROOT::VecOps::Sum(RP_pz)")
            .Define("RP_e_sum", "ROOT::VecOps::Sum(RP_e)")
            .Define(
                "sprime_over_s",
                " (RP_e_sum>0) ? ((RP_e_sum*RP_e_sum - (RP_px_sum*RP_px_sum+RP_py_sum*RP_py_sum+RP_pz_sum*RP_pz_sum)) / (91.1876*91.1876)) : 0.0"
            )
        )

        # Optional selection to emulate “pole” conditions; disabled by default (sprime_cut=0).
        if self.sprime_cut > 0.0:
            d = d.Filter(f"sprime_over_s>={self.sprime_cut}", "s'/s pole selection")

        return d

    def output(self):
        """
        Columns to Snapshot into the stage-1 slim tree.
        The output filename is set by the CLI (--output).
        """
        return [
            # AFB observable
            "cosTheta", "isForward",
            # jet-level diagnostics
            "jet_Q", "jet_nch", "bjet_idx", "other_idx", "DeltaQ",
            # truth helper for dilution studies
            "cth_truth",
            # ISR proxy
            "sprime_over_s",
        ]
