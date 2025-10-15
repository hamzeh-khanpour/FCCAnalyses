# analysis_stage1.py
# Stage-1 of the b-quark A_FB analysis (exclusive Z->bb, no backgrounds) using FCCAnalyses.
# Implements:
#  - Durham/ee-kT reclustering with N=2
#  - Charged-only JC3 jet charge via utils_bafb (kappa configurable; track p_min threshold)
#  - Deterministic primary jet = higher energy
#  - Signed polar angle: signed_cos = sign(Q_primary) * cos(theta_primary)
#  - Fiducial cuts on BOTH jets: |cos(theta)| < abs_costheta_max, |p| > jet_p_min
#  - Optional opposite-sign requirement on jet charges
#  - s'/s proxy with optional lower cut
#  - Truth hooks for closure (cth_truth)
#
# Conventions:
#  - +z is the electron-beam direction.
#  - Lab ~ CM at the Z pole; ISR can be controlled via s'/s if desired.
#
# Requires utils_bafb.py providing: bafb_jet_charge (5-arg overload), bafb_jet_nch, bafb_opposite_sign,
# bafb_choose_bjet (kept for diagnostics), bafb_is_b, bafb_sum_mask, bafb_argmax_masked.

import os
import ROOT
from utils_bafb import *  # C++ helpers are injected there

output_dir = "outputs/bafb/stage1"


class Analysis:
    def __init__(self, args):
        # -------------------- configurable knobs --------------------
        self.kappa = float(getattr(args, "kappa", 0.3))                 # JC3 exponent (default 0.3)
        self.track_p_min = float(getattr(args, "track_p_min", 0.5))     # GeV, threshold on track |p| in JC3
        self.abs_costheta_max = float(getattr(args, "abs_costheta_max", 0.97))
        self.jet_p_min = float(getattr(args, "jet_p_min", 20.0))        # GeV, per-jet |p|
        self.opposite_sign = bool(getattr(args, "opposite_sign", False))
        self.sprime_cut = float(getattr(args, "sprime_cut", 0.0))       # visible s'/s cut; 0 disables
        self.nThreads = int(getattr(args, "n_threads", 4))              # respected by FCCAnalyses managed mode

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

        # Print provenance summary (the framework controls threads)
        try:
            act = ROOT.ROOT.GetThreadPoolSize() if hasattr(ROOT.ROOT, "GetThreadPoolSize") else -1
        except Exception:
            act = -1
        print("[stage1] cfg: "
              f"kappa={self.kappa}  track_p_min={self.track_p_min} GeV  "
              f"|cosθ|<{self.abs_costheta_max}  p_jet>{self.jet_p_min} GeV  "
              f"opp_sign={int(self.opposite_sign)}  s'/s>={self.sprime_cut}  threads~{act}")

    def analyzers(self, df):
        """
        RDataFrame pipeline:
          - ReconstructedParticles -> pseudojets -> ee_kT exclusive 2-jet clustering
          - Per-jet charge (charged-only JC3, kappa, track_p_min)
          - Fiducials on BOTH jets
          - Optional opposite-sign requirement
          - Primary jet (higher energy)
          - Signed cos(theta)
          - Truth hooks & s'/s proxy
        """
        # make local copies for f-strings
        kappa = self.kappa
        pmin = self.track_p_min
        cth_max = self.abs_costheta_max
        pjet_min = self.jet_p_min

        d = (
            df
            # ---- EDM4hep aliases ----
            .Alias("RP", "ReconstructedParticles")
            .Alias("MC", "Particle")

            # ---- RP kinematics / charges ----
            .Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(RP)")
            .Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(RP)")
            .Define("RP_pz", "FCCAnalyses::ReconstructedParticle::get_pz(RP)")
            .Define("RP_e",  "FCCAnalyses::ReconstructedParticle::get_e(RP)")
            .Define("RP_p",  "FCCAnalyses::ReconstructedParticle::get_p(RP)")
            .Define("RP_q",  "FCCAnalyses::ReconstructedParticle::get_charge(RP)")

            # ---- recluster exclusive 2 jets in e+e- mode (Durham/ee-kT) ----
            .Define("fj_pseudojets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
            .Define("fj_clusters",    "JetClustering::clustering_ee_kt(2, 2, 1, 0)(fj_pseudojets)")
            .Define("jets",           "FCCAnalyses::JetClusteringUtils::get_pseudoJets(fj_clusters)")

            # per-jet angles, momenta, energies
            .Define("jet_theta", "FCCAnalyses::JetClusteringUtils::get_theta(jets)")
            .Define("jet_px",    "FCCAnalyses::JetClusteringUtils::get_px(jets)")
            .Define("jet_py",    "FCCAnalyses::JetClusteringUtils::get_py(jets)")
            .Define("jet_pz",    "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
            .Define("jet_e",     "FCCAnalyses::JetClusteringUtils::get_e(jets)")
            .Define("jet_constits", "FCCAnalyses::JetClusteringUtils::get_constituents(fj_clusters)")

            # cos(theta) and |p| for each jet
            .Define("cos0", " (jet_theta.size()>0) ? std::cos(jet_theta[0]) : -2.f")
            .Define("cos1", " (jet_theta.size()>1) ? std::cos(jet_theta[1]) : -2.f")
            .Define("p0",   " (jet_px.size()>0) ? std::sqrt(jet_px[0]*jet_px[0]+jet_py[0]*jet_py[0]+jet_pz[0]*jet_pz[0]) : 0.f")
            .Define("p1",   " (jet_px.size()>1) ? std::sqrt(jet_px[1]*jet_px[1]+jet_py[1]*jet_py[1]+jet_pz[1]*jet_pz[1]) : 0.f")

            # ---- jet charge (charged-only JC3 with p_min) ----
            .Define("jet_Q", f"bafb_jet_charge(RP_q, RP_p, jet_constits, {kappa}, {pmin})")
            .Define("Q0",    " (jet_Q.size()>0) ? jet_Q[0] : 0.f")
            .Define("Q1",    " (jet_Q.size()>1) ? jet_Q[1] : 0.f")
            .Define("DeltaQ"," Q0 - Q1")

            # charged multiplicity (can add pmin variant later if desired)
            .Define("jet_nch","bafb_jet_nch(RP_q, jet_constits)")

            # ---- fiducial cuts on BOTH jets ----
            .Filter(f"fabs(cos0) < {cth_max} && fabs(cos1) < {cth_max}", "|cos(theta_jet)| acceptance on both jets")
            .Filter(f"(p0 > {pjet_min}) && (p1 > {pjet_min})", "jet |p| > p_min on both jets")

            # ---- optional opposite-sign requirement ----
        )
        if self.opposite_sign:
            d = d.Filter("bafb_opposite_sign(jet_Q)", "Opposite jet-charge signs (purity boost)")

        # ---- primary jet = higher energy ----
        d = (
            d.Define("primary_idx", " (jet_e.size()>1) ? ((jet_e[0] >= jet_e[1]) ? 0 : 1) : 0 ")
             .Define("cos_primary"," (primary_idx==0) ? cos0 : cos1 ")
             .Define("Q_primary",  " (primary_idx==0) ? Q0   : Q1 ")
             .Define("signed_cos", " (Q_primary>0.f ? 1.f : (Q_primary<0.f ? -1.f : 0.f)) * cos_primary ")
        )

        # ---- truth helpers (closure) ----
        d = (
            d.Define("MC_pdg",   "FCCAnalyses::MCParticle::get_pdg(MC)")
             .Define("MC_theta", "FCCAnalyses::MCParticle::get_theta(MC)")
             .Define("MC_p",     "FCCAnalyses::MCParticle::get_p(MC)")
             .Define("is_b",     "bafb_is_b(MC_pdg)")
             .Define("n_b",      "bafb_sum_mask(is_b)")
             .Define("b_idx",    "bafb_argmax_masked(MC_p, is_b)")
             .Define("cth_truth","(n_b>0 && b_idx >= 0 && b_idx < (int)MC_theta.size()) ? std::cos(MC_theta[b_idx]) : -2.f")
        )

        # ---- s'/s proxy from reconstructed four-momentum ----
        d = (
            d.Define("RP_px_sum","ROOT::VecOps::Sum(RP_px)")
             .Define("RP_py_sum","ROOT::VecOps::Sum(RP_py)")
             .Define("RP_pz_sum","ROOT::VecOps::Sum(RP_pz)")
             .Define("RP_e_sum", "ROOT::VecOps::Sum(RP_e)")
             .Define("sprime_over_s",
                     " (RP_e_sum>0) ? ((RP_e_sum*RP_e_sum - (RP_px_sum*RP_px_sum+RP_py_sum*RP_py_sum+RP_pz_sum*RP_pz_sum)) / (91.1876*91.1876)) : 0.0")
        )
        if self.sprime_cut > 0.0:
            d = d.Filter(f"sprime_over_s >= {self.sprime_cut}", "s'/s pole selection")

        # ---- config echo as scalars (to snapshot) ----
        d = (
            d.Define("cfg_kappa",            f"double({self.kappa})")
             .Define("cfg_track_p_min",      f"double({self.track_p_min})")
             .Define("cfg_abs_costheta_max", f"double({self.abs_costheta_max})")
             .Define("cfg_jet_p_min",        f"double({self.jet_p_min})")
             .Define("cfg_opposite_sign",    f"int({1 if self.opposite_sign else 0})")
             .Define("cfg_sprime_cut",       f"double({self.sprime_cut})")
        )

        return d

    def output(self):
        """
        Columns to Snapshot into the stage-1 slim tree.
        The output filename is set by the CLI (--output).
        """
        return [
            # Primary signed observable + components
            "signed_cos", "cos_primary", "Q_primary", "primary_idx",
            # Per-jet monitoring
            "Q0", "Q1", "DeltaQ", "cos0", "cos1", "p0", "p1", "jet_nch",
            # Truth for closure
            "cth_truth",
            # ISR proxy
            "sprime_over_s",
            # Config echo
            "cfg_kappa", "cfg_track_p_min", "cfg_abs_costheta_max", "cfg_jet_p_min", "cfg_opposite_sign", "cfg_sprime_cut",
        ]
