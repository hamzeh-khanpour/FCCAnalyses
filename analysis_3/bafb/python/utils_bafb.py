# utils_bafb.py
# Helpers for jet charge, AFB, and simple corrections used by the b-AFB pipeline.
# NOTE: JC3 here uses ONLY charged constituents; neutrals are ignored in both
# numerator and denominator. Default kappa is passed by Stage-1; recommended kappa=0.3.
# Momentum units are assumed to be GeV.
#
# Robust C++ helpers are injected for RDataFrame Defines to avoid JIT type erasure.

import math
import ROOT

# --- C++ helpers injected for RDataFrame Defines -----------------------------------------------
ROOT.gInterpreter.Declare(r"""
#include <vector>
#include <cmath>
#include <cstddef>
#include <ROOT/RVec.hxx>

using Vfloat     = ROOT::VecOps::RVec<float>;
using Vint       = ROOT::VecOps::RVec<int>;
using VVintRVec  = ROOT::VecOps::RVec< ROOT::VecOps::RVec<int> >;

// ============================================================================
// Core implementations (charged-only JC3)
// Q_J = sum_i (q_i * |p_i|^kappa) / sum_i (|p_i|^kappa), with selection:
//   - track is CHARGED: rp_q[idx] != 0
//   - |p| >= pmin (pmin can be 0)
// ============================================================================

template <class Constituents>
Vfloat bafb_jet_charge_impl_charged(const Vfloat& rp_q,
                                    const Vfloat& rp_p,
                                    const Constituents& constituents,
                                    const double kappa,
                                    const double pmin)
{
  Vfloat out; out.reserve(constituents.size());
  for (const auto& idxs : constituents){
    double num=0., den=0.;
    for (const auto& idx : idxs){
      if ((std::size_t)idx >= rp_q.size() || (std::size_t)idx >= rp_p.size()) continue; // bounds
      const float q = rp_q[idx];
      if (q == 0.f) continue;                       // skip neutrals (charged-only JC3)
      const double p = std::fabs((double)rp_p[idx]);
      if (p < pmin) continue;                       // p threshold
      const double w = std::pow(p, kappa);
      num += (double)q * w;
      den += w;
    }
    out.emplace_back( den>0.0 ? static_cast<float>(num/den) : 0.f );
  }
  return out;
}

// -------------------- Public overloads (keep 4-arg name/signature) -----------------------------

// 4-arg version (backward compatible): behaves as pmin=0.0
Vfloat bafb_jet_charge(const Vfloat& rp_q,
                       const Vfloat& rp_p,
                       const VVintRVec& constituents,
                       const double kappa)
{
  return bafb_jet_charge_impl_charged(rp_q, rp_p, constituents, kappa, 0.0);
}

Vfloat bafb_jet_charge(const Vfloat& rp_q,
                       const Vfloat& rp_p,
                       const std::vector<std::vector<int>>& constituents,
                       const double kappa)
{
  return bafb_jet_charge_impl_charged(rp_q, rp_p, constituents, kappa, 0.0);
}

// 5-arg versions with pmin threshold
Vfloat bafb_jet_charge(const Vfloat& rp_q,
                       const Vfloat& rp_p,
                       const VVintRVec& constituents,
                       const double kappa,
                       const double pmin)
{
  return bafb_jet_charge_impl_charged(rp_q, rp_p, constituents, kappa, pmin);
}

Vfloat bafb_jet_charge(const Vfloat& rp_q,
                       const Vfloat& rp_p,
                       const std::vector<std::vector<int>>& constituents,
                       const double kappa,
                       const double pmin)
{
  return bafb_jet_charge_impl_charged(rp_q, rp_p, constituents, kappa, pmin);
}

// -------------------- Charged multiplicity helpers --------------------------------------------

// Counts charged constituents (rp_q != 0) per jet
ROOT::VecOps::RVec<int> bafb_jet_nch(const Vfloat& rp_q,
                                     const std::vector<std::vector<int>>& constituents)
{
  ROOT::VecOps::RVec<int> out; out.reserve(constituents.size());
  for (const auto& idxs : constituents){
    int n=0;
    for (const auto& idx : idxs){
      if ((std::size_t)idx < rp_q.size() && rp_q[idx]!=0.f) ++n;
    }
    out.emplace_back(n);
  }
  return out;
}

// Same, but require |p| >= pmin (and charged)
ROOT::VecOps::RVec<int> bafb_jet_nch(const Vfloat& rp_q,
                                     const VVintRVec& constituents,
                                     const Vfloat& rp_p,
                                     const double pmin)
{
  ROOT::VecOps::RVec<int> out; out.reserve(constituents.size());
  for (const auto& idxs : constituents){
    int n=0;
    for (const auto& idx : idxs){
      if ((std::size_t)idx >= rp_q.size() || (std::size_t)idx >= rp_p.size()) continue;
      if (rp_q[idx]==0.f) continue;
      if (std::fabs((double)rp_p[idx]) < pmin) continue;
      ++n;
    }
    out.emplace_back(n);
  }
  return out;
}

ROOT::VecOps::RVec<int> bafb_jet_nch(const Vfloat& rp_q,
                                     const std::vector<std::vector<int>>& constituents,
                                     const Vfloat& rp_p,
                                     const double pmin)
{
  ROOT::VecOps::RVec<int> out; out.reserve(constituents.size());
  for (const auto& idxs : constituents){
    int n=0;
    for (const auto& idx : idxs){
      if ((std::size_t)idx >= rp_q.size() || (std::size_t)idx >= rp_p.size()) continue;
      if (rp_q[idx]==0.f) continue;
      if (std::fabs((double)rp_p[idx]) < pmin) continue;
      ++n;
    }
    out.emplace_back(n);
  }
  return out;
}

// -------------------- Small utilities ----------------------------------------------------------

// Return 1 if first two jets have opposite jet-charge signs, else 0.
// (If <2 jets, returns 0.)
int bafb_opposite_sign(const Vfloat& qjets) {
  if (qjets.size() < 2) return 0;
  return (qjets[0]*qjets[1] < 0.f) ? 1 : 0;
}

// Pick the b-jet index using the conventional rule: more negative jet charge => b.
// If equal, choose the one with larger |Q| (and if still equal, index 0).
int bafb_choose_bjet(const Vfloat& qjets){
  if (qjets.empty()) return -1;
  if (qjets.size()==1) return 0;
  if (qjets[0]==qjets[1]){
    // tie-break: larger |Q| wins; if still equal, index 0
    return (std::abs(qjets[0]) >= std::abs(qjets[1])) ? 0 : 1;
  }
  return (qjets[0] < qjets[1]) ? 0 : 1;
}

// --- Truth helpers (robust replacements for inline-JIT lambdas) -------------------------------

// Make a 0/1 mask for b-quarks from PDG IDs (abs(pdg)==5).
Vint bafb_is_b(const Vint& pdg){
  Vint out; out.resize(pdg.size());
  for (std::size_t i=0;i<pdg.size();++i) out[i] = (std::abs(pdg[i])==5) ? 1 : 0;
  return out;
}

// Sum of an integer mask (avoids JIT confusion with Sum on unknown types).
int bafb_sum_mask(const Vint& m){
  int s=0; for (auto v : m) s+=v; return s;
}

// Argmax over "vals" but only where mask!=0. Returns -1 if no valid elements.
int bafb_argmax_masked(const Vfloat& vals, const Vint& mask){
  if (vals.size()==0 || vals.size()!=mask.size()) return -1;
  int idx=-1; float best=-1e30f;
  for (std::size_t i=0;i<vals.size();++i){
    if (!mask[i]) continue;
    if (idx<0 || vals[i]>best){ best=vals[i]; idx=(int)i; }
  }
  return idx;
}
""")
# ------------------------------------------------------------------------------------------------


def afb_from_counts(nf: int, nb: int):
    """
    Compute A_FB and its (binomial) statistical uncertainty from integer counts.

    A = (N_F - N_B) / (N_F + N_B)
    Var(A) ~ (1 - A^2) / (N_F + N_B)

    Returns:
        (A_FB, sigma_A)
    """
    nf = float(nf)
    nb = float(nb)
    n = nf + nb
    if n <= 0:
        return 0.0, 0.0
    a = (nf - nb) / n
    err = math.sqrt(max(0.0, (1.0 - a * a) / n))
    return a, err


def apply_purity_dilution(a_raw: float,
                          purity: float,
                          omega: float,
                          err_raw: float | None = None,
                          err_purity: float | None = None,
                          err_omega: float | None = None):
    """
    Correct the observed asymmetry A_FB for sample purity (P) and charge-confusion
    dilution (1 - 2ω):

        A_true = A_raw / [ P * (1 - 2ω) ]

    If uncertainties are provided, propagate them approximately in quadrature.
    Args:
        a_raw:      observed (raw) asymmetry
        purity:     P in [0,1] (set to 1.0 if not using a b-tag yet)
        omega:      wrong-sign probability ω in [0,1]
        err_raw:    (optional) statistical error on a_raw
        err_purity: (optional) uncertainty on P
        err_omega:  (optional) uncertainty on ω
    Returns:
        (A_true, sigma_A_true)  (sigma may be None if err_raw is None)
    """
    # scale may be small; keep numerically safe
    scale = (1.0 - 2.0 * omega) * max(0.0, purity)
    if scale == 0.0:
        return 0.0, 0.0

    a_corr = a_raw / scale

    if err_raw is None:
        return a_corr, None

    # dA/dAraw = 1/scale
    da = abs(err_raw / scale)

    # dA/dP = -A_raw / scale^2  => fractional piece ~ |A_raw| / (scale * P)
    if err_purity and purity > 0:
        da = math.hypot(da, abs(a_raw) * err_purity / abs(scale * purity))

    # dA/dω = +2*A_raw / scale^2 => fractional piece ~ |2*A_raw| / scale
    if err_omega:
        da = math.hypot(da, abs(2.0 * a_raw) * err_omega / abs(scale))

    return a_corr, da


def afb_corrected(a_raw: float, sig_raw: float, P: float, omega: float,
                  sig_P: float | None = None, sig_omega: float | None = None):
    """
    Convenience wrapper used by stage-2 to emit A_FB^corr with propagated error.

        A_FB^corr = a_raw / [ P * (1 - 2ω) ]

    Parameters are identical to apply_purity_dilution; see that docstring.
    """
    return apply_purity_dilution(a_raw, P, omega, sig_raw, sig_P, sig_omega)
