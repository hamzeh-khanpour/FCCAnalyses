# AFBb @ Z Pole — Exclusive `Z → b b̄` (No Backgrounds)

> Forward–backward asymmetry of b quarks, (A_{FB}^b), measured at the Z pole using an **exclusive** (Z\to b\bar b) sample with **no backgrounds**.
> Pipeline: **Stage-1** (feature building) → **Stage-2** (fit & reporting).

---

## Table of Contents

* [Physics Overview](#physics-overview)
* [Repository Layout](#repository-layout)
* [Environment](#environment)
* [Quick Start](#quick-start)
* [Commands](#commands)

  * [Run everything](#run-everything)
  * [Run Stage-1 only](#run-stage1-only)
  * [Run Stage-2 only](#run-stage2-only)
* [Configuration](#configuration)

  * [Stage-1 knobs](#stage1-knobs)
  * [Stage-2 knobs](#stage2-knobs)
* [Outputs](#outputs)
* [Method Details](#method-details)

  * [Signed angle and dilution](#signed-angle-and-dilution)
  * [Fit model](#fit-model)
  * [Acceptance folding (optional)](#acceptance-folding-optional)
  * [Expected precision](#expected-precision)
* [Validation Checklist](#validation-checklist)
* [Troubleshooting](#troubleshooting)
* [Reproducibility](#reproducibility)
* [Citations & Acknowledgements](#citations--acknowledgements)

---

## Physics Overview

At the Z pole (unpolarized beams), the polar angle of the **b quark** relative to the **electron-beam** direction ((+\hat z)) follows
[
\frac{dN}{d\cos\theta}\ \propto\ (1+\cos^2\theta);+;\frac{8}{3},A_{FB}^b,\cos\theta .
]

We:

1. Reconstruct two jets (Durham, **N=2**).
2. Compute **jet charge** (Q_J(\kappa)) (JC3, charged tracks only).
3. Choose a **primary** jet, compute (\cos\theta_{\rm primary}) w.r.t. (+\hat z) (electron beam).
4. Build the **signed** observable (x=\operatorname{sign}(Q_{\rm primary})\cdot \cos\theta_{\rm primary}).
5. Fit the signed distribution to extract (A_{FB}^b).

**Sign convention:** (+\hat z) along the **electron** beam; lab (\approx) CM at (\sqrt{s}\simeq m_Z).
**Assumption:** **exclusive** (Z\to b\bar b), **no backgrounds**.

---

## Repository Layout

```
analysis_3/bafb/python/
├─ analysis_stage1.py      # Stage-1: clustering, jet charge, signed cosθ, truth hooks, s'/s
├─ analysis_stage2.py      # Stage-2: binned fit, acceptance folding, plots, JSON/TXT
├─ analysis_plots.py       # Optional monitors / overlays
├─ utils_bafb.py           # C++ helpers (jet charge, cosθ) + Python utilities
├─ run.sh                  # End-to-end runner (Stage-1 → Stage-2)
└─ my_bbbar_files.txt      # Example input file list (xrootd/paths)
```

**Outputs**

```
outputs/bafb/
├─ stage1/
│  └─ bafb_stage1.root
└─ stage2/
   ├─ afbb_stage2.root          # contains TDirectory "AFBb"
   ├─ plots/
   │  ├─ AFBb_signed_fit.png
   │  └─ AFBb_signed_fit.pdf
   └─ results/
      ├─ AFBb_fit.json
      └─ AFBb_fit.txt
```

---

## Environment

This analysis is designed for the **Key4hep** stack with **FCCAnalyses**.

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
# (Optionally choose the release matching your campaign.)
```

---

## Quick Start

```bash
chmod +x analysis_3/bafb/python/run.sh

# Use the provided list or create your own.
./analysis_3/bafb/python/run.sh --input-file-list my_bbbar_files.txt --n-threads 100
```

* **Stage-1** produces `outputs/bafb/stage1/bafb_stage1.root`
* **Stage-2** produces fit results, plots, and `outputs/bafb/stage2/afbb_stage2.root`

---

## Commands

### Run everything

```bash
./analysis_3/bafb/python/run.sh --input-file-list my_bbbar_files.txt --n-threads 100
```

### Run Stage-1 only

```bash
fccanalysis run analysis_3/bafb/python/analysis_stage1.py \
  --input-file-list my_bbbar_files.txt \
  --n-threads 100 \
  --output outputs/bafb/stage1/bafb_stage1.root
```

### Run Stage-2 only

```bash
fccanalysis run analysis_3/bafb/python/analysis_stage2.py \
  --input outputs/bafb/stage1/bafb_stage1.root \
  --n-threads 1 \
  --nbins 20 \
  --omega 0.20 --omega_sigma 0.00 \
  --sprime_min -1.0 \
  --fold_acceptance ""
```

> Use `--n-threads 1` for deterministic logs; multithreading isn’t required in Stage-2.

---

## Configuration

### Stage-1 knobs

* `--kappa 0.3` : JC3 exponent (default `0.3`)
* `--track-p-min 0.5` : track momentum threshold in GeV (default `0.5`)
* `--opposite-sign 0|1` : require opposite jet-charge signs (default `0`)
* `--sprime-cut <x>` : minimal (s'/s) (default `0.0`)
* `--sqrts 91.1876` : (\sqrt{s}) used in (s'/s) denominator (default `m_Z`)

### Stage-2 knobs

* `--nbins <int>` : signed (\cos\theta) bins (default `20`)
* `--omega <float>` : wrong-sign prior mean if truth absent (default `0.20`)
* `--omega_sigma <float>` : Gaussian prior width; if `>0`, profile (\omega)
* `--sprime_min <float>` : minimal (s'/s) at Stage-2 selection (default `-1.0` → off)
* `--fold_acceptance <path.json>` : fold (\varepsilon(x)) into the template (optional)

---

## Outputs

* **ROOT**: `outputs/bafb/stage2/afbb_stage2.root` with a real directory `AFBb/`:

  * `AFBb/h_signed_cosTheta`
  * `AFBb/h_fit_overlay`
  * monitors (no slashes in object names)
* **Plots**: `outputs/bafb/stage2/plots/AFBb_signed_fit.{png,pdf}`
* **Results**: `outputs/bafb/stage2/results/AFBb_fit.{json,txt}`

---

## Method Details

### Signed angle and dilution

* **Jet charge (JC3, charged only)**
  [
  Q_J(\kappa)=\frac{\sum_{i\in J} q_i,|\vec p_i|^\kappa}{\sum_{i\in J} |\vec p_i|^\kappa},\quad
  \kappa=0.3\ \text{(default)}.
  ]
* **Signed observable**
  (x=\operatorname{sign}(Q_{\rm primary})\cdot \cos\theta_{\rm primary}), where the **primary** jet is the higher-energy jet.
* **Dilution**
  [
  A_{FB}^{\rm meas}=(1-2\omega),A_{FB}^{\rm true},\qquad
  A_{FB}^{\rm true}=\frac{A_{FB}^{\rm meas}}{1-2\omega}.
  ]
  If truth tags exist, (\omega) is measured; otherwise use `--omega` prior (optionally profiled with `--omega_sigma`).

### Fit model

We fit the signed distribution with
[
f(x\mid A_{FB}^{\rm true},\omega)\ \propto\ (1+x^2)+\frac{8}{3}(1-2\omega),A_{FB}^{\rm true},x,
\quad x\in[-1,1],
]
using a **binned Poisson** likelihood (default 20 bins).

### Acceptance folding (optional)

Supply a JSON:

```json
{
  "edges": [-1.0, -0.9, ..., 1.0],
  "eff":   [0.98, 0.99, ..., 0.97]
}
```

with `len(edges)=nbins+1`, `len(eff)=nbins`. The template is multiplied by `eff` per bin and renormalized.

### Expected precision

Statistical (no dilution):
[
\sigma(A_{FB}) \simeq \frac{\sqrt{1-A_{FB}^2}}{\sqrt{N}}.
]
With dilution (\omega), inflate by (1/(1-2\omega)).
Example: (N=10^6), (A\approx0.1), (\omega=0.20\Rightarrow(1-2\omega)=0.6) → (\sigma(A)\approx 0.0017).

---

## Validation Checklist

* **Truth vs reco closure:** compare (A_{FB}^{\rm true}) from truth-tag vs reconstructed signed histogram.
* **(\omega(|x|)) shape:** monitor wrong-sign fraction vs (|x|).
* **(\kappa) scan:** (0.2)–(0.5); choose default on precision × robustness.
* **Axis robustness:** (if enabled) thrust vs Durham axis signing.
* **ISR window:** vary (s'/s) cuts and check stability vs efficiency.
* **Toys / bootstrap:** pull width (\sim1); coverage at the quoted error.

---

## Troubleshooting

* **`CPyCppyy_NoneType has no attribute 'Draw'`**
  ROOT histogram **names must not contain “/”**. We create a real `AFBb/` directory in the ROOT file and use **slash-free object names**.
  Also, histograms drawn **after closing** the file are detached via:

  ```python
  h.SetDirectory(0)
  ```

* **Overlay text shows `Â±`**
  Canvas text uses Unicode `±` (`\u00B1`). TXT files use ASCII `+/-`.

* **Empty histograms**
  Relax `--sprime_min` or fiducials; verify Stage-1 outputs exist and contain `signed_cos`.

* **Acceptance JSON rejected**
  Ensure `edges`/`eff` lengths match the configured `nbins`.

---

## Reproducibility

* Deterministic RDataFrame flow (no RNG).
* Stage-2 writes a **provenance** block in JSON (includes Stage-1 `cfg_*` echoes when present).

---

## Citations & Acknowledgements

* **FCCAnalyses**, **Key4hep** software stack.
* **FastJet** for jet clustering.
* FCC-ee community inputs on (A_{FB}^b) methodology.

> If you use this analysis or parts of it, please cite FCCAnalyses and relevant FCC-ee documentation, and reference this repository.
