## AFBb @ Z Pole — Exclusive `Z → b b̄` (No Backgrounds)

**Goal:** measure the forward–backward asymmetry of b quarks, (A_{FB}^b), at the Z pole using an **exclusive** (Z\to b\bar b) sample with **no backgrounds**.
**Flow:** Stage-1 (feature building) → Stage-2 (binned Poisson fit, plots, JSON/TXT).

---

## Quick Start

Running analysis scripts is done with the `fccanalysis` command shipped in the Key4hep stack.

```sh
# 1) Environment
source /cvmfs/sw.hsf.org/key4hep/setup.sh

# 2) End-to-end run (Stage-1 → Stage-2)
chmod +x analysis_3/bafb/python/run.sh
./analysis_3/bafb/python/run.sh --input-file-list my_bbbar_files.txt --n-threads 100
```

Expected outputs:

```
outputs/bafb/
├─ stage1/
│  └─ bafb_stage1.root
└─ stage2/
   ├─ afbb_stage2.root            # contains TDirectory "AFBb" with histograms
   ├─ plots/
   │  ├─ AFBb_signed_fit.png
   │  └─ AFBb_signed_fit.pdf
   └─ results/
      ├─ AFBb_fit.json
      └─ AFBb_fit.txt
```

---

## Pre-generated Samples

### Access

To have read access to the FCC pre-generated samples, one needs to be subscribed to
the following e-group (with owner approval): `fcc-eos-access`.

### Example file list

Point the runner to a text file with ROOT paths (XRootD or local):

```
root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91/events_000083138.root
root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91/events_000101935.root
...
```

Use it via:

```sh
./analysis_3/bafb/python/run.sh --input-file-list my_bbbar_files.txt --n-threads 100
```

---

## Physics Overview (short)

At the Z pole (unpolarized beams), the b-quark polar angle w.r.t. the **electron-beam** direction ((+\hat z)) follows

```math
\frac{dN}{d\cos\theta} \propto (1+\cos^2\theta) + \frac{8}{3}\,A_{FB}^b\,\cos\theta.
```

We:

1. Cluster jets with **Durham, N=2**.
2. Compute **jet charge** (Q_J(\kappa)) (JC3, charged tracks).
3. Choose a **primary** jet (higher energy) and build the **signed** observable
   (x=\mathrm{sign}(Q_{\rm primary})\cdot\cos\theta_{\rm primary}).
4. Fit the signed distribution to extract (A_{FB}^b), correcting for jet-charge dilution (\omega):
   (A_{FB}^{\rm meas}=(1-2\omega),A_{FB}^{\rm true}).

**Convention:** (+\hat z) = electron beam; lab (\approx) CM at (\sqrt{s}\simeq m_Z).
**Assumption:** exclusive (Z\to b\bar b), **no backgrounds**.

---

## Repository Layout

```
analysis_3/bafb/python/
├─ analysis_stage1.py      # Stage-1: clustering, jet charge, signed cosθ, truth hooks, s'/s
├─ analysis_stage2.py      # Stage-2: binned fit, acceptance folding (optional), plots, JSON/TXT
├─ analysis_plots.py       # Optional monitors / overlays
├─ utils_bafb.py           # C++ helpers (jet charge, cosθ) + Python utilities
├─ run.sh                  # End-to-end runner
└─ my_bbbar_files.txt      # Example input file list
```

---

## Run Stage-1 only

```sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh
fccanalysis run analysis_3/bafb/python/analysis_stage1.py \
  --input-file-list my_bbbar_files.txt \
  --n-threads 100 \
  --output outputs/bafb/stage1/bafb_stage1.root
```

**Key Stage-1 options (may vary by version):**

* `--kappa 0.3` : JC3 exponent (default `0.3`)
* `--track-p-min 0.5` : track (p_{\min}) [GeV] (default `0.5`)
* `--opposite-sign 0|1` : require opposite jet-charge signs (default `0`)
* `--sprime-cut <x>` : minimal (s'/s) (default `0.0`)
* `--sqrts 91.1876` : (\sqrt{s}) for (s'/s) denominator (default (m_Z))

---

## Run Stage-2 only

```sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh
fccanalysis run analysis_3/bafb/python/analysis_stage2.py \
  --input outputs/bafb/stage1/bafb_stage1.root \
  --n-threads 1 \
  --nbins 20 \
  --omega 0.20 --omega_sigma 0.00 \
  --sprime_min -1.0 \
  --fold_acceptance ""
```

**Stage-2 options:**

* `--nbins <int>` : signed (\cos\theta) bins (default `20`)
* `--omega <float>` : wrong-sign prior mean (if truth absent) (default `0.20`)
* `--omega_sigma <float>` : Gaussian prior width; if `>0`, profile (\omega)
* `--sprime_min <float>` : minimal (s'/s) (default `-1.0` → off)
* `--fold_acceptance <json>` : path to acceptance file with

  ```json
  { "edges":[...], "eff":[...] }  // len(edges)=nbins+1, len(eff)=nbins
  ```

---

## Outputs (Stage-2)

* **ROOT:** `outputs/bafb/stage2/afbb_stage2.root`
  TDirectory **`AFBb/`** contains:

  * `h_signed_cosTheta`
  * `h_fit_overlay`
  * (monitors) e.g. `h_Q_primary`, `resp_truth_vs_signed`, etc.
* **Plots:** `outputs/bafb/stage2/plots/AFBb_signed_fit.{png,pdf}`
* **Results:** `outputs/bafb/stage2/results/AFBb_fit.{json,txt}` (includes provenance)

---

## Troubleshooting

* **Null hist / canvas draw errors:**
  We avoid slashes in histogram **names** and write under a real `AFBb/` directory.
  Histograms used after file close are detached with:

  ```cpp
  h->SetDirectory(0);
  ```
* **Overlay text shows `Â±`:**
  The canvas uses Unicode ± (`\u00B1`); TXT uses `+/-`.

---

## Validation (quick checklist)

* Truth vs reco closure on (A_{FB}) (if truth branches present)
* (\omega(|x|)) monitoring
* (\kappa) scan (0.2–0.5)
* ISR window scan via `--sprime_min`
* (Optional) acceptance folding JSON consistency (`edges/eff` lengths)

---

## License / Citation

Please cite FCCAnalyses, the Key4hep stack, and FastJet where appropriate.
If you use this analysis or parts of it, reference this repository and branch.
