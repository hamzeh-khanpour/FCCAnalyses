# analysis_stage2.py
# Stage-2: build histograms and A_FB from the Stage-1 slim tree.
# - Filters out sentinel cosTheta values
# - Computes N_F / N_B, A_FB (raw) and (optionally) corrected using purity P and charge confusion ω
# - Writes ROOT histograms and a JSON summary (with provenance: threads, note)

import os, json, math, ROOT

output_dir = "outputs/bafb/stage2"

class Analysis:
    def __init__(self, args):
        os.makedirs(output_dir, exist_ok=True)
        # Let FCCAnalyses manage threads; we just record what was requested.
        self.n_threads = int(getattr(args, "n_threads", 4))
        self.input_path = getattr(args, "input", None)

        # record basic config for provenance
        self.cfg = {
            "n_threads": self.n_threads,
            "input": self.input_path,
            "note": "Purity (P) and omega (ω) are placeholders; update them from MC studies for corrected A_FB."
        }
        print(f"[stage2] (framework-managed threads), input={self.input_path}")

    def analyzers(self, df):
        d = (
            df
            .Filter("cosTheta>-1.5f")          # drop sentinel
            .Define("isF","cosTheta>0.f")
            .Define("isB","cosTheta<0.f")
            .Define("abscth","std::abs(cosTheta)")
        )

        # Histos & counters
        self.h_cos = d.Histo1D(("h_cos","signed cos#theta;cos#theta;Events",20,-1,1),"cosTheta")
        self.h_abs = d.Histo1D(("h_abs","|cos#theta|;|cos#theta|;Events",20,0,1),"abscth")

        self.nF = d.Sum("isF")
        self.nB = d.Sum("isB")

        # Optional truth histogram (present only if stage-1 stored it)
        self.h_cth_truth = None
        try:
            self.h_cth_truth = d.Filter("cth_truth>-1.5f").Histo1D(("h_cth_truth","truth cos#theta_{b};cos#theta_{b};Events",20,-1,1),"cth_truth")
        except Exception:
            pass

        return d

    def output(self):
        # Pull numbers
        nf, nb = int(self.nF.GetValue()), int(self.nB.GetValue())
        n = nf + nb
        if n <= 0:
            a_raw, da_raw = 0.0, 0.0
        else:
            a_raw = (nf - nb) / n
            da_raw = math.sqrt(max(0.0, (1.0 - a_raw * a_raw) / n))

        # Placeholder corrections (update with MC-derived values)
        purity = 1.0
        omega  = 0.0
        scale  = purity * (1.0 - 2.0 * omega)
        if scale != 0.0:
            a_corr = a_raw / scale
            da_corr = abs(da_raw / scale)
        else:
            a_corr, da_corr = 0.0, 0.0

        # Save ROOT outputs
        fout = ROOT.TFile(os.path.join(output_dir,"bafb_outputs.root"),"RECREATE")
        self.h_cos.Write(); self.h_abs.Write()
        if self.h_cth_truth: self.h_cth_truth.Write()
        fout.Close()

        # Write JSON summary (with provenance)
        meta = {
            "N_F": nf, "N_B": nb, "N_tot": n,
            "A_FB_raw": a_raw, "A_FB_raw_err": da_raw,
            "purity_assumed": purity, "omega_assumed": omega,
            "A_FB_corr": a_corr, "A_FB_corr_err": da_corr
        }
        meta.update(self.cfg)
        with open(os.path.join(output_dir,"summary.json"),"w") as f:
            json.dump(meta, f, indent=2)

        # Minimal list to satisfy managed mode snapshot API (not used further)
        return ["cosTheta","isF","isB","abscth"]
