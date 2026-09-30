#!/usr/bin/env python3
"""Stage 12. What can be tested against data the docking never saw.

The brief's held-out test is CATNAP neutralisation data against the
protein-binder arm. That arm was not run, and CATNAP (hiv.lanl.gov/catnap) is
served through a web form that was not downloaded here, so the CATNAP
retrospective and the twelve-virus panel coverage were not run. This script
writes that statement to results/holdout_not_run.md and does two smaller tests
that the small-molecule arm can support.

1. Potency association. Among the candidates that carry a ChEMBL potency,
   Spearman correlation between Vina score and -log10(potency in nM). The
   potencies come from many assays (different viruses, cells, readouts), and
   ChEMBL was searched with the same papers that defined the candidate set, so
   this is not a blind test in the strict sense: the potency values were never
   used by the docking, and that is all it is. A residual version removes the
   linear trend of score on heavy-atom count, because Vina score grows with
   size and so does the potency of medicinal-chemistry series that were
   optimised by adding groups.

2. Redocking of the temsavir crystal pose. Docking temsavir into the 5U7O
   receptor it came out of is the easiest case there is. If the top pose does
   not overlap the crystal ligand, the box or the preparation is wrong. Reported
   as centroid distance and as the fraction of crystal heavy atoms within 2 A of
   a pose heavy atom (atom-matched RMSD needs a graph mapping that the
   coordinate file does not carry).
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def pose_xyz(f):
    xyz = []
    for l in Path(f).read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")) and not l[77:79].strip().startswith("H"):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if l.startswith("ENDMDL"):
            break
    return np.array(xyz)


def main():
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args()
    R = ROOT / "results"
    sc = pd.read_csv(R / "dock_scores.tsv", sep="\t")
    lab = pd.read_csv(DATA / "candidates" / "labels.tsv", sep="\t")
    d = sc.merge(lab, on="id")
    cand = d[(d.role != "random_null") & d.nM.notna() & (d.nM > 0)].copy()
    cand["pot"] = -np.log10(cand.nM)
    rows = []
    for (s, a), g in cand.groupby(["state", "arm"]):
        b = np.polyfit(g.heavy_atoms, g.score, 1)
        res = g.score - np.polyval(b, g.heavy_atoms)
        rows.append(dict(state=s, arm=a, n=len(g), spearman_score_vs_pot=round(float(spearmanr(g.score, g.pot)[0]), 3),
                         spearman_size_adjusted=round(float(spearmanr(res, g.pot)[0]), 3),
                         spearman_heavy_vs_pot=round(float(spearmanr(g.heavy_atoms, g.pot)[0]), 3)))
    pd.DataFrame(rows).to_csv(R / "potency_association.tsv", sep="\t", index=False)

    ref = np.array(json.loads((DATA / "panel" / "boxes.json").read_text())["reference_ligand_xyz"])
    t = lab[lab.role == "control:temsavir_83J"].id.iloc[0]
    out = []
    for arm in ("glycan_free", "glycosylated"):
        f = DATA / "poses" / f"{t}__S2_closed_holo__{arm}.pdbqt"
        if not f.exists():
            continue
        p = pose_xyz(f)
        cover = float(np.mean(cKDTree(p).query(ref)[0] < 2.0))
        out.append(dict(arm=arm, centroid_distance_A=round(float(np.linalg.norm(p.mean(0) - ref.mean(0))), 2),
                        crystal_atoms_within_2A_of_pose=round(cover, 2), crystal_heavy_atoms=len(ref), pose_heavy_atoms=len(p)))
    pd.DataFrame(out).to_csv(R / "redock_temsavir_S2.tsv", sep="\t", index=False)

    (R / "holdout_not_run.md").write_text(
        "# Held-out tests not run\n\n"
        "The CATNAP retrospective test and the twelve-virus panel coverage were not run. Both belong to the "
        "protein-binder arm, which has no candidates (stage 5 was never run), and the CATNAP export is a web-form "
        "download that was not obtained. No CATNAP record count or download date exists for this run.\n\n"
        "What was run instead is in potency_association.tsv (docking score against ChEMBL potency) and "
        "redock_temsavir_S2.tsv (temsavir back into its own crystal receptor).\n")


if __name__ == "__main__":
    main()
