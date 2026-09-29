#!/usr/bin/env python3
"""Stage 11. Rank within the small-molecule modality, place the controls, read
the null floor. The only stage that opens labels.tsv.

There is one modality here, so there is no cross-modality comparison to
refuse. The script still writes a rank percentile and a distance from the
control and never a pooled score.

Ranking pool. Env-directed candidates plus the two controls. The null molecules
are not in the pool, because a pool that mixes them in would make every
percentile depend on how many nulls were drawn. They enter separately as the
floor: (a) the AUC for separating the pool from the null on Vina score, and
(b) where each control falls inside the null score distribution.

Percentile. Rank 1 is the lowest (most negative) Vina score. Percentile is
rank / pool size, so 5 % means the top twentieth. Ties are broken by id, which
is arbitrary.

Conservation. Two rankings are written for each state and arm. 'unfiltered' is
the ranking to answer question 1 (does the pipeline find a known inhibitor at
all). 'filtered' removes escape-prone poses first, as the brief requires for
the list you would order from. Reporting both is what makes it visible when the
filter removes a control.

Glycan effect. A candidate is 'lost to glycans' if its best Vina score is worse
by more than GLYCAN_LOSS units in the glycosylated arm than in the glycan-free
arm. The 1.0 is arbitrary and chosen by eye, and the table also carries the raw
difference so a reader can pick a different line. Only glycans that the
deposition resolved are present (a handful of residues, see receptor_summary.json),
so this counts a lower bound.
"""
import argparse, itertools, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
GLYCAN_LOSS = 1.0
SEED = 20260929


def auc_lower_better(pos, neg):
    # probability a random pool member scores lower (better) than a random null molecule
    pos, neg = np.asarray(pos), np.asarray(neg)
    return float(np.mean([(p < neg).mean() + 0.5 * (p == neg).mean() for p in pos]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.parse_args()
    R = ROOT / "results"
    sc = pd.read_csv(R / "dock_scores.tsv", sep="\t")
    lab = pd.read_csv(ROOT / "data" / "candidates" / "labels.tsv", sep="\t")
    cons = pd.read_csv(R / "conservation_contacts.tsv", sep="\t")
    d = sc.merge(lab[["id", "role", "chembl", "nM"]], on="id").merge(
        cons[["id", "state", "arm", "escape_prone", "mean_entropy_bits", "n_contacts"]], on=["id", "state", "arm"], how="left")
    d["escape_prone"] = d.escape_prone.fillna(False).astype(bool)
    pool = d[d.role != "null"]
    null = d[d.role == "null"]
    rng = np.random.default_rng(SEED)

    long, ctl, funnel, floor = [], [], [], []
    for (state, arm), g in pool.groupby(["state", "arm"]):
        n = g.copy().sort_values(["score", "id"])
        n["rank_unfiltered"] = np.arange(1, len(n) + 1)
        n["pctl_unfiltered"] = 100 * n.rank_unfiltered / len(n)
        f = n[~n.escape_prone].copy()
        f["rank_filtered"] = np.arange(1, len(f) + 1)
        f["pctl_filtered"] = 100 * f.rank_filtered / len(f)
        n = n.merge(f[["id", "rank_filtered", "pctl_filtered"]], on="id", how="left")
        t = n[n.role.str.startswith("control:temsavir")]
        ctl_score = float(t.score.iloc[0]) if len(t) else np.nan
        n["delta_vs_temsavir"] = n.score - ctl_score
        long.append(n)
        gn = null[(null.state == state) & (null.arm == arm)]
        for _, r in n[n.role.str.startswith("control")].iterrows():
            inside_null = 100 * float((gn.score < r.score).mean()) if len(gn) else np.nan
            ctl.append(dict(state=state, arm=arm, control=r.role.split(":")[1], score=r.score, heavy_atoms=r.heavy_atoms,
                            pool_size=len(n), rank_unfiltered=int(r.rank_unfiltered), pctl_unfiltered=round(r.pctl_unfiltered, 1),
                            escape_prone=bool(r.escape_prone), rank_filtered=r.rank_filtered, pctl_filtered=r.pctl_filtered,
                            pctl_within_null=round(inside_null, 1), null_size=len(gn)))
        funnel.append(dict(state=state, arm=arm, docked_pool=len(n), removed_by_conservation=int(n.escape_prone.sum()),
                           remaining=len(f), null_docked=len(gn), null_removed_by_conservation=int(gn.escape_prone.sum())))
        boots = [auc_lower_better(rng.choice(n.score.values, len(n)), rng.choice(gn.score.values, len(gn))) for _ in range(300)] if len(gn) else []
        floor.append(dict(state=state, arm=arm, auc_pool_vs_null=round(auc_lower_better(n.score, gn.score), 3) if len(gn) else np.nan,
                          auc_ci95_lo=round(float(np.percentile(boots, 2.5)), 3) if boots else np.nan,
                          auc_ci95_hi=round(float(np.percentile(boots, 97.5)), 3) if boots else np.nan,
                          spearman_score_vs_heavy_atoms_pool=round(float(spearmanr(n.score, n.heavy_atoms)[0]), 3),
                          spearman_score_vs_heavy_atoms_null=round(float(spearmanr(gn.score, gn.heavy_atoms)[0]), 3) if len(gn) else np.nan,
                          median_score_pool=round(float(n.score.median()), 2), median_score_null=round(float(gn.score.median()), 2) if len(gn) else np.nan))
    L = pd.concat(long)
    L.to_csv(R / "ranks_long.tsv", sep="\t", index=False)
    pd.DataFrame(ctl).to_csv(R / "control_ranks.tsv", sep="\t", index=False)
    pd.DataFrame(funnel).to_csv(R / "funnel.tsv", sep="\t", index=False)
    pd.DataFrame(floor).to_csv(R / "null_floor.tsv", sep="\t", index=False)

    # rank correlation between states, per arm, on candidates docked in both
    rc = []
    for arm, g in pool.groupby("arm"):
        w = g.pivot(index="id", columns="state", values="score").dropna()
        for a, b in itertools.combinations(w.columns, 2):
            rho = spearmanr(w[a], w[b])[0]
            rc.append(dict(arm=arm, state_a=a, state_b=b, n=len(w), spearman=round(float(rho), 3)))
    pd.DataFrame(rc).to_csv(R / "state_rank_correlation.tsv", sep="\t", index=False)

    # divergence between states for each control: top decile in one state and bottom half in another
    pv = L.pivot_table(index=["id", "role", "arm"], columns="state", values="pctl_unfiltered").reset_index()
    sc_cols = [c for c in pv.columns if c.startswith("S")]
    pv["state_disagreement"] = pv[sc_cols].apply(lambda r: bool((r <= 10).any() and (r > 50).any()), axis=1)
    pv.to_csv(R / "state_percentiles.tsv", sep="\t", index=False)

    # glycan effect
    w = pool.pivot_table(index=["id", "state"], columns="arm", values="score").dropna().reset_index()
    w["delta_glycosylated_minus_free"] = w.glycosylated - w.glycan_free
    w["lost_to_glycans"] = w.delta_glycosylated_minus_free > GLYCAN_LOSS
    w = w.merge(lab[["id", "role"]], on="id")
    w.to_csv(R / "glycan_effect.tsv", sep="\t", index=False)
    print(pd.DataFrame(ctl).to_string(), flush=True)


if __name__ == "__main__":
    main()
