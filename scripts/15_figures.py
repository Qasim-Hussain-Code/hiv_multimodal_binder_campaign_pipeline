#!/usr/bin/env python3
"""Stage 15. Figures for the README, written to figures/ from results/ tables.

Colours are an Okabe-Ito subset so the two glycan arms and the three roles stay
distinguishable without colour vision. Every figure reads only files in
results/, so it can be regenerated without re-docking.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
R, F = ROOT / "results", ROOT / "figures"
ARM_C = {"glycan_free": "#0072B2", "glycosylated": "#D55E00"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def main():
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args()
    F.mkdir(exist_ok=True)
    ctl = pd.read_csv(R / "control_ranks.tsv", sep="\t")
    fun = pd.read_csv(R / "funnel.tsv", sep="\t")
    L = pd.read_csv(R / "ranks_long.tsv", sep="\t")
    nl = pd.read_csv(R / "null_floor.tsv", sep="\t")
    sc = pd.read_csv(R / "dock_scores.tsv", sep="\t")
    lab = pd.read_csv(ROOT / "data" / "candidates" / "labels.tsv", sep="\t")
    rc = pd.read_csv(R / "state_rank_correlation.tsv", sep="\t")
    gl = pd.read_csv(R / "glycan_effect.tsv", sep="\t")
    states = sorted(ctl.state.unique())

    # 1 controls
    fig, axs = plt.subplots(1, 2, figsize=(7.5, 3), sharey=True)
    for ax, name in zip(axs, ["temsavir_83J", "BMS378806_83G"]):
        for j, arm in enumerate(["glycan_free", "glycosylated"]):
            d = ctl[(ctl.control == name) & (ctl.arm == arm)].set_index("state").reindex(states)
            x = np.arange(len(states)) + (j - 0.5) * 0.25
            ax.scatter(x, d.pctl_unfiltered, color=ARM_C[arm], label=arm, s=28, zorder=3)
            fl = d[d.escape_prone == True]
            ax.scatter(np.arange(len(states))[[states.index(s) for s in fl.index]] + (j - 0.5) * 0.25, fl.pctl_unfiltered,
                       facecolors="none", edgecolors="k", s=70, zorder=2)
        ax.axhline(50, color="0.7", lw=0.8, ls="--")
        ax.set_xticks(range(len(states))); ax.set_xticklabels([s.split("_")[0] for s in states])
        ax.set_title(name.replace("_", " "), fontsize=9); ax.set_ylim(0, 105); ax.invert_yaxis()
    axs[0].set_ylabel("rank percentile in pool (top = 0)")
    axs[1].legend(frameon=False, fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 0.75))
    fig.text(0.5, -0.02, "Dashed line: median of the pool. Lower on the axis is a worse rank. No pose was flagged escape-prone.", ha="center", fontsize=7)
    fig.tight_layout(); fig.savefig(F / "fig1_positive_controls.png", bbox_inches="tight"); plt.close(fig)

    # 2 funnel
    fig, ax = plt.subplots(figsize=(7.5, 3))
    lbl = [f"{s.split('_')[0]}\n{a[:5]}" for s, a in zip(fun.state, fun.arm)]
    ax.bar(lbl, fun.docked_pool, color="0.8", label="docked pool")
    ax.bar(lbl, fun.remaining, color="#0072B2", label="after conservation filter")
    for i, r in fun.iterrows():
        ax.text(i, r.docked_pool + 1, f"-{r.removed_by_conservation}", ha="center", fontsize=7)
    ax.set_ylabel("compounds"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(F / "fig2_funnel.png"); plt.close(fig)

    # 3 state correlation
    fig, axs = plt.subplots(1, 2, figsize=(7, 3))
    for ax, arm in zip(axs, ["glycan_free", "glycosylated"]):
        m = pd.DataFrame(np.eye(len(states)), index=states, columns=states)
        for r in rc[rc.arm == arm].itertuples():
            m.loc[r.state_a, r.state_b] = m.loc[r.state_b, r.state_a] = r.spearman
        im = ax.imshow(m.values, vmin=-0.2, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(states))); ax.set_yticks(range(len(states)))
        ax.set_xticklabels([s.split("_")[0] for s in states]); ax.set_yticklabels([s.split("_")[0] for s in states])
        for i in range(len(states)):
            for j in range(len(states)):
                ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", color="w" if m.values[i, j] < 0.6 else "k", fontsize=7)
        ax.set_title(arm, fontsize=9)
    fig.tight_layout(); fig.savefig(F / "fig3_state_correlation.png"); plt.close(fig)

    # 4 null floor
    d = sc.merge(lab[["id", "role"]], on="id")
    d["grp"] = np.where(d.role == "random_null", "null", "Env-directed")
    fig, axs = plt.subplots(1, len(states), figsize=(8, 2.8), sharey=True)
    for ax, s in zip(np.atleast_1d(axs), states):
        g = d[(d.state == s) & (d.arm == "glycosylated")]
        for k, (grp, c) in enumerate([("null", "0.6"), ("Env-directed", "#009E73")]):
            v = g[g.grp == grp].score
            ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, medianprops=dict(color="k"))
            ax.scatter(np.random.default_rng(1).normal(k, 0.06, len(v)), v, s=6, color=c, alpha=0.7)
        for role, m in [("control:temsavir_83J", "*"), ("control:BMS378806_83G", "D")]:
            v = g[g.role == role].score
            ax.scatter([1.28], v, marker=m, color="#D55E00", s=40, zorder=4)
        n_out = int((g.score > 15).sum())
        ax.set_ylim(-15, 15)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["null", "Env"]); ax.set_title(f"{s.split('_')[0]} ({n_out} above +15 not shown)", fontsize=8)
    np.atleast_1d(axs)[0].set_ylabel("Vina score (Vina units)")
    fig.tight_layout(); fig.savefig(F / "fig4_null_floor.png"); plt.close(fig)

    # 5 glycan effect
    fig, ax = plt.subplots(figsize=(6, 3))
    for s in states:
        v = gl[(gl.state == s)].delta_glycosylated_minus_free
        ax.hist(v, bins=np.arange(-3, 3.1, 0.25), histtype="step", label=s.split("_")[0], lw=1.2)
    ax.axvline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("Vina score, glycosylated minus glycan-free (positive = worse with glycans)")
    ax.set_ylabel("compounds"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(F / "fig5_glycan_effect.png"); plt.close(fig)


if __name__ == "__main__":
    main()
