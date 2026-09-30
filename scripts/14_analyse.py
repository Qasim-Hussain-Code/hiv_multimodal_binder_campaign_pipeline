#!/usr/bin/env python3
"""Stage 14. Funnels, glycan counts, cost model, timing distributions and the
numbers file that the README is written from.

Every figure quoted in the README should be found in results/summary_numbers.json
or in one of the tables this and the earlier stages write. If a number is not
in one of them it does not go in the README.

Cost model. The unit costs in config/cost_assumptions.tsv are placeholders set
by the author, not quotes. Bench hours are reported separately from currency
and labour is not priced.
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"


def q(x, p):
    return round(float(np.percentile(x, p)), 2)


def main():
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args()
    sc = pd.read_csv(R / "dock_scores.tsv", sep="\t")
    lab = pd.read_csv(ROOT / "data" / "candidates" / "labels.tsv", sep="\t")
    ctl = pd.read_csv(R / "control_ranks.tsv", sep="\t")
    fun = pd.read_csv(R / "funnel.tsv", sep="\t")
    gl = pd.read_csv(R / "glycan_effect.tsv", sep="\t")
    con = pd.read_csv(R / "conservation_contacts.tsv", sep="\t")
    ol = pd.read_csv(R / "order_list.tsv", sep="\t")
    cs = json.loads((R / "candidate_summary.json").read_text())

    # timing and memory distributions, per state and arm
    t = sc.groupby(["state", "arm"]).agg(n=("wall_s", "size"), wall_median_s=("wall_s", "median"),
        wall_p05_s=("wall_s", lambda x: q(x, 5)), wall_p95_s=("wall_s", lambda x: q(x, 95)), wall_max_s=("wall_s", "max"),
        rss_median_mb=("peak_rss_mb", "median"), rss_p95_mb=("peak_rss_mb", lambda x: q(x, 95)), rss_max_mb=("peak_rss_mb", "max")).round(1)
    t.to_csv(R / "timing_distributions.tsv", sep="\t")

    # glycan effect per state
    g = gl.merge(lab[["id"]], on="id")
    gs = g[g.role != "control:x"].groupby("state").agg(compounds=("id", "size"), lost_to_glycans=("lost_to_glycans", "sum"),
        median_delta=("delta_glycosylated_minus_free", "median"), max_delta=("delta_glycosylated_minus_free", "max")).round(2)
    gc = con[con.arm == "glycosylated"].groupby("state").agg(poses_touching_glycan=("n_glycan_residues_in_contact", lambda x: int((x > 0).sum())))
    gs = gs.join(gc)
    gs.to_csv(R / "glycan_summary.tsv", sep="\t")

    # cost of the order list
    a = pd.read_csv(ROOT / "config" / "cost_assumptions.tsv", sep="\t").set_index("item")["value"].astype(float)
    n = len(ol)
    v = a["viruses_per_compound"]
    cost = pd.DataFrame([
        dict(item="compounds ordered", value=n, unit="compounds"),
        dict(item="purchase", value=n * a["compound_purchase"], unit="USD"),
        dict(item="assay reagents", value=n * v * a["assay_reagents_per_compound"], unit="USD"),
        dict(item="total currency", value=n * (a["compound_purchase"] + v * a["assay_reagents_per_compound"]), unit="USD"),
        dict(item="bench hours", value=n * v * a["bench_hours_per_compound"], unit="hours"),
        dict(item="cost per compound tested", value=a["compound_purchase"] + v * a["assay_reagents_per_compound"], unit="USD"),
        dict(item="bench hours per compound tested", value=v * a["bench_hours_per_compound"], unit="hours"),
    ])
    cost.to_csv(R / "cost_model.tsv", sep="\t", index=False)

    # stage wall clock and peak memory from run_timed
    st = []
    for f in sorted((ROOT / "logs").glob("*.timing")):
        parts = f.read_text().split()
        st.append(dict(stage=f.stem, wall_s=parts[0], peak_rss_kb=parts[1] if len(parts) > 1 else "NA"))
    pd.DataFrame(st).to_csv(R / "stage_timing.tsv", sep="\t", index=False)

    dsz = sum(p.stat().st_size for p in (ROOT / "data").rglob("*") if p.is_file()) / 1e6
    summ = dict(candidate_set=cs, dockings=len(sc), disk_data_mb=round(dsz, 1),
                controls=ctl.to_dict("records"), funnel=fun.to_dict("records"),
                total_dock_wall_h=round(float(sc.wall_s.sum()) / 3600, 2),
                n_ordered=n)
    (R / "summary_numbers.json").write_text(json.dumps(summ, indent=1, default=str))
    print(t.to_string()); print(gs.to_string()); print(cost.to_string(index=False))


if __name__ == "__main__":
    main()
