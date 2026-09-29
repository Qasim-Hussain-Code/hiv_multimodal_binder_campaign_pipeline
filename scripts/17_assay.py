#!/usr/bin/env python3
"""Stage 17. Tier B or C only. Ingest assay results and compare them with the
pre-registration.

At tier A this script refuses. No assay was run for this repository and no
file named on the command line is read at tier A.

At tier B or C it expects a TSV with columns id, virus, ic50_uM (blank if no
curve), curve_ok (1 or 0). It evaluates the criteria in results/preregistration.md
exactly as written there: no criterion is edited here after the fact.
"""
import argparse, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", help="assay results TSV")
    ap.add_argument("--tier", default=None)
    a = ap.parse_args()
    conf = {}
    p = ROOT / "project.conf"
    if p.exists():
        for l in p.read_text().splitlines():
            if "=" in l and not l.startswith("#"):
                k, v = l.split("=", 1); conf[k] = v.strip('"')
    tier = a.tier or conf.get("BENCH_TIER", "A")
    if tier == "A":
        sys.exit("tier A: no wet lab. There is nothing to ingest. This is the intended end of the repository at tier A.")
    if not a.results:
        sys.exit("pass --results")
    cfg = (ROOT / "config" / "bench_tier.conf").read_text()
    if 'BIOSAFETY_APPROVAL_REF=""' in cfg:
        sys.exit("BIOSAFETY_APPROVAL_REF is empty in config/bench_tier.conf. Record the approval before ingesting results.")
    ol = pd.read_csv(ROOT / "results" / "order_list.tsv", sep="\t")
    r = pd.read_csv(a.results, sep="\t")
    bg = r[r.virus.str.contains("BG505", case=False)].merge(ol[["id", "kind"]], on="id")
    bg["hit"] = (bg.ic50_uM < 10) & (bg.curve_ok == 1)
    anchor = bg[bg.kind == "anchor temsavir_83J"]
    valid = bool(len(anchor) and (anchor.ic50_uM < 1).any())
    top = bg[bg.kind == "ranked candidate"]; nul = bg[bg.kind == "null negative"]
    print(f"anchor valid (temsavir IC50 < 1 uM on BG505): {valid}")
    if not valid:
        sys.exit("assay validity criterion failed. No other criterion is evaluated.")
    print(f"ranked hits: {int(top.hit.sum())} of {len(top)}; null hits: {int(nul.hit.sum())} of {len(nul)}")
    print("supports:", bool(top.hit.sum() >= 3 and nul.hit.sum() <= 1))
    print("refutes (criterion 1 or 2):", bool(top.hit.sum() == 0 or top.hit.mean() <= (nul.hit.mean() if len(nul) else 0)))


if __name__ == "__main__":
    main()
