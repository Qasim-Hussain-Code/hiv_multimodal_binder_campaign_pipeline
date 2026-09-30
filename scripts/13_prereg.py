#!/usr/bin/env python3
"""Stage 13. Write and hash the pre-registration.

Run this after stages 11 and 12 and before any assay result exists. It writes
results/order_list.tsv and results/preregistration.md, records the git commit
the tree was at, the SHA-256 of every input table and the date, and refuses to
overwrite an existing preregistration. The file is committed alone with the
message add_preregistration and is not amended. A correction goes in a new file
that names this one and its commit hash.

Order list rule (within the small-molecule modality, no pooled score):
  1. glycosylated arm only, since that is the arm nearer a real virion;
  2. keep candidates that are not escape-prone in at least 3 of the 4 states;
  3. sort by the median unfiltered rank percentile across the 4 states;
  4. take the top TOP_N.
Median across states is a rank aggregation inside one modality and is not a
cross-modality score. Two anchors (temsavir, BMS-378806) and N_NULL random
property-matched null molecules are added as assay controls: the anchors
check that the assay can see a known inhibitor, the nulls give the assay's own
false-positive floor. The states in the median are chosen by availability of
structures, and TOP_N and the 3-of-4 rule are arbitrary.
"""
import argparse, datetime, hashlib, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
TOP_N, N_NULL = 10, 5
INPUTS = ["results/dock_scores.tsv", "results/candidate_labels.tsv", "results/candidates_anonymous.tsv",
          "results/conservation_contacts.tsv", "results/conservation_entropy.tsv", "results/ranks_long.tsv",
          "results/control_ranks.tsv", "results/null_floor.tsv", "results/funnel.tsv",
          "results/state_rank_correlation.tsv", "results/glycan_effect.tsv", "results/potency_association.tsv",
          "results/redock_temsavir_S2.tsv", "config/panel.tsv", "config/boxes.json", "config/candidate_sources.tsv",
          "config/cost_assumptions.tsv"]


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True).stdout.strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--seed", type=int, default=20260929)
    a = ap.parse_args()
    out = R / "preregistration.md"
    if out.exists():
        sys.exit(f"{out} exists. The preregistration is never overwritten. Write a correction file that references it.")
    L = pd.read_csv(R / "ranks_long.tsv", sep="\t")
    states = sorted(L.state.unique())
    g = L[L.arm == "glycosylated"]
    cand = g[g.role == "candidate"]
    piv = cand.pivot(index="id", columns="state", values="pctl_unfiltered")
    esc = cand.pivot(index="id", columns="state", values="escape_prone").astype(float)
    ok = (esc.fillna(1).sum(axis=1) <= 1)
    sel = piv[ok].median(axis=1).sort_values().head(TOP_N)
    info = cand.drop_duplicates("id").set_index("id")
    lab = pd.read_csv(ROOT / "data" / "candidates" / "labels.tsv", sep="\t").set_index("id")
    smi = pd.read_csv(R / "candidates_anonymous.tsv", sep="\t").set_index("id").smiles
    rows = []
    for i, (cid, med) in enumerate(sel.items(), 1):
        rows.append(dict(order=i, id=cid, kind="ranked candidate", median_pctl_across_states=round(float(med), 1),
                         chembl=info.loc[cid, "chembl"], smiles=smi[cid]))
    for role in ("control:temsavir_83J", "control:BMS378806_83G"):
        cid = lab[lab.role == role].index[0]
        rows.append(dict(order="", id=cid, kind="anchor " + role.split(":")[1], median_pctl_across_states=round(
            float(g[g.id == cid].pivot(index="id", columns="state", values="pctl_unfiltered").median(axis=1).iloc[0]), 1),
            chembl=lab.loc[cid, "chembl"], smiles=smi[cid]))
    nulls = lab[lab.role == "random_null"].sample(N_NULL, random_state=a.seed)
    for cid, r in nulls.iterrows():
        rows.append(dict(order="", id=cid, kind="null negative", median_pctl_across_states="", chembl=r.chembl, smiles=smi[cid]))
    ol = pd.DataFrame(rows)
    ol.to_csv(R / "order_list.tsv", sep="\t", index=False)

    cost = pd.read_csv(ROOT / "config" / "cost_assumptions.tsv", sep="\t")
    n_tests = len(ol)
    hashes = "\n".join(f"| {p} | `{sha(p)}` |" for p in INPUTS)
    dirty = git("status", "--porcelain")
    table = "\n".join(f"| {r['order']} | {r['id']} | {r['kind']} | {r['median_pctl_across_states']} | {r['chembl']} |" for r in rows)
    text = f"""# Pre-registration: small-molecule arm, HIV-1 Env, tier A

Written {datetime.date.today().isoformat()} by scripts/13_prereg.py, before any assay result exists.
Tree commit: `{git('rev-parse', 'HEAD')}`. Uncommitted changes at the moment of writing: {len(dirty.splitlines())} paths.
Bench tier: A (no wet lab). Nothing below has been tested at a bench.

Scope. Only the small-molecule modality was ranked. The peptide, designed-binder and aptamer arms were not
run (results/modality_*_not_run.md), so this document makes no prediction about them and no cross-modality
statement is made.

## What is predicted

The {TOP_N} ranked candidates below are the compounds this ranking would order first. The ranking is a Vina score
(AutoDock Vina 1.2.7, empirical scoring function, ranking key only) into the temsavir pocket of four deposited Env
states, glycosylated arm, median unfiltered rank percentile across states, after removing candidates that are
escape-prone in two or more states. The score is not an affinity and not a potency estimate.

| order | id | kind | median percentile across states | ChEMBL id |
|---|---|---|---|---|
{table}

Ordering {n_tests} compounds in total: {TOP_N} ranked candidates, 2 anchors, {N_NULL} null negatives.
SMILES are in results/order_list.tsv.

## Assay and readout

Env-pseudotyped single-round neutralisation in TZM-bl cells (Sarzotti-Kelsoe et al., J Immunol Methods 2014,
409:131-146), luciferase readout, IC50 by four-parameter fit, duplicate wells, on BG505 and at least three further
viruses from the global panel (deCamp et al., J Virol 2014) chosen to include one non-clade-A virus. Containment
requirements must be confirmed with an institutional biosafety committee before this is run. This repository does
not supply a classification. No approval reference exists at tier A.

## Result that would support the ranking

All of the following, stated before any measurement:
1. The temsavir anchor gives a curve with IC50 below 1 uM on BG505 (assay validity; if this fails the assay, not the
   ranking, is in doubt and no other criterion is evaluated).
2. At least 3 of the {TOP_N} ranked candidates reach IC50 below 10 uM on BG505 at 50 percent neutralisation with a
   full dose-response curve.
3. No more than 1 of the {N_NULL} null negatives reaches the same threshold.
The 3 of {TOP_N}, the 10 uM line and the 1 uM anchor line are arbitrary and were fixed by eye.

## Result that would refute the ranking

Any of the following, with the anchor valid:
1. 0 of the {TOP_N} ranked candidates reach IC50 below 10 uM on BG505.
2. The hit rate among the ranked candidates is not higher than the hit rate among the null negatives.
3. Every ranked candidate that neutralises BG505 fails to neutralise every panel virus from a different clade at
   the same threshold (potency without breadth).
A refutation on criterion 1 or 2 means Vina scoring into this pocket, in these four states, does not enrich for
neutralising compounds. It would not say anything about temsavir or about the pocket.

## Cost and bench time

Stated in results/cost_model.tsv with its assumptions (config/cost_assumptions.tsv). Those values are placeholders
set by the author, not vendor quotes.

## Inputs, by SHA-256

| file | sha256 |
|---|---|
{hashes}

## Corrections

None. A correction is a new file that names this one and the commit above.
"""
    out.write_text(text)
    print(f"[prereg] wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
