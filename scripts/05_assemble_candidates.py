#!/usr/bin/env python3
"""Stage 05. Assemble the small-molecule candidate set, the two positive
controls and the null set, and write manifests.

Blinding. The docking stage reads only data/candidates/candidates.tsv, which
carries a shuffled anonymous id and a SMILES. Which ids are controls, which are
null draws and which ChEMBL potency belongs to which id lives in labels.tsv,
read only by stage 11. The control is scored by the same function as every
other row for that reason, and the docking code has no way to tell it apart.

Sources. The brief names stage 4's LIT-PCBA-calibrated screen as the source of
candidates. That screen is calibrated on 15 non-Env targets and holds no Env
compounds, so it cannot supply them. The candidate set is ChEMBL compounds
annotated to an HIV-1 Env target whose source paper is on the include list in
config/candidate_sources.tsv. ChEMBL's Env annotation is noisy (NNRTI and
integrase compounds sit under gp160 targets), which is why the include list is
by paper and not by target id.

Null. Molecules drawn at random from ChEMBL (300 <= MW <= 600, no more than one
Lipinski violation, InChIKey absent from the Env-annotated export), then
sub-sampled so the heavy-atom histogram matches the candidates in 5-atom bins.
Matching matters because Vina scores grow with heavy-atom count, so an
unmatched null would look worse than the candidates for no biological reason.
"""
import argparse, json, random, sys
from pathlib import Path
import pandas as pd, requests
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"
HEAVY_MAX = 55
PROPS = {"molecule_properties__mw_freebase__range": "300,600", "molecule_properties__num_ro5_violations__lte": 1}


def ikey14(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToInchiKey(m)[:14] if m else None


def ccd_smiles(cid):
    d = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{cid}", timeout=60).json()
    return d["rcsb_chem_comp_descriptor"]["SMILES_stereo"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--n-null", type=int, default=100)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = DATA / "candidates"; out.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(DATA / "chembl_env_raw.tsv", sep="\t")
    src = pd.read_csv(ROOT / "config" / "candidate_sources.tsv", sep="\t")
    inc = set(src[src.decision == "include"].doc)
    n0 = raw.mol.nunique()
    mols = raw.smiles.map(Chem.MolFromSmiles)
    raw["mw"] = [Descriptors.MolWt(m) if m else None for m in mols]
    raw["heavy"] = [m.GetNumHeavyAtoms() if m else None for m in mols]

    excl = []
    in_inc = set(raw[raw.doc.isin(inc)].mol)
    for mol in sorted(set(raw.mol) - in_inc):
        excl.append((mol, "source paper not on include list"))
    kept = raw[raw.doc.isin(inc)].copy()
    big = set(kept[(kept.heavy > HEAVY_MAX) | (kept.mw > 700)].mol)
    excl += [(m, f"heavy atoms > {HEAVY_MAX} or MW > 700") for m in sorted(big)]
    kept = kept[~kept.mol.isin(big)]
    kept["ik"] = kept.smiles.map(ikey14)
    # One row per structure; potency is the most potent reported value, a
    # minimum over heterogeneous assays (cell type, virus, readout differ).
    kept = kept.sort_values("nM").drop_duplicates("ik")

    controls = {"temsavir_83J": ccd_smiles("83J"), "BMS378806_83G": ccd_smiles("83G")}
    cik = {k: ikey14(v) for k, v in controls.items()}
    rows = [dict(smiles=r.smiles, role="candidate", chembl=r.mol, nM=r.nM, assay_type=r.type, doc=r.doc, ik=r.ik)
            for r in kept.itertuples()]
    by_ik = {r["ik"]: r for r in rows}
    in_chembl = {}
    for k, s in controls.items():
        in_chembl[k] = cik[k] in by_ik
        if cik[k] in by_ik:
            by_ik[cik[k]]["role"] = "control:" + k
            by_ik[cik[k]]["smiles"] = s      # coordinates-derived SMILES from the CCD, same molecule
        else:
            rows.append(dict(smiles=s, role="control:" + k, chembl="", nM=float("nan"), assay_type="", doc="CCD", ik=cik[k]))

    known = set(raw.smiles.map(ikey14).dropna()) | {r["ik"] for r in rows}
    total = requests.get(CHEMBL + "molecule.json", params={**PROPS, "limit": 1}, timeout=60).json()["page_meta"]["total_count"]
    pool = []
    while len(pool) < 3000:
        off = rng.randrange(0, total - 100)
        r = requests.get(CHEMBL + "molecule.json", params={**PROPS, "limit": 100, "offset": off}, timeout=120).json()
        for m in r["molecules"]:
            st = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if st and "." not in st:
                mm = Chem.MolFromSmiles(st)
                if mm and mm.GetNumHeavyAtoms() <= HEAVY_MAX and ikey14(st) not in known:
                    pool.append((m["molecule_chembl_id"], st, mm.GetNumHeavyAtoms()))
    b = lambda h: h // 5
    want = {}
    for r in rows:
        k = b(Chem.MolFromSmiles(r["smiles"]).GetNumHeavyAtoms()); want[k] = want.get(k, 0) + 1
    tot = sum(want.values()); rng.shuffle(pool); nulls = []
    for k, n in sorted(want.items()):
        nulls += [p for p in pool if b(p[2]) == k][:round(a.n_null * n / tot)]
    for cid, st, h in nulls:
        rows.append(dict(smiles=st, role="null", chembl=cid, nM=float("nan"), assay_type="", doc="ChEMBL_random", ik=ikey14(st)))

    rng.shuffle(rows)
    for i, r in enumerate(rows):
        r["id"] = f"C{i+1:04d}"
    df = pd.DataFrame(rows)
    (ROOT / "results").mkdir(exist_ok=True)
    df[["id", "smiles"]].to_csv(out / "candidates.tsv", sep="\t", index=False)
    df[["id", "role", "chembl", "nM", "assay_type", "doc", "ik"]].to_csv(out / "labels.tsv", sep="\t", index=False)
    df[["id", "smiles"]].to_csv(ROOT / "results" / "candidates_anonymous.tsv", sep="\t", index=False)
    df[["id", "role", "chembl", "nM", "assay_type", "doc"]].to_csv(ROOT / "results" / "candidate_labels.tsv", sep="\t", index=False)
    pd.DataFrame(excl, columns=["chembl", "reason"]).to_csv(ROOT / "results" / "excluded_candidates.tsv", sep="\t", index=False)
    summ = dict(env_annotated_molecules=int(n0), excluded=len({m for m, _ in excl}),
                candidates=int((df.role == "candidate").sum()), controls=int(df.role.str.startswith("control").sum()),
                null=int((df.role == "null").sum()), null_pool=len(pool), seed=a.seed, control_found_in_chembl_env=in_chembl)
    (ROOT / "results" / "candidate_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ), file=sys.stderr)


if __name__ == "__main__":
    main()
