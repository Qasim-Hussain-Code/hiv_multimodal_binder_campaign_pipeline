#!/usr/bin/env python3
"""Stage 02. Fetch the Env structures, the ChEMBL Env-annotated activity export
and a set of UniProt Env sequences, each stamped with the date it was pulled.

Nothing here is taken on trust from the project brief: resolution, deposition
date, method and ligand identity are read from RCSB and from the coordinate
file itself. Clade is not a field in the RCSB entry, so it is written to
config/panel.tsv by hand in stage 03 with the source of each assignment.
"""
import argparse, gzip, json, sys, time, datetime, hashlib
from pathlib import Path
import requests, gemmi

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PANEL_IDS = ["4TVP", "5U7O", "5U7M", "6U0L", "8Z7N"]
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"
ENV_TARGETS = ["CHEMBL1293311", "CHEMBL2396504", "CHEMBL3520", "CHEMBL4105927",
               "CHEMBL5057", "CHEMBL5826", "CHEMBL6180", "CHEMBL6181"]
WATER_LIKE = {"HOH", "DOD"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fetch_structures(log):
    d = DATA / "structures"; d.mkdir(parents=True, exist_ok=True)
    meta = {}
    for pid in PANEL_IDS:
        f = d / f"{pid}.cif.gz"
        if not f.exists():
            r = requests.get(f"https://files.rcsb.org/download/{pid}.cif.gz", timeout=120)
            r.raise_for_status(); f.write_bytes(r.content)
        e = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=60).json()
        st = gemmi.read_structure(str(f))
        ligs = sorted({r.name for ch in st[0] for r in ch
                       if r.het_flag == "H" and r.name not in WATER_LIKE and r.entity_type != gemmi.EntityType.Branched})
        meta[pid] = dict(
            title=e["struct"]["title"],
            resolution_A=e["rcsb_entry_info"]["resolution_combined"][0],
            deposited=e["rcsb_accession_info"]["deposit_date"][:10],
            method=e["exptl"][0]["method"],
            non_glycan_ligands=ligs,
            sha256=sha(f), mb=round(f.stat().st_size / 1e6, 2))
        log(f"{pid}: {meta[pid]['resolution_A']} A, {meta[pid]['deposited']}, ligands {ligs}")
    return meta


def fetch_chembl(log):
    out = DATA / "chembl_env_raw.tsv"
    import pandas as pd
    rows = []
    for t in ENV_TARGETS:
        url, p = CHEMBL + "activity.json", {"target_chembl_id": t, "limit": 1000, "standard_units": "nM"}
        while url:
            r = requests.get(url, params=p, timeout=120).json(); p = None
            for a in r["activities"]:
                if (a["standard_type"] in ("IC50", "EC50", "Ki", "Kd") and a["standard_value"]
                        and a["canonical_smiles"] and a["standard_relation"] == "="):
                    rows.append(dict(target=t, mol=a["molecule_chembl_id"], smiles=a["canonical_smiles"],
                                     type=a["standard_type"], nM=float(a["standard_value"]),
                                     assay=a["assay_chembl_id"], doc=a["document_chembl_id"]))
            url = "https://www.ebi.ac.uk" + r["page_meta"]["next"] if r["page_meta"]["next"] else None
    df = pd.DataFrame(rows)
    docs = sorted(df.doc.unique())
    info = {}
    for i in range(0, len(docs), 50):
        r = requests.get(CHEMBL + "document.json", params={"document_chembl_id__in": ",".join(docs[i:i+50]), "limit": 100}, timeout=60).json()
        for x in r["documents"]:
            info[x["document_chembl_id"]] = ((x["title"] or "").replace("\t", " "), x["year"])
    df["doc_title"] = df.doc.map(lambda k: info.get(k, ("", None))[0])
    df["doc_year"] = df.doc.map(lambda k: info.get(k, ("", None))[1])
    df.to_csv(out, sep="\t", index=False)
    v = requests.get(CHEMBL + "status.json", timeout=30).json()
    log(f"ChEMBL {v.get('chembl_db_version')}: {len(df)} activity rows, {df.mol.nunique()} molecules")
    return dict(chembl_version=v.get("chembl_db_version"), rows=len(df), molecules=int(df.mol.nunique()))


def fetch_uniprot(log):
    # Unreviewed and reviewed Env entries, length-filtered to full gp160.
    # LANL alignments need a web form and cannot be scripted reliably, so this
    # is the stand-in and stage 10 says so in its output.
    out = DATA / "uniprot_env.fasta"
    q = "(gene:env) AND (organism_id:11676) AND (length:[600 TO 900]) AND (fragment:false)"
    r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                     params={"query": q, "format": "fasta", "size": 500}, timeout=300)
    r.raise_for_status()
    txt = r.text
    n = txt.count(">")
    out.write_text(txt)
    log(f"UniProt: {n} env sequences (length 600-900, non-fragment)")
    return dict(uniprot_query=q, sequences=n)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--skip-uniprot", action="store_true")
    a = ap.parse_args()
    DATA.mkdir(exist_ok=True)
    (ROOT / "logs").mkdir(exist_ok=True)
    log = lambda m: print(f"[fetch] {m}", file=sys.stderr)
    t0 = time.time()
    rec = dict(date=datetime.date.today().isoformat())
    rec["structures"] = fetch_structures(log)
    rec["chembl"] = fetch_chembl(log)
    if not a.skip_uniprot:
        rec["uniprot"] = fetch_uniprot(log)
    rec["elapsed_s"] = round(time.time() - t0, 1)
    (ROOT / "logs" / "02_fetch.json").write_text(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
