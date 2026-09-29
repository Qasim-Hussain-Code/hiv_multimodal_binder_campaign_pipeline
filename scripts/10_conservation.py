#!/usr/bin/env python3
"""Stage 10. Per-position entropy, contact residues per pose, escape-prone flags.

Alignment source. The brief asks for a Los Alamos HIV sequence database
alignment. LANL alignments come from a web form and there is no scriptable
download, so this stage uses 3,000+ UniProt HIV-1 env sequences (non-fragment,
600-900 residues; the query is in logs/02_fetch.json), each aligned pairwise
(local, BLOSUM62) to the BG505 reference and read off column by column. That
is weaker than a curated multiple alignment: pairwise alignment to one
reference mis-places residues inside indel-rich loops (V1, V2), which inflates
entropy exactly where entropy is already high. The V1/V2/V3 flag does not
depend on the alignment, only on HXB2 position ranges. UniProt is also a
convenience sample of what was deposited, not of the global epidemic, and its
clade composition was not measured here.

Reference numbering is the author numbering of 5U7O (HXB2-based). Other states
are mapped to it by aligning their gp120 and gp41 chains, because CH119 (8Z7N)
does not share BG505 numbering.

Thresholds are arbitrary and named as such: a contact is 'hypervariable' if
its column entropy exceeds the 75th percentile of all reference positions, and
a pose is escape-prone if any contact lies in V1 (131-157), V2 (158-196) or V3
(296-331), or if at least one third of contacts are hypervariable. The one
third and the 75th percentile were chosen by eye before looking at any
candidate. Contacts are receptor residues with a heavy atom within 4.5 A of a
ligand heavy atom in the top-ranked Vina pose.
"""
import argparse, json, math, sys
from collections import Counter
from pathlib import Path
import numpy as np, pandas as pd, gemmi
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
VREG = {"V1": (131, 157), "V2": (158, 196), "V3": (296, 331)}
CONTACT_A = 4.5
HYPER_PCTL = 75
HYPER_FRAC = 1 / 3


def poly_seq(chain):
    p = chain.get_polymer()
    out = []
    for r in p:
        aa = gemmi.find_tabulated_residue(r.name)
        out.append((r.seqid.num, r.seqid.icode.strip(), aa.one_letter_code.upper() if aa else "X"))
    return out


def aligner():
    al = PairwiseAligner(); al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10; al.extend_gap_score = -0.5
    return al


def pair_map(al, ref_seq, q_seq):
    """map index in q_seq -> index in ref_seq via the best local alignment"""
    aln = al.align(ref_seq, q_seq)[0]
    m = {}
    for (r0, r1), (q0, q1) in zip(*aln.aligned):
        for k in range(r1 - r0):
            m[q0 + k] = r0 + k
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scores", default=str(ROOT / "results" / "dock_scores.tsv"))
    a = ap.parse_args()
    al = aligner()

    # reference: 5U7O gp120 (entity 4) then gp41 (entity 1), modelled residues only
    ref = gemmi.read_structure(str(DATA / "structures" / "5U7O.cif.gz")); ref.setup_entities()
    def chain_of(st, ent):
        return next(c for c in st[0] if c[0].entity_id == ent and c[0].entity_type == gemmi.EntityType.Polymer)
    rs = poly_seq(chain_of(ref, "4")) + poly_seq(chain_of(ref, "1"))
    ref_str = "".join(x[2] for x in rs)
    ref_label = [(n, ic) for n, ic, _ in rs]

    # column entropy from pairwise alignments
    cols = [Counter() for _ in rs]
    n_seq = 0
    for rec in SeqIO.parse(DATA / "uniprot_env.fasta", "fasta"):
        s = str(rec.seq).replace("*", "X")
        m = pair_map(al, ref_str, s)
        for qi, ri in m.items():
            if s[qi] != "X":
                cols[ri][s[qi]] += 1
        n_seq += 1
    ent, cov = [], []
    for c in cols:
        n = sum(c.values()); cov.append(n / n_seq)
        ent.append(-sum(v / n * math.log2(v / n) for v in c.values()) if n >= 0.5 * n_seq else float("nan"))
    ent = np.array(ent)
    thr = float(np.nanpercentile(ent, HYPER_PCTL))
    df = pd.DataFrame(dict(hxb2=[n for n, _ in ref_label], icode=[i for _, i in ref_label], aa=list(ref_str),
                           entropy_bits=ent, coverage=cov))
    df.to_csv(ROOT / "results" / "conservation_entropy.tsv", sep="\t", index=False)
    (ROOT / "results" / "conservation_summary.json").write_text(json.dumps(
        dict(sequences_aligned=n_seq, positions=len(rs), positions_with_entropy=int(np.isfinite(ent).sum()),
             hypervariable_threshold_bits=round(thr, 3), hyper_percentile=HYPER_PCTL, hyper_fraction_rule=HYPER_FRAC,
             median_entropy_bits=round(float(np.nanmedian(ent)), 3),
             alignment="UniProt env sequences aligned pairwise to BG505 (not a LANL alignment)"), indent=1))
    print(f"[conservation] {n_seq} sequences, threshold {thr:.2f} bits", file=sys.stderr)

    # per-state numbering map onto the reference
    boxes = json.loads((ROOT / "config" / "boxes.json").read_text())
    maps = {}
    for state, b in boxes.items():
        st = gemmi.read_structure(str(DATA / "structures" / f"{b['pdb']}.cif.gz")); st.setup_entities()
        ents = b["env_entities"]
        q = poly_seq(chain_of(st, ents[0])) + poly_seq(chain_of(st, ents[1]))
        q_str = "".join(x[2] for x in q)
        pm = pair_map(al, ref_str, q_str)
        maps[state] = {(q[i][0], q[i][1]): pm[i] for i in pm}

    # contacts
    if not Path(a.scores).exists():
        sys.exit("no dock scores yet")
    sc = pd.read_csv(a.scores, sep="\t")
    rows = []
    rec_cache = {}
    for state, arm in sc[["state", "arm"]].drop_duplicates().itertuples(index=False):
        st = gemmi.read_structure(str(DATA / "receptors" / f"{state}__{arm}.pdb"))
        xyz, res = [], []
        for ch in st[0]:
            for r in ch:
                for at in r:
                    if at.element.name != "H":
                        xyz.append([at.pos.x, at.pos.y, at.pos.z])
                        res.append((r.name, r.seqid.num, r.seqid.icode.strip(), r.het_flag))
        rec_cache[(state, arm)] = (cKDTree(np.array(xyz)), res)
    for r in sc.itertuples(index=False):
        pose = DATA / "poses" / f"{r.id}__{r.state}__{r.arm}.pdbqt"
        lx = []
        for l in pose.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")) and not l[77:79].strip().startswith("H"):
                lx.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            if l.startswith("ENDMDL"):
                break
        tree, res = rec_cache[(r.state, r.arm)]
        idx = set()
        for hits in tree.query_ball_point(np.array(lx), CONTACT_A):
            idx.update(hits)
        prot = {(res[i][1], res[i][2]) for i in idx if res[i][3] == "A"}
        n_gly = len({(res[i][1], res[i][2], res[i][0]) for i in idx if res[i][3] != "A"})
        mp = maps[r.state]
        ents_ = [ent[mp[k]] for k in prot if k in mp and np.isfinite(ent[mp[k]])]
        pos = [ref_label[mp[k]][0] for k in prot if k in mp]
        inv = {n: any(lo <= p <= hi for p in pos) for n, (lo, hi) in VREG.items()}
        fh = float(np.mean([e > thr for e in ents_])) if ents_ else float("nan")
        esc = bool(any(inv.values()) or (ents_ and fh >= HYPER_FRAC))
        rows.append(dict(id=r.id, state=r.state, arm=r.arm, n_contacts=len(prot), n_mapped=len(ents_),
                         mean_entropy_bits=float(np.mean(ents_)) if ents_ else float("nan"),
                         frac_hypervariable=fh, in_V1=inv["V1"], in_V2=inv["V2"], in_V3=inv["V3"],
                         escape_prone=esc, n_glycan_residues_in_contact=n_gly))
    pd.DataFrame(rows).to_csv(ROOT / "results" / "conservation_contacts.tsv", sep="\t", index=False)
    print(f"[conservation] {len(rows)} poses scored, {sum(x['escape_prone'] for x in rows)} escape-prone", file=sys.stderr)


if __name__ == "__main__":
    main()
