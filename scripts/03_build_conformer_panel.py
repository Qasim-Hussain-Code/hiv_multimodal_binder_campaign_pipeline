#!/usr/bin/env python3
"""Stage 03. Assemble the conformational panel and place one docking box per
state at the temsavir pocket.

Panel logic. Every state is a deposited coordinate set. No State 1 structure
exists (Lu et al. 2019 placed the solved structures in States 2 and 3), so the
panel cannot contain the conformation the virion most often presents. That is
a limitation of the data and this script does not try to fix it.

Box placement. The box centre is the centroid of temsavir (CCD 83J) in 5U7O,
carried into every other frame by a C-alpha superposition of the gp120
protomer that bears the pocket. The alternative was to define the pocket by
residue number, which breaks for the CH119 strain (8Z7N) whose numbering and
sequence differ from BG505. Superposition needs only that the fold aligns.

The pocket is induced. In the closed, unliganded trimer (4TVP) the cavity that
temsavir occupies in 5U7O is not open to the same degree, so docking into it
is the case stage 1 taught to distrust.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, gemmi

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# state, pdb, gp120 entity id (as gemmi numbers them), env entities kept,
# clade, clade source, conformation, CD4 bound
PANEL = [
 ("S1_closed_apo",   "4TVP", "1", ["1", "2"], "A", "BG505, clade A (Pancera 2014)", "closed, pocket not liganded", "no"),
 ("S2_closed_holo",  "5U7O", "4", ["1", "4"], "A", "BG505, clade A (Pancera 2017)", "closed, temsavir-bound (83J)", "no"),
 ("S3_open_cd4",     "6U0L", "1", ["1", "5"], "A", "BG505, clade A (RCSB title)", "asymmetrically open, CD4 and E51 bound", "yes"),
 ("S4_cd4_cladeCRF", "8Z7N", "1", ["1", "2"], "CRF07_BC?", "strain CH119; subtype NOT read from a database field, from deposition title 'Asia CRFs'. Verify.", "intermediate open, CD4 bound", "yes"),
]
REF_PDB, REF_LIGAND = "5U7O", "83J"


def gp120_polymer(st, ent_id):
    # gemmi entity ids, not chain names: the same protein sits under different
    # chain letters in each deposition.
    ch = next(c for c in st[0] if c[0].entity_id == ent_id and c[0].entity_type == gemmi.EntityType.Polymer)
    return ch.get_polymer()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.parse_args()
    ref = gemmi.read_structure(str(DATA / "structures" / f"{REF_PDB}.cif.gz")); ref.setup_entities()
    ref_poly = gp120_polymer(ref, "4")
    lig = [r for ch in ref[0] for r in ch if r.name == REF_LIGAND]
    pts = [a.pos for r in lig[:1] for a in r]
    ref_centre = gemmi.Position(np.mean([p.x for p in pts]), np.mean([p.y for p in pts]), np.mean([p.z for p in pts]))
    ref_ligand_xyz = [[a.pos.x, a.pos.y, a.pos.z] for a in lig[0]]

    out, boxes = [], {}
    for state, pdb, gp, keep, clade, csrc, conf, cd4 in PANEL:
        st = gemmi.read_structure(str(DATA / "structures" / f"{pdb}.cif.gz")); st.setup_entities()
        poly = gp120_polymer(st, gp)
        sup = gemmi.calculate_superposition(ref_poly, poly, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
        # sup.transform maps `poly` onto `ref_poly`; the box centre is wanted in the poly frame
        inv = sup.transform.inverse()
        v = inv.apply(gemmi.Vec3(ref_centre.x, ref_centre.y, ref_centre.z))
        c = gemmi.Position(v.x, v.y, v.z)
        near = sum(1 for ch in st[0] for r in ch for a in r if a.pos.dist(c) < 8.0 and r.het_flag == "A")
        meta = json.loads((ROOT / "logs" / "02_fetch.json").read_text())["structures"][pdb]
        boxes[state] = dict(pdb=pdb, centre=[round(c.x, 3), round(c.y, 3), round(c.z, 3)],
                            gp120_entity=gp, env_entities=keep, superposition_rmsd=round(sup.rmsd, 2),
                            aligned_ca=sup.count, protein_atoms_within_8A=near)
        out.append([state, pdb, meta["resolution_A"], meta["deposited"], meta["method"], clade, csrc, conf, cd4,
                    ";".join(meta["non_glycan_ligands"]), sup.count, round(sup.rmsd, 2)])
        print(f"[panel] {state} {pdb}: {sup.count} CA aligned, rmsd {sup.rmsd:.2f} A, {near} protein atoms within 8 A of the box centre", file=sys.stderr)
    (ROOT / "config").mkdir(exist_ok=True)
    hdr = ["state", "pdb", "resolution_A", "deposited", "method", "clade", "clade_source", "conformation",
           "cd4_bound", "non_glycan_ligands", "aligned_ca_to_5U7O", "ca_rmsd_to_5U7O"]
    (ROOT / "config" / "panel.tsv").write_text("\n".join(["\t".join(hdr)] + ["\t".join(map(str, r)) for r in out]) + "\n")
    (DATA / "panel").mkdir(exist_ok=True)
    (DATA / "panel" / "boxes.json").write_text(json.dumps(dict(boxes=boxes, reference_ligand_xyz=ref_ligand_xyz), indent=1))
    (ROOT / "config" / "boxes.json").write_text(json.dumps(boxes, indent=1))


if __name__ == "__main__":
    main()
