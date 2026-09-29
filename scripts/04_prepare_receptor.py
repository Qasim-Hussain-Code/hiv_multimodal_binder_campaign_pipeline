#!/usr/bin/env python3
"""Stage 04. Build one cropped receptor per state and per glycan arm.

Glycan arms. The brief asks for modelled glycans. No glycan modelling tool
was run: the machine had 2 GB of free disk and no GlycoSHIELD or Rosetta
glycan library, and a defensible model of the shield needs one or the other.
What this stage does instead is keep the sugar residues that the deposition
itself resolved (the 'glycosylated' arm) or strip them (the 'glycan_free'
arm). Resolved glycans are a fraction of the real shield, so the glycosylated
arm is a lower bound on steric occlusion, and the README says so.

Receptor content. Env chains only (gp120 and gp41 of all three protomers, via
the biological assembly). CD4, Fabs and other partners are removed in every
state, including the CD4-bound ones, because the question is what a small
molecule sees on the trimer and not what a CD4-bound complex excludes. The
alternative, keeping CD4 as an obstacle, would penalise every compound in
S3 and S4 for a partner it would not meet on a free virion.

Cropping. Residues with any atom within CROP_A of the box centre are kept
whole. The crop bounds pdbqt size and Vina start-up time. It cuts the chain,
so termini are artefacts of the crop and sit at least CROP_A minus half the
box away from the search space.
"""
import argparse, json, subprocess, sys, shutil
from pathlib import Path
import gemmi

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CROP_A = 24.0
OBABEL = shutil.which("obabel") or "obabel"


def build(pdb, box, arm):
    st = gemmi.read_structure(str(DATA / "structures" / f"{pdb}.cif.gz")); st.setup_entities()
    m = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    c = gemmi.Position(*box["centre"])
    keep = set(box["env_entities"])
    out = gemmi.Structure(); out.spacegroup_hm = "P 1"; out.cell = gemmi.UnitCell(1, 1, 1, 90, 90, 90)
    om = gemmi.Model("1")
    letters = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    n_prot = n_gly = 0
    for ch in m:
        if len(ch) == 0:
            continue
        r0 = ch[0]
        is_poly = r0.entity_type == gemmi.EntityType.Polymer and r0.entity_id in keep
        is_gly = r0.entity_type == gemmi.EntityType.Branched
        if not (is_poly or (is_gly and arm == "glycosylated")):
            continue
        new = gemmi.Chain(next(letters))
        for r in ch:
            if not any(a.pos.dist(c) < CROP_A for a in r):
                continue
            r2 = gemmi.Residue(); r2.name = r.name; r2.seqid = r.seqid; r2.het_flag = r.het_flag
            for a in r:
                if a.altloc not in ("\0", "A"):
                    continue
                a2 = gemmi.Atom(); a2.name = a.name; a2.element = a.element; a2.pos = a.pos
                a2.occ = 1.0; a2.b_iso = a.b_iso
                r2.add_atom(a2)
            if len(r2):
                new.add_residue(r2)
                if is_poly: n_prot += 1
                else: n_gly += 1
        if len(new):
            om.add_chain(new)
    out.add_model(om)
    return out, n_prot, n_gly


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.parse_args()
    boxes = json.loads((ROOT / "config" / "boxes.json").read_text())
    rd = DATA / "receptors"; rd.mkdir(parents=True, exist_ok=True)
    summary = {}
    for state, box in boxes.items():
        for arm in ("glycan_free", "glycosylated"):
            f = rd / f"{state}__{arm}"
            if f.with_suffix(".pdbqt").exists():
                continue
            st, npro, ngly = build(box["pdb"], box, arm)
            st.write_pdb(str(f.with_suffix(".pdb")))
            subprocess.run([OBABEL, str(f.with_suffix(".pdb")), "-O", str(f.with_suffix(".pdbqt")), "-xr", "-p", "7.4"],
                           check=True, capture_output=True)
            n_atoms = sum(1 for l in f.with_suffix(".pdbqt").read_text().splitlines() if l.startswith(("ATOM", "HETATM")))
            summary[f"{state}__{arm}"] = dict(protein_residues=npro, glycan_residues=ngly, pdbqt_atoms=n_atoms,
                                              kb=round(f.with_suffix(".pdbqt").stat().st_size / 1e3))
            print(f"[receptor] {state} {arm}: {npro} protein residues, {ngly} glycan residues, {n_atoms} atoms", file=sys.stderr)
    (ROOT / "results").mkdir(exist_ok=True)
    if summary:
        (ROOT / "results" / "receptor_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
