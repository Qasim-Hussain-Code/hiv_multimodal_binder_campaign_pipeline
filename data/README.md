# data/

Nothing else under data/ is tracked. Everything here is rebuilt by scripts 02 to 05. Full hashes and counts are in logs/02_fetch.json.

## Inputs from earlier stages of the roadmap

Only conventions were imported, not tables. The three stage repositories that exist are listed with the commit that was current on 2026-09-29, the date they were read.

| Stage | Repository | Commit | What this repository took from it |
|---|---|---|---|
| 1, docking | vina_gnina_pose_benchmark_pipeline | 313851a | Vina settings and the ligand and receptor preparation conventions. Code was not vendored. |
| 2, structure prediction | colabfold_boltz_structure_confidence_pipeline | ed49138 | Nothing. It would re-predict designed binders and there are none. |
| 4, virtual screening | vina_litpcba_virtual_screening_pipeline | dc95aaa | Nothing. It is calibrated on 15 LIT-PCBA targets, none of them Env, so it holds no Env candidates. |

Stages 3 (stability simulation), 5 (BindCraft and RFdiffusion) and 6 (aptamers) were never run and have no repository. Stage 4 has no peptide arm.

## Downloaded on 2026-09-29

| File | Source | Content |
|---|---|---|
| structures/4TVP, 5U7O, 5U7M, 6U0L, 8Z7N .cif.gz | files.rcsb.org | mmCIF coordinates. Resolution, deposition date and method are read from data.rcsb.org and written to logs/02_fetch.json. |
| chembl_env_raw.tsv | ChEMBL_37 web API | 680 activity rows, 564 molecules, for eight targets annotated to HIV-1 Env or gp160. Equality-relation IC50, EC50, Ki and Kd in nM. |
| uniprot_env.fasta | rest.uniprot.org | 3,361 non-fragment env sequences of 600 to 900 residues, organism 11676. |

Not downloaded: any Los Alamos alignment, and any CATNAP export. Both need a web form.

## Built by the pipeline

panel/boxes.json (box centres), receptors/ (cropped pdb and pdbqt per state and glycan arm), candidates/ (candidates.tsv with anonymous ids, labels.tsv with the roles, read only by stages 11 to 14), ligands/ (pdbqt), poses/ (top Vina pose per docking).
