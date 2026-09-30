# Peptide arm: not run

No peptide candidates were ranked.

The brief asks this arm to import stage 4's peptide screen. Stage 4
(vina_litpcba_virtual_screening_pipeline) screens small molecules on 15
LIT-PCBA targets and has no peptide arm. Nothing was carried forward, so the
positive control for this modality, enfuvirtide (T-20), was not scored and no
percentile exists for it.

Docking a peptide with a small-molecule scoring function that was never fitted
on peptides would give a number and no meaning. That is the reason the arm was
left empty and not filled with a table.

What would have to be true for the arm to be worth running: a peptide docking
or design pipeline validated on a held-out set of peptide-protein complexes,
with a null arm of random sequences of matched length, run against a gp41
heptad-repeat model whose conformational state has been checked against the
literature. Enfuvirtide's approval details and structural basis were not
verified for this repository and are left out on purpose.
