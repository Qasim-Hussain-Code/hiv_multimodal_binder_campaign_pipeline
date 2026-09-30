# Aptamer arm: not run

Stage 6 was never run, so no aptamer candidate exists. The control for this
modality, UCLA1 (derived from B40), was not scored.

The repository stanford-rna-fold predicts RNA 3D structure from sequence for a
folding competition. It does not dock or rank aptamers against a protein and it
was not used.

An aptamer arm would be worth running if there were a sequence-to-affinity
model validated on held-out aptamer-protein pairs and a nucleic-acid null
drawn from the same length and composition. Neither exists in this project. An
empty modality reported as empty is better than a padded one.
