# Pre-registration: small-molecule arm, HIV-1 Env, tier A

Written 2026-10-01 by scripts/13_prereg.py, before any assay result exists.
Tree commit: `88cb501d81b059061ff0246469018e4bd8f16313`. Uncommitted changes at the moment of writing: 1 paths.
Bench tier: A (no wet lab). Nothing below has been tested at a bench.

Scope. Only the small-molecule modality was ranked. The peptide, designed-binder and aptamer arms were not
run (results/modality_*_not_run.md), so this document makes no prediction about them and no cross-modality
statement is made.

## What is predicted

The 10 ranked candidates below are the compounds this ranking would order first. The ranking is a Vina score
(AutoDock Vina 1.2.7, empirical scoring function, ranking key only) into the temsavir pocket of four deposited Env
states, glycosylated arm, median unfiltered rank percentile across states, after removing candidates that are
escape-prone in two or more states. The score is not an affinity and not a potency estimate.

| order | id | kind | median percentile across states | ChEMBL id |
|---|---|---|---|---|
| 1 | C0070 | ranked candidate | 7.2 | CHEMBL4865734 |
| 2 | C0047 | ranked candidate | 12.3 | CHEMBL1645256 |
| 3 | C0024 | ranked candidate | 18.1 | CHEMBL1645265 |
| 4 | C0116 | ranked candidate | 19.6 | CHEMBL1645111 |
| 5 | C0151 | ranked candidate | 19.6 | CHEMBL4579482 |
| 6 | C0148 | ranked candidate | 19.6 | CHEMBL4547005 |
| 7 | C0097 | ranked candidate | 22.5 | CHEMBL1645263 |
| 8 | C0103 | ranked candidate | 23.9 | CHEMBL4563477 |
| 9 | C0102 | ranked candidate | 24.6 | CHEMBL1645291 |
| 10 | C0013 | ranked candidate | 25.4 | CHEMBL1645287 |
|  | C0101 | anchor temsavir_83J | 88.4 | CHEMBL3301620 |
|  | C0117 | anchor BMS378806_83G | 55.8 | CHEMBL337301 |
|  | C0168 | null negative |  | CHEMBL1970028 |
|  | C0076 | null negative |  | CHEMBL3985764 |
|  | C0057 | null negative |  | CHEMBL3552224 |
|  | C0048 | null negative |  | CHEMBL3985901 |
|  | C0033 | null negative |  | CHEMBL3623166 |

Ordering 17 compounds in total: 10 ranked candidates, 2 anchors, 5 null negatives.
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
2. At least 3 of the 10 ranked candidates reach IC50 below 10 uM on BG505 at 50 percent neutralisation with a
   full dose-response curve.
3. No more than 1 of the 5 null negatives reaches the same threshold.
The 3 of 10, the 10 uM line and the 1 uM anchor line are arbitrary and were fixed by eye.

## Result that would refute the ranking

Any of the following, with the anchor valid:
1. 0 of the 10 ranked candidates reach IC50 below 10 uM on BG505.
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
| results/dock_scores.tsv | `f321fd5a5dfa48fe8eca630ffb93e775f37cc266143f6ff49749f8314f14faec` |
| results/candidate_labels.tsv | `74cc8adb1189d6f9770a413b073b49e8045b6663dbf2ffc3a0cc3166c0fd3824` |
| results/candidates_anonymous.tsv | `9ff7eeacc7eca73c354b456995a1b405ea3feadb800f55708ed875d769a6f541` |
| results/conservation_contacts.tsv | `4dc32d87c32c7348b6f32b8074bf2d91833362d4e02688600f6a9be78e6ca73a` |
| results/conservation_entropy.tsv | `8967365af8b68291a699d15a21cf03f51012d8f30db26eed83c17152fb53ee03` |
| results/ranks_long.tsv | `e790520eda35dd71b1429feaf7d80de42822e8ed7514150e1c9a2374e744a833` |
| results/control_ranks.tsv | `5bf1cbe6aa6629a57f4a2dd6a3fefba840a1ebebc87a1f6a459be460370af1a8` |
| results/null_floor.tsv | `09b0d1bb698eea324c828f0dd1332852ce6284cd43a17ecca1fa423f621adb0b` |
| results/funnel.tsv | `5f6798abf9c7aa2c31269917ac49176912e1c09ff09a53783e29cc1036369639` |
| results/state_rank_correlation.tsv | `b5590ad7069566cf0e3d4737c903afe17b6af41c74944d74d6ba53de1d89134a` |
| results/glycan_effect.tsv | `2e40b6d53e5555a89652d84973ad4364668f8d256571205c8f93aaaea661042e` |
| results/potency_association.tsv | `661c8b1cfd405cc0ca3af0ffd867e564c446ce09d6f1c7029be798d5f87858be` |
| results/redock_temsavir_S2.tsv | `f0949fbfdaf3b43465233a3507772a7c416f18cb1e97b3d8fbbd7ba26486a6af` |
| config/panel.tsv | `24b23dbc06192ee0a493cd2d11482b2ae8fff12417a14997d2913315b745703f` |
| config/boxes.json | `a35eb5b3d411c9f2d94b09cceeab9019988ce77f109b00203cb51a0b537d092d` |
| config/candidate_sources.tsv | `6576be80b880d8f65c43ef475c44e862badb8630d8e2216547e7cde9a760cdaa` |
| config/cost_assumptions.tsv | `dd2541c7be0f4ca198678943bc5974db0d6138a34478028091c69654ed78a650` |

## Corrections

None. A correction is a new file that names this one and the commit above.
