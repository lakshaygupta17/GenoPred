# GenoPred V2 data design

## Unit of prediction

The first supported task is gene-level ranking: estimate whether a human gene has credible
pancreatic ductal adenocarcinoma driver evidence. It is not a patient-level mutation predictor.

## Coordinate policy

The canonical build is GRCh38. Tracks supplied as hg19/GRCh37 must be lifted over or regenerated
from raw reads before they are combined with GRCh38 annotations. Unknown builds are blocked.

## Labels

Labels use three states:

- `positive`: PAAD-specific evidence meeting a documented curation rule.
- `unlabelled`: insufficient evidence in either direction.
- `reliable_negative`: adequately measured and passing explicit negative-selection rules.

The full OncoKB pan-cancer gene list must never be treated as PAAD-specific ground truth.
Mutation-derived labels must not be evaluated with mutation-frequency features from the same
cohort; use an independent cohort, orthogonal evidence, or exclude the circular feature.

Stage 4 operationalizes this as a positive-unlabelled table. A positive requires IntOGen PAAD
driver support in at least two qualifying external families among CPTAC, ICGC PACA-AU and ICGC
PACA-CA. TCGA-PAAD evidence is recorded but cannot qualify because TCGA mutation prevalence is a
planned feature. PCAWG is support-only because it overlaps ICGC. No reliable negatives are emitted
until an independently reviewed negative-evidence rule exists.

## Evaluation

One final test partition remains untouched until the model and hyperparameters are frozen.
Repeated evaluation of multiple models on the same test genes is not cross-validation. Results
must distinguish transductive gene ranking from inductive generalization.

Stage 6's fixed three-positive test partition was evaluated in the original run and subsequent
correction/reproducibility reruns. It is no longer untouched. Corrected v2 results are exploratory;
Stage 7 cannot tune against this partition. Repartitioning the same examined genes does not create
independent validation. Final performance claims require genuinely held-out or orthogonal evidence.
