import json,hashlib,sys
from pathlib import Path
import pandas as pd
import numpy as np
root=Path(sys.argv[1]);out=root/'data/processed/msk_driver_audit_v1';raw=root/'data/raw/validation/msk_pdac_2024'
a=json.loads((out/'audit_summary.json').read_text());e=json.loads((out/'provisional_evaluation.json').read_text())
g=pd.read_csv(out/'gene_evidence_audit.csv');p=g[g.primary_new_evidence]
assert len(p)==80 and p.gene_id.is_unique and (~p.existing_positive).all() and p.in_IMPACT341.all() and p.patients_sequence_alteration.ge(3).all()
for name,item in json.loads((raw/'source_manifest.json').read_text()).items(): assert hashlib.sha256((raw/name).read_bytes()).hexdigest()==item['sha256']
f=np.load(root/'data/processed/external_validation_v1/frozen_predictions.npz');original=np.load(root/'data/processed/stage6_benchmark_v1/all_scores.npz')
idx=[list(original['model_names']).index(n) for n in list(f['model_names'][:4])+['feature_06']]
np.testing.assert_array_equal(f['scores'],original['scores'][:,idx])
rows='\n'.join(f"| {x['model']} | {x['observed_label_AP']:.4f} | {x['recall_at_50']:.1%} | {x['recall_at_100']:.1%} |" for x in e['results'])
doc=f'''# MSK pancreatic driver-evidence audit and provisional benchmark

Completed 4 November 2026. This extends the Stage 8 evidence assessment and does not change
the completed Stage 1–6 dataset, labels or fitted models. Stage 7 remains pending.

## Answer to the dataset question

The public MSK pancreatic cohort supplies **80 additional genes with recurrent sequence-alteration
evidence** under the rules below. This is substantially more than 17, but these are provisional
external evidence positives, not 80 newly experimentally confirmed pancreatic cancer drivers.
We ran a frozen-model comparison and preserved the original training labels.

## Sources and what was downloaded

- [Nature Medicine paper](https://www.nature.com/articles/s41591-024-03362-3), published January 2025,
  study ID `pdac_msk_2024` (the dataset ID retains 2024).
- [Author supplementary workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-024-03362-3/MediaObjects/41591_2024_3362_MOESM2_ESM.xlsx):
  Table 1 cohort metadata and Table 2 somatic alterations.
- [Official IMPACT341 panel definition](https://www.cbioportal.org/api/gene-panels/IMPACT341?projection=DETAILED).
- [HGNC approved and previous symbols](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt).

Raw files and SHA-256 checksums are in `data/raw/validation/msk_pdac_2024/`.
The manuscript describes oncogenic/likely-oncogenic OncoKB annotations as drivers. The table
does not carry a per-row oncogenicity field, so we inherit the publication's annotation claim;
we have not independently reannotated each variant. Broad copy-number events can involve
neighbouring genes, which is another reason not to use every Table 2 gene as a driver label.

## Counts and exclusions

| Step | Count |
|---|---:|
| Patients, one selected sample per patient | 2,336 |
| Published alteration rows | 9,024 |
| Original distinct gene-symbol strings across all alteration types | 396 |
| Distinct symbols after approved/unique previous-symbol normalization | 390 |
| Normalized all-type symbols matching existing positives | 14 |
| Normalized all-type symbols outside existing positives | 376 |
| Genes with eligible sequence-alteration classes | 209 |
| Sequence-evidence genes in at least three patients, across all panels | 107 |
| Genes on the common IMPACT341 panel | 341 |
| Common-panel genes uniquely matched to GenoPred and excluding existing positives | 326 |
| New recurrent sequence-evidence positives within that universe | **80** |
| Remaining unlabelled comparison genes | 246 |

The initial 382 figure referred to raw symbols outside the original list. It is not the final
positive count. Aliases, alteration types, panel coverage, recurrence and gene mapping matter.
Historical names are resolved only through exact approved HGNC symbols or unambiguous previous
symbols. Ambiguous aliases and current symbols absent from GENCODE v36 are retained in the
audit and excluded rather than guessed. This is conservative and may omit usable evidence.

The primary evidence set excludes copy-number alterations, fusions, flank and splice-region
annotations. Eligible sequence classes include missense, nonsense, frameshift, splice-site,
in-frame insertion/deletion, translation-start and nonstop mutations. A gene must occur in at
least three distinct patients. The three-patient threshold follows the paper's gene-enrichment
eligibility rule and was not selected using model performance. The common panel avoids treating
unassayed genes as unlabelled comparators. Panel presence does not prove equal callability or
detection sensitivity in every tumour.

Example new evidence genes include BRCA1, BRCA2, PALB2, BRAF, PIK3CA, PTEN, FBXW7, ARID2,
SMAD2, SMAD3 and TGFBR1. The complete 80-gene list is saved, not manually entered.

## Frozen-model results

Models were not retrained or tuned. The original four-model shortlist and mutation-prevalence
control were evaluated. Eligibility was settled and `configs/msk_evidence_evaluation_v1.json`
was saved before model-score comparison. All 17 original positives were excluded, not just
the original training subset. Evaluation is restricted to the 326 eligible panel genes.

| Model | Observed-label average precision | Recall at top 50 | Recall at top 100 |
|---|---:|---:|---:|
{rows}

The evidence-positive fraction is 80/326 = 24.54%, a prevalence reference for AP, not an exact
finite-sample random AP estimate. Recall uses fractional credit at tied cutoffs. For example,
nnPU recovers 38 of the 80 evidence genes in its top 100. Gradient boosting has the highest
AP point estimate, while nnPU has higher top-100 recall. No pairwise superiority or universal
winner claim follows from these results. Values are not comparable directly to the old
genome-wide three-positive test, because the universe and label definition differ.

Four one-sided mean-rank enrichment tests used 2,000 gene-label permutations with Holm
adjustment. Adjusted p is approximately 0.002 for nnPU, gradient boosting and LightGBM, and
0.0055 for the decision tree. These are tests of enrichment against exchangeable labels,
not estimates of clinical accuracy. Dependencies between genes, panel ascertainment,
expression and mutation-rate confounding are not controlled by these permutations.

## Can we call this independent validation?

**Only provisionally as external evidence, not a certified independent gold standard.**

1. The study is an MSK patient cohort rather than the TCGA/CPTAC/ICGC cohort families used in
   GenoPred. However, cross-study donor identities have not been reconciled. Distinct study
   names and different public patient-ID formats do not prove zero overlap.
2. OncoKB is a curated cancer knowledge base. Its annotations can draw on biology and studies
   related to existing training evidence. New clinical observations do not automatically
   mean independent knowledge-based labels or pancreatic-specific causation.
3. The candidate genes were present as unlabelled examples during PU training. This is a
   transductive label-confirmation exercise, not performance on entirely unseen genes.
4. The other 246 genes remain unlabelled, not proven negatives. AP is observed-label AP;
   specificity, clinical diagnostic accuracy and true precision cannot be inferred.
5. The supplementary workbook has no explicit per-variant driver-confidence column. Review
   or licensed reannotation is needed before upgrading provisional evidence labels.
6. These outcomes have now been seen. Future choices based on them require a fresh evaluation
   source rather than reusing this benchmark as an untouched final test.

## Mouse benchmark status

[Mann et al., PNAS 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341075/) reports 543 mouse
pancreatic candidate cancer genes from a Sleeping Beauty mutagenesis screen. This is a separate
experimental endpoint and must not be pooled with human clinical labels.

The supplementary files could not be retrieved in this session: PNAS returned HTTP 403,
PMC supplied verification HTML instead of Excel files, and Europe PMC's supplementary API
reported the article unavailable through that service. Those responses were not accepted as
data. No mouse candidate count after human-orthologue mapping and no mouse evaluation is
claimed. To complete that benchmark, obtain the original supplementary gene tables through
an authorized download, map one-to-one human orthologues with a pinned source, remove original
positive genes and predefine the background universe before examining model scores.

## Deliverables and reproduction

- `configs/msk_driver_audit_v1.json`: initial eligibility protocol; symbol normalization uses
  pinned HGNC approved/unique previous symbols as documented here and in the saved code.
- `configs/msk_evidence_evaluation_v1.json`: frozen evaluation rules and input hashes.
- `data/processed/msk_driver_audit_v1/gene_evidence_audit.csv`: all candidate/panel genes and reasons.
- `eligible_new_sequence_evidence.csv`: the 80 eligible genes and donor counts.
- `symbol_normalization.csv`: original-to-normalized gene symbols.
- `audit_summary.json`: evidence counts and unresolved independence limitations.
- `provisional_model_comparison.csv`, `provisional_evaluation.json`: exploratory results.

Run `scripts/audit_msk_driver_evidence.py PROJECT_ROOT` with Python containing openpyxl,
pandas and numpy. Run `scripts/evaluate_msk_evidence.py PROJECT_ROOT` with the project venv.
Run `scripts/summarize_msk_audit.py PROJECT_ROOT` to verify and regenerate this report.
The evaluation rejects changed evidence-table or prediction hashes. The audit script does
not modify original training labels. The saved source manifest pins downloaded bytes.

Verification passed: all cohort sample IDs join, selected samples and patients are unique,
all 80 evidence genes are uniquely mapped and meet recurrence/panel/exclusion rules,
source hashes match, and frozen predictions exactly equal their original benchmark columns.

## Next scientific step

Resolve variant-level driver annotation and cohort overlap before treating these 80 genes as
high-confidence independent test positives. Preserve this provisional result regardless of
whether later stricter filtering improves or reduces performance. Mouse assessment remains
pending source acquisition. The existing 17-positive training dataset remains unchanged.
'''
(root/'docs/MSK_DRIVER_EVIDENCE_AUDIT.md').write_text(doc,encoding='utf-8')
readme=root/'README.md';s=readme.read_text(encoding='utf-8'); marker='### 2026-11-04 — MSK driver-evidence audit and provisional benchmark'
if marker not in s:
 s+='\n'+marker+'\n\n- Downloaded and pinned the MSK 2,336-patient supplementary tables, common IMPACT341 panel and HGNC symbol mapping.\n- Identified 80 new recurrent sequence-evidence genes among 326 eligible panel genes after excluding original positives. These are provisional evidence labels, not 80 independently confirmed PDAC drivers.\n- Evaluated frozen models without retraining: observed-label AP 0.4179 nnPU, 0.3490 decision tree, 0.4365 gradient boosting and 0.4056 LightGBM. nnPU top-100 recall is 47.5%.\n- Saved the full audit, exclusions, source hashes, frozen protocol and results in `docs/MSK_DRIVER_EVIDENCE_AUDIT.md` and `data/processed/msk_driver_audit_v1/`.\n- Cross-study donor overlap and variant-level annotation remain unresolved. Mouse screen preparation is blocked by unavailable supplement downloads. No fully independent gold-standard claim is made.\n- Original labels and models are unchanged; Stage 7 remains pending.\n'
 readme.write_text(s,encoding='utf-8')
paths=list(out.glob('*'))+[root/'docs/MSK_DRIVER_EVIDENCE_AUDIT.md']+[root/'scripts'/name for name in ['audit_msk_driver_evidence.py','evaluate_msk_evidence.py','summarize_msk_audit.py']]
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file() and p.name!='artifact_manifest.json'}
(out/'artifact_manifest.json').write_text(json.dumps(hashes,indent=2))
print('Verified evidence rules, source hashes and unchanged model predictions; saved report and README.')
