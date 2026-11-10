"""Build final exploratory results snapshot from verified saved artifacts."""
import json,hashlib,sys,platform,importlib.metadata
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);out=root/'data/processed/final_results_v1';out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
cvdir=root/'data/processed/final_robustness_v1';cv=pd.read_csv(cvdir/'summary.csv').set_index('model')
mskdir=root/'data/processed/msk_driver_audit_v1';msk=pd.read_csv(mskdir/'provisional_model_comparison.csv').set_index('model')
autdir=root/'data/processed/autopsy_evidence_reconciled_v1';aut=pd.read_csv(autdir/'model_comparison.csv').set_index('model');ar=json.loads((autdir/'report.json').read_text())
crdir=root/'data/processed/external_validation_v1';crispr=pd.read_csv(crdir/'model_comparison.csv').set_index('model')
auditdir=root/'data/processed/final_evidence_audit_v1';audit=json.loads((auditdir/'audit_summary.json').read_text())
manifest=json.loads((cvdir/'manifest.json').read_text())
for name,digest in manifest['outputs'].items():assert sha(cvdir/name)==digest,name
assert manifest['source_sha256']==sha(root/'scripts/final_robustness.py')
folds=json.loads((cvdir/'splits.json').read_text());original_y=np.load(root/'data/processed/stage5/y.npy')
for repeat in range(5):
 selected=[r for r in folds if r['repeat']==repeat]
 assert len(selected)==3
 np.testing.assert_array_equal(np.sort(np.concatenate([r['test_idx'] for r in selected])),np.arange(len(original_y)))
 for fold in selected:
  assert not set(fold['train_idx']) & set(fold['test_idx'])
  assert len(fold['train_idx'])+len(fold['test_idx'])==len(original_y)
  assert (original_y[fold['test_idx']]==1).sum() in [5,6]
oof=np.load(cvdir/'oof_scores.npz');assert oof['scores'].shape==(5,19962,13) and np.isfinite(oof['scores']).all()
np.testing.assert_array_equal(oof['labels'],original_y)
assert manifest['model_roundtrips']==60 and manifest['all_exact']
spec=[('nnPU','nnpu','reference_nnpu'),('LightGBM','lightgbm','pu_lightgbm'),('Gradient boosting','gradient_boosting','pu_gradient_boosting'),('Shallow tree','decision_tree','pu_decision_tree'),('Mutation baseline','mutation_control','mutation_prevalence_control')]
table=[]
for display,c,e in spec:
 table.append({'model':display,'nested_AP_mean':float(cv.loc[c,'ap_mean']),'nested_AP_repeat_min':float(cv.loc[c,'ap_min']),'nested_AP_repeat_max':float(cv.loc[c,'ap_max']),'nested_recall100':float(cv.loc[c,'recall100_mean']),'nested_recall500':float(cv.loc[c,'recall500_mean']),'MSK_provisional_AP':float(msk.loc[e,'observed_label_AP']),'MSK_recall100':float(msk.loc[e,'recall_at_100']),'autopsy_secondary_AP':float(aut.loc[e,'observed_label_AP']),'autopsy_recall500':float(aut.loc[e,'recall_at_500']),'CRISPR_primary_rho':float(crispr.loc[e,'rho'])})
combined=pd.DataFrame(table);combined.to_csv(out/'final_comparison.csv',index=False)
plt.rcParams.update({'font.size':10,'axes.titlesize':11})
fig,ax=plt.subplots(1,3,figsize=(13,4.5),layout='constrained');yy=np.arange(5);labels=combined.model
vals=combined.nested_AP_mean
ax[0].errorbar(vals,yy,xerr=np.vstack([vals-combined.nested_AP_repeat_min,combined.nested_AP_repeat_max-vals]),fmt='o',capsize=4,color='#146c80')
ax[1].barh(yy,combined.MSK_provisional_AP,color='#146c80');ax[2].barh(yy,combined.autopsy_secondary_AP,color='#146c80')
titles=['Nested internal evaluation\n17 original positives','MSK provisional evidence\n80 positives / 326 genes','Autopsy secondary evidence\n15 positives / 11,925 genes']
for i,a in enumerate(ax):
 a.set_yticks(yy,labels);a.invert_yaxis();a.set_title(titles[i]);a.set_xlabel('Observed-label average precision');a.grid(axis='x',alpha=.2);a.set_axisbelow(True)
ax[0].set_xlim(0,.65);ax[1].set_xlim(0,.55);ax[2].set_xlim(0,.075)
fig.suptitle('GenoPred V2: results depend on the evaluation question',fontsize=14)
fig.savefig(root/'docs/figures/final_model_results.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.7),layout='constrained');groups=['Full features','Without mutation','Expression only','Missingness only','Shuffled labels'];keys=['','__no_mutation','__expression_only','__missing_flags_only','__shuffled_labels'];yy=np.arange(5)
for offset,m,color in [(-.17,'nnpu','#146c80'),(.17,'gradient_boosting','#d27b38')]:
 ax.barh(yy+offset,[cv.loc[m+k,'ap_mean'] for k in keys],height=.3,label='nnPU' if m=='nnpu' else 'Gradient boosting',color=color)
ax.set_yticks(yy,groups);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlabel('Mean observed-label AP (log scale)');ax.set_title('Mutation features account for most of the internal signal');ax.legend();ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
fig.savefig(root/'docs/figures/final_feature_ablations.png',dpi=180);plt.close(fig)
internal_rows='\n'.join(f"| {r.model} | {r.nested_AP_mean:.4f} | {r.nested_AP_repeat_min:.4f}–{r.nested_AP_repeat_max:.4f} | {r.nested_recall100:.1%} | {r.nested_recall500:.1%} |" for r in combined.itertuples())
external_rows='\n'.join(f"| {r.model} | {r.MSK_provisional_AP:.4f} | {r.MSK_recall100:.1%} | {r.autopsy_secondary_AP:.4f} | {r.autopsy_recall500:.1%} | {r.CRISPR_primary_rho:.4f} |" for r in combined.itertuples())
ablations='\n'.join(f"| {label} | {cv.loc['nnpu'+key,'ap_mean']:.4f} | {cv.loc['gradient_boosting'+key,'ap_mean']:.4f} |" for label,key in zip(groups,keys))
repeats=pd.read_csv(cvdir/'repeat_metrics.csv').pivot(index='repeat',columns='model',values='ap');delta=repeats.nnpu-repeats.mutation_control
doc=f'''# GenoPred V2 — final exploratory results for report writing

**Analysis snapshot: 4 November 2026.** Stages 1–6 are complete. This package adds repeated
internal evaluation and external evidence assessments. It does not establish clinical utility,
fully independent gold-standard validation, or completion of Stage 7 transfer learning.

## Main conclusion

**nnPU is the strongest primary comparator in the completed repeated internal evaluation.**
Its mean observed-label average precision is **0.5517**, compared with **0.5025 for mutation
prevalence alone**, 0.4490 for LightGBM, 0.3160 for gradient boosting and 0.1365 for the shallow
tree. This reverses the apparent advantage of the shallow tree on the earlier three-positive
test split. The earlier result remains an exploratory historical result, not the headline.

Most internal predictive signal depends on mutation features. External evidence is mixed:
gradient boosting has the largest AP point estimate on provisional MSK labels, nnPU has the
largest AP on the small source-reconciled autopsy assessment, and CRISPR survival associations
are weak. The project supports **candidate prioritization and methodological comparison**,
not diagnosis, clinical decision-making or validated discovery of novel drivers.

## Report scope, objectives and question

Suggested title: **GenoPred V2: Reproducible Multi-Omic Feature Integration and
Positive–Unlabelled Learning for Pancreatic Cancer Driver-Gene Prioritization.**

Research question: How reliably can pancreatic molecular features rank known and independently
reported candidate driver genes, and how do linear nnPU, tree methods and simple biological
baselines compare under leakage-controlled evaluation?

Objectives:

1. Build a reproducible gene-level pancreatic dataset with aligned identifiers and assay QC.
2. Construct conservative positive/unlabelled labels without treating absent evidence as negative.
3. Compare ranking procedures and quantify sensitivity to data splits and feature groups.
4. Assess consistency with external clinical-genomic and functional evidence, retaining its limits.

The project contains 19,962 protein-coding genes, 10 biological measurements and 10 missingness
indicators. Labels remain the original 17 IntOGen-supported positives and 19,945 unlabelled genes.
No external genes were added to training. The implemented reference nnPU is **linear** non-negative
PU risk minimization, not a deep neural network. Breast-to-pancreas transfer is not implemented
in this completed core study and must not appear as a completed methodology in the report.

## 1. Repeated nested evaluation: primary internal results

Five fixed random seeds each define three stratified outer folds. Each gene appears in the
outer test set once per repeat. An outer fold contains five or six known positives and about
6,650 unlabelled genes. For each model, three inner folds select hyperparameters by mean
observed-label AP. Imputation and scaling are fitted only on the relevant training partition.
Outer test features and labels do not participate in fitting or selection.

nnPU searches the original three class priors (0.001, 0.003, 0.01) and two L2 penalties
(0.0001, 0.01), with 600 epochs. Each tree family searches depth 2 or 3 using six PU bags per
fit and ten sampled unlabelled genes per training positive. This repeated procedure uses one
fixed seed per fit, unlike the earlier three-seed ensemble benchmark. It compares procedures,
not repeated evaluation of the exact historical fitted model files. The external assessments
continue to use those original frozen predictions.

The table averages metrics within each repeat's three outer folds and then across five repeats.
The range is across repeat means, **not a confidence interval**. These repeats share genes,
labels and feature sources and are not independent biological replications. The design was
created after earlier screening, so this is retrospective robustness evidence.

| Model | Mean observed-label AP | Range of repeat means | Mean top-100 recall | Mean top-500 recall |
|---|---:|---:|---:|---:|
{internal_rows}

Top-k recall is measured within each outer test fold, not among all 19,962 genes. Tied score
cutoffs receive fractional credit rather than arbitrary gene-ID ordering. AP uses known positives
versus unlabelled comparators; it is neither accuracy nor true precision against established
negatives. Do not write “55.2% accuracy.”

nnPU exceeds the mutation baseline AP in all five repeat means; the mean difference is
{delta.mean():.4f}, ranging from {delta.min():.4f} to {delta.max():.4f}. This is a descriptive
comparison. No independently powered superiority claim or formal confidence interval is made.

![Model comparison](figures/final_model_results.png)

The three panels use different gene universes, label prevalences and x-axis limits.
Compare models within a panel, not absolute AP values across panels.

## 2. Feature dependence and negative controls

Ablations reuse settings selected on the full-feature inner folds and refit on outer training
genes only. They diagnose dependence under fixed settings; they are not separately optimized
models and cannot establish that an assay is universally uninformative.

| Diagnostic | nnPU mean AP | Gradient boosting mean AP |
|---|---:|---:|
{ablations}

![Feature ablations](figures/final_feature_ablations.png)

Removing mutation prevalence and its missingness flag greatly reduces AP. Expression-only
and missingness-only models are near the observed-positive prevalence (about 0.00085 overall).
Shuffled-label results are low and variable, as expected with few positives. The defensible
interpretation is that mutation features dominate this dataset. **The study has not demonstrated
a robust benefit of combining all assays over mutation evidence alone.** The small nnPU advantage
over the mutation baseline warrants further study rather than a strong multi-omic superiority claim.

## 3. External assessments: keep the endpoints separate

| Model | MSK AP | MSK top-100 recall | Autopsy AP | Autopsy top-500 recall | CRISPR correlation |
|---|---:|---:|---:|---:|---:|
{external_rows}

### MSK: 80 provisional evidence genes among 326 eligible panel genes

The published 2,336-patient supplementary dataset contains 9,024 alteration rows. The audit
retains sequence-alteration evidence in at least three distinct patients, restricts the universe
to consistently assayed IMPACT341 genes, requires unique gene mapping and excludes all 17
original positives. This leaves 80 provisional evidence positives and 246 unlabelled comparators.

The completed 80-gene review includes supporting alterations, patient counts, source annotation,
pancreatic context, overlap concerns and an inclusion decision. Of 924 sequence-alteration rows
for these genes, 917 match public cBioPortal sample/gene/protein records. Seven remain unmatched
under exact protein-string matching, including splice/frameshift representations and one missing
protein description. They are flagged, not silently corrected. The API contains no explicit
per-variant driver annotation fields for this cohort. **Zero genes were promoted to independently
confirmed driver labels** by this audit.

The paper describes OncoKB oncogenic/likely-oncogenic annotation, but the workbook lacks a
per-row confidence field. The study's pancreatic context is verified at cohort level, not
gene-specific causation. Patient identity across TCGA/CPTAC/ICGC and this cohort is unresolved.
The 80-gene review is a reproducible provenance audit, not an expert literature adjudication of
all 80 genes. Keep its benchmark provisional. Known database annotations can share underlying
knowledge with the training evidence even when clinical cohorts differ.

Gradient boosting has the highest MSK AP point estimate; nnPU retrieves 38/80 evidence genes
in its top 100. No paired superiority test establishes a winner. The original genes were present
as unlabelled training examples, so this is transductive external label confirmation rather than
generalization to genes never encountered during training.

### Additional autopsy cohort: 15 new sequence-evidence genes

An official Europe PMC supplementary archive was acquired for the 2026 Cancer Discovery paper
on TGF-beta signalling in advanced pancreatic cancer. Its cohort contains 30 Johns Hopkins
research-autopsy donors and 210 tumour samples represented in Table S3.

An important source discrepancy was detected. Applying LiFD > 0 to the full supplement produces
1,098 genes, not the publication's displayed 32-gene set. That initial broad extraction and its
metrics are retained only for audit and **excluded from headline driver validation**.

Restricting to the explicit Figure 1C gene list reproduces the reported **792 positive sequence
records**. These span 28 genes; CDH1, GATA6, MUC16 and MYC have no LiFD-positive sequence rows
in that set. After removing the 13 overlapping original positives, **15 new evidence genes**
remain among 11,925 uniquely mapped, mutation-ascertained background genes. Nine overlap the
MSK 80; six do not: BMPR1B, DNMT3A, DNMT3B, KMT2E, PCDH15 and SMARCB1. Do not add these
cohorts' positive counts as if they were disjoint.

The new set is: {', '.join(ar['new_evidence_symbols'])}.

The supplement contains 64,543 data rows versus 64,544 reported in the text. The selected rows
contain 98 unique gene/protein strings versus 99 variants reported in the paper; distinct genomic
variants can share protein representation, so these are not necessarily equivalent counts.
These discrepancies remain documented. Source reconciliation occurred after the broad extraction
had been scored, with a saved amendment based on the published gene set, not on model performance.
Therefore this is a **secondary, source-reconciled assessment**, not a pristine preregistered test.

nnPU AP is 0.0581 and it retrieves 8/15 positives in the top 500. The background positive fraction
is about 0.00126, but background genes are not confirmed negatives. LiFD itself integrates public
databases and computational predictors, so this is additional computational evidence, not an
orthogonal experimental driver validation. Cross-study donor overlap is still not certified.

### CRISPR survival: weak orthogonal support

The prior Sanger Project Score assessment evaluated 15,767 eligible unlabelled non-common-essential
genes across 33 pancreatic adenocarcinoma donors, with 856 other-cancer donors for specificity.
Known PANC-1 donor-related models were excluded. Primary correlations range from 0.0229 to
0.0695 for the four models. None of their top-100 lists is enriched against random gene sets.
Pancreatic specificity is weak. The tree's nominal top-500 expands to 1,684 genes through ties.

Cell-survival dependency is not the same as cancer-driver activity. Tumour suppressors need not
behave as knockout dependencies. These findings neither validate all driver predictions nor
prove every candidate wrong. They limit claims of functional support.

## 4. Why the models disagree

Among originally unlabelled genes, nnPU score correlates 0.662 with mutation prevalence and
0.366 with expression. Gradient boosting correlations are 0.352 and 0.526; LightGBM correlations
are 0.256 and 0.624. These are associations, not causal feature importance.

nnPU and gradient boosting share 53 genes in their respective top 100; nnPU and LightGBM share
82. The shallow tree's tied top-100 threshold includes 103 genes in this whole-unlabelled-universe
comparison. A saved candidate-review table combines ranks, mutation/CNV/expression features,
missingness, MSK evidence and CRISPR outcomes. Use it to illustrate agreement and disagreement,
not to label individual genes as newly discovered drivers.

## 5. Scope limits and unfinished external work

- Original labels are computationally supported and incomplete; there are no reliable negatives.
- Only 17 original positives constrain evaluation precision and biological coverage.
- Prior model screening consumed the old test split; the new nested analysis remains retrospective.
- Dataset-level biological feature processing predates gene folds. The new separation concerns
  gene-level fitting/preprocessing, not independent patient-cohort generalization.
- External genomic label sets have unresolved annotation/overlap issues and have now been seen.
- The mouse Sleeping Beauty screen was not evaluated. Modern and original PNAS supplement URLs
  returned HTTP 403; PMC returned verification HTML and Europe PMC did not supply those files.
  No verified mouse gene table or human-orthologue mapping was available. No result is invented.
- Licensed variant reannotation, donor-level provenance reconciliation and expert gene-specific
  causal review remain future work. Neither a fully independent gold-standard test nor clinical
  validation has been completed.
- Stage 7 transfer learning is outside this frozen core-study report and remains unimplemented.

## 6. Verification and reproducibility

All **58 automated tests pass**, including new checks that perturbing held-out features and labels
cannot change nnPU or gradient-boosting fitting, and that tied scores receive fair top-k credit.
All **60 nested-CV model bundles** reproduce their saved predictions exactly. Every gene is
held out exactly once per repeat; all outer partitions are disjoint and cover the full gene set.
All 5 x 19,962 x 13 stored out-of-fold scores are finite. The original labels remain unchanged.

The earlier 60-row benchmark is preserved. Do not rerun or overwrite old experiments to make
their results agree with the new evaluation. Protocols, fold assignments, inner selections,
models, raw source hashes and output manifests are saved for audit.

| Artifact | Purpose |
|---|---|
| `data/processed/final_results_v1/final_comparison.csv` | Main numerical table |
| `data/processed/final_robustness_v1/` | Nested folds, model bundles, selections, OOF scores, ablations |
| `data/processed/final_evidence_audit_v1/msk_80_gene_quality_review.csv` | Per-gene evidence audit |
| `data/processed/final_evidence_audit_v1/msk_80_gene_variants.csv.gz` | Supporting variants and match flags |
| `data/processed/final_evidence_audit_v1/top_candidate_review.csv` | Model disagreement and biological evidence |
| `data/processed/autopsy_evidence_reconciled_v1/` | Secondary reconciled autopsy evidence and results |
| `configs/final_robustness_v1.json` | Nested-evaluation protocol saved before training |
| `configs/autopsy_reconciliation_v1.json` | Transparent source interpretation amendment |

Reproduction commands from the project folder:

```powershell
.\\.venv\\Scripts\\python.exe scripts/final_robustness.py .
.\\.venv\\Scripts\\python.exe scripts/build_final_results.py .
.\\.venv\\Scripts\\python.exe -m pytest -q
```

Raw Excel-reading audits require openpyxl in addition to the project scientific dependencies.
This run used the bundled openpyxl reader with the project numpy/pandas/scipy environment,
without modifying installed dependencies. Pure model evaluation and report generation use the
project venv. All sources are retained locally for research; source licensing still applies.

## 7. What to say in the report

Recommended result statement:

> GenoPred V2 provides a reproducible pancreatic driver-gene prioritization pipeline. In repeated
> nested gene-level evaluation using 17 conservative positive labels, linear nnPU achieved the
> highest mean observed-label average precision among the shortlisted procedures. Mutation
> prevalence explained most of the observed signal. External computational evidence showed
> enrichment with dataset-dependent model rankings, while functional dependency support was
> weak. These results support exploratory candidate prioritization and motivate stronger
> independent validation, rather than establishing clinical utility or novel driver discovery.

Use nnPU as the primary method for presenting this completed study, mutation prevalence as the
essential baseline, and LightGBM/gradient boosting as substantive alternative models. Report
the shallow tree and its split sensitivity openly. This is a presentation choice based on the
completed analyses, not a new claim of an independently validated best model.

## Primary sources

- [MSK study and supplementary-data description](https://pmc.ncbi.nlm.nih.gov/articles/PMC11835752/).
- [MSK public study](https://www.cbioportal.org/study/summary?id=pdac_msk_2024).
- [Sanger CRISPR dataset methods](https://depmap.sanger.ac.uk/documentation/datasets/wg-crispr-knockout/).
- [Autopsy study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13223546/).
- [Mouse screen, currently unassessed](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341075/).
'''
(root/'docs/FINAL_RESULTS.md').write_text(doc,encoding='utf-8')
outline='''# Report writing plan

Use FINAL_RESULTS.md as the frozen numerical source. Do not describe this project as cancer diagnosis or completed breast-to-pancreas transfer learning.

1. **Introduction:** pancreatic cancer driver prioritization, incomplete labels and why unlabelled genes cannot be assumed negative.
2. **Objectives:** reproducible molecular data integration, conservative evidence labels, model comparison and evaluation of limitations.
3. **Data and preprocessing:** source cohorts, gene identity, genome-coordinate QC, patient-level aggregation, 10 biological features plus missingness, provenance and corrected V1 defects.
4. **Methods:** linear nnPU, PU tree bagging, gradient boosting, LightGBM, mutation baseline, original exploratory screen, repeated nested gene-level evaluation and fixed-setting ablations.
5. **Results:** lead with nested evaluation; then ablations, MSK provisional evidence, secondary source-reconciled autopsy assessment and weak CRISPR support. Each table must name its gene universe, positive count and endpoint.
6. **Discussion:** why the old shallow-tree advantage was split-sensitive, mutation dominance, model disagreement and the difference between driver evidence and survival dependency.
7. **Limitations:** 17 positives, no true negatives, retrospective screening, external annotation and donor-overlap uncertainty, computational rather than experimentally verified labels, incomplete mouse assessment.
8. **Conclusion:** reproducible exploratory prioritization, not validated novel drivers or clinical accuracy.
9. **Future work:** licensed variant annotation, independent donor provenance, expert review, untouched functional/driver evidence and separately preregistered transfer-learning studies.
10. **Appendices:** full 60-experiment historical screen, all fold metrics, per-gene/variant audit, source discrepancy log, test and artifact verification.

Suggested figures: pipeline diagram; final_model_results.png; final_feature_ablations.png; existing external_functional_validation.png. Keep original single-split numbers in the exploratory-screen section or appendix.

Required wording: use observed-label AP rather than accuracy; provisional MSK evidence rather than confirmed new drivers; descriptive repeat range rather than confidence interval; external computational evidence rather than certified independent ground truth.
'''
(root/'docs/REPORT_WRITING_PLAN.md').write_text(outline,encoding='utf-8')
readme=root/'README.md';s=readme.read_text(encoding='utf-8');marker='### 2026-11-04 — Final exploratory evaluation snapshot'
if marker not in s:
 s+='\n'+marker+'\n\n- Completed 5 repeats x 3 outer folds with 3-fold inner model selection for nnPU, shallow PU trees, gradient boosting and LightGBM; added mutation baseline, feature ablations and shuffled-label controls.\n- Mean observed-label AP: nnPU 0.5517, mutation baseline 0.5025, LightGBM 0.4490, gradient boosting 0.3160, shallow tree 0.1365. Mutation features dominate; the historical single-split tree result is not the final headline.\n- Audited all 80 provisional MSK evidence genes: 917/924 supporting variant rows match public records; explicit per-variant driver fields and cross-cohort donor verification remain unavailable.\n- Added a separate 30-donor autopsy evidence assessment: author-selected sequence records yield 15 new genes, of which 9 overlap MSK. Documented source-count discrepancies and quarantined the overly broad LiFD extraction.\n- All 58 tests pass, 60 nested model round trips are exact, and all repeated test partitions cover each gene once.\n- Final results, figures and report outline: `docs/FINAL_RESULTS.md` and `docs/REPORT_WRITING_PLAN.md`. Original labels/models are unchanged. Stage 7 is not implemented; strict independent and clinical validation remain unresolved; mouse evaluation remains unavailable.\n'
 readme.write_text(s,encoding='utf-8')
paths=[root/'docs/FINAL_RESULTS.md',root/'docs/REPORT_WRITING_PLAN.md',root/'docs/figures/final_model_results.png',root/'docs/figures/final_feature_ablations.png',out/'final_comparison.csv']
for folder in [auditdir,autdir]:paths.extend(p for p in folder.iterdir() if p.is_file())
paths.extend(root/'scripts'/p for p in ['final_robustness.py','final_evidence_audit.py','autopsy_evidence.py','reconcile_autopsy_evidence.py','build_final_results.py'])
versions={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','lightgbm','matplotlib']}
(out/'manifest.json').write_text(json.dumps({'date':'2026-11-04','python':platform.python_version(),'versions':versions,'tests_passed':58,'nested_model_roundtrips':60,'nested_artifacts_verified':len(manifest['outputs']),'hashes':{str(p.relative_to(root)):sha(p) for p in paths}},indent=2))
print(combined.to_string(index=False));print('Verified nested artifacts, partitions and finite OOF scores. Final report, figures, README and manifest saved.')
