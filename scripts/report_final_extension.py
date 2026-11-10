"""Verify and summarize the completed final extension experiments."""
import sys,json,hashlib,itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(sys.argv[1]);out=root/'data/processed/final_extension_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((out/'manifest.json').read_text());plan=json.loads((root/'configs/final_extension_v1.json').read_text())
assert manifest['source_sha256']==sha(root/'scripts/final_extension.py')==plan['source_sha256']
for name,digest in manifest['outputs'].items():assert sha(out/name)==digest,name
assert manifest['model_roundtrips']==248 and manifest['all_exact']
assert sha(root/'data/processed/stage5/x.npy')==plan['input_x_sha256']
assert sha(root/'data/processed/stage5/y.npy')==plan['input_y_sha256']
y=np.load(root/'data/processed/stage5/y.npy');splits=json.loads((out/'cv_splits.json').read_text())
for repeat in range(5):
 fs=[s for s in splits if s['repeat']==repeat];assert len(fs)==3
 np.testing.assert_array_equal(np.sort(np.concatenate([s['test'] for s in fs])),np.arange(len(y)))
 for s in fs:assert not set(s['train']) & set(s['test'])
fold=pd.read_csv(out/'feature_fold_metrics.csv');assert len(fold)==180
summary=pd.read_csv(out/'feature_summary.csv');rep=pd.read_csv(out/'feature_repeat_metrics.csv')
lopo=pd.read_csv(out/'lopo_results.csv');assert len(lopo)==68 and lopo.symbol.nunique()==17
held=np.load(out/'lopo_heldout_unlabelled.npy');assert len(held)==1995 and (y[held]==-1).all()
stability=pd.read_csv(out/'stability_summary.csv');prior=pd.read_csv(out/'prior_summary.csv');assert len(prior)==10
groups=['mutation','mutation_cnv','mutation_expression','full'];models=['nnpu','gradient_boosting','lightgbm']
glabels={'mutation':'Mutation only','mutation_cnv':'Mutation + CNV','mutation_expression':'Mutation + expression','full':'Full features'}
mlabels={'nnpu':'nnPU','gradient_boosting':'Gradient boosting','lightgbm':'LightGBM'}
paired=[]
for model in models:
 p=rep[rep.model.eq(model)].pivot(index='repeat',columns='group',values='ap')
 for group in groups[1:]:
  diff=p[group]-p.mutation
  paired.append({'model':model,'comparison':group+' minus mutation','mean_difference':float(diff.mean()),'repeat_min':float(diff.min()),'repeat_max':float(diff.max()),'positive_repeats':int((diff>0).sum())})
pd.DataFrame(paired).to_csv(out/'paired_feature_differences.csv',index=False)
lo_summary=lopo.groupby(['model','group']).agg(mean_single_positive_ap=('ap','mean'),median_average_rank=('average_rank','median'),mean_top100_credit=('recall100','mean'),mean_top500_credit=('recall500','mean')).reset_index()
lo_summary.to_csv(out/'lopo_summary.csv',index=False)
# All candidate weights are based on a fixed held-out subset, not the whole genome.
candidates=pd.read_csv(out/'candidate_stability.csv')
candidates.sort_values(['model','mean_top100_weight','mean_rank'],ascending=[True,False,True]).groupby('model',sort=False).head(30).to_csv(out/'stable_candidate_shortlists.csv',index=False)
plt.rcParams.update({'font.size':10})
fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
for ax,model in zip(axes,models):
 s=summary[summary.model.eq(model)].set_index('group').loc[groups];v=s.ap_mean
 ax.errorbar(v,np.arange(4),xerr=np.vstack([v-s.ap_min,s.ap_max-v]),fmt='o',capsize=4,color='#146c80');ax.set_yticks(range(4),[glabels[g] for g in groups]);ax.invert_yaxis();ax.set_xlim(0,.8);ax.set_title(mlabels[model]);ax.set_xlabel('Observed-label average precision');ax.grid(axis='x',alpha=.2)
fig.suptitle('Independently tuned feature sets: mean and range across five repeat means')
fig.savefig(root/'docs/figures/final_extension_feature_comparison.png',dpi=180);plt.close(fig)
rank=lopo.pivot(index='symbol',columns=['model','group'],values='average_rank');rank=rank.sort_values(('nnpu','full'))
fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
for model,group,marker in [('nnpu','full','o'),('nnpu','mutation','s'),('gradient_boosting','full','^'),('lightgbm','full','D')]:
 ax.scatter(rank[(model,group)],np.arange(len(rank)),label=mlabels[model]+' / '+glabels[group],marker=marker,s=30,alpha=.8)
ax.set_yticks(range(len(rank)),rank.index);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlabel('Held-out positive average rank among 1,996 test genes (lower is better)');ax.set_title('Leave-one-positive-out: difficulty varies by gene');ax.grid(axis='x',alpha=.2);ax.legend(fontsize=8)
fig.savefig(root/'docs/figures/final_extension_lopo.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
for group in ['mutation','full']:
 s=prior[prior.group.eq(group)].sort_values('prior');ax.plot(s.prior,s.ap,marker='o',label=glabels[group])
ax.set_xscale('log');ax.set_xlabel('Assumed positive-class prior');ax.set_ylabel('Mean observed-label AP');ax.set_title('nnPU prior sensitivity: first repeat, three outer folds');ax.legend();ax.grid(alpha=.2)
fig.savefig(root/'docs/figures/final_extension_prior.png',dpi=180);plt.close(fig)
table='\n'.join(f"| {mlabels[r.model]} | {glabels[r.group]} | {r.ap_mean:.4f} | {r.ap_min:.4f}–{r.ap_max:.4f} | {r.recall100:.1%} |" for r in summary.itertuples())
diffrows='\n'.join(f"| {mlabels[r['model']]} | {r['comparison']} | {r['mean_difference']:+.4f} | {r['repeat_min']:+.4f} to {r['repeat_max']:+.4f} | {r['positive_repeats']}/5 |" for r in paired)
lorows='\n'.join(f"| {mlabels[r.model]} / {glabels[r.group]} | {r.median_average_rank:.1f} | {r.mean_single_positive_ap:.4f} | {r.mean_top100_credit:.1%} | {r.mean_top500_credit:.1%} |" for r in lo_summary.itertuples())
generows='\n'.join('| '+name+' | '+' | '.join(f'{row[(m,g)]:.1f}' for m,g in [('nnpu','mutation'),('nnpu','full'),('gradient_boosting','full'),('lightgbm','full')])+' |' for name,row in rank.iterrows())
stabrows='\n'.join(f"| {r.model} | {r.weighted_jaccard_mean:.3f} | {r.weighted_jaccard_min:.3f}–{r.weighted_jaccard_max:.3f} | {r.genes_full_top100_all17} |" for r in stability.itertuples())
priorrows='\n'.join(f"| {glabels[r.group]} | {r.prior:g} | {r.ap:.4f} | {r.recall100:.1%} | {r.recall500:.1%} |" for r in prior.itertuples())
best=summary.loc[summary.ap_mean.idxmax()];nn=summary[summary.model.eq('nnpu')].set_index('group');delta=nn.loc['full','ap_mean']-nn.loc['mutation','ap_mean']
nnpaired=rep[rep.model.eq('nnpu')].pivot(index='repeat',columns='group',values='ap')
expression_delta=nnpaired.full-nnpaired.mutation_expression
doc=f'''# Final extension experiments — 4 November 2026

This is an addendum to FINAL_RESULTS.md. It answers the remaining feature-value, gene-level,
stability and class-prior questions. All results are exploratory analyses on the existing
dataset. No external labels were added, no new external-test claims are made, and the earlier
analysis snapshots are preserved.

## Main numerical findings

- Highest mean AP in this feature comparison: **{mlabels[best.model]} with {glabels[best.group].lower()}, {best.ap_mean:.4f}**.
- Independently tuned nnPU mutation-only AP: **{nn.loc['mutation','ap_mean']:.4f}**.
- Independently tuned full-feature nnPU AP: **{nn.loc['full','ap_mean']:.4f}**, a difference of **{delta:+.4f}**.
- All 17 original positive genes were individually held out and evaluated.
- Candidate stability was measured over 17 positive-deletion fits on the same 1,995 held-out unlabelled genes.
- Five nnPU priors were examined with L2 tuned separately inside training folds.

AP means observed-label average precision, not diagnostic accuracy or precision against proven
negatives. Differences and ranges below are descriptive, not independent confidence intervals.

**Interpretation:** full-feature nnPU improves mean AP over mutation-only learning by {delta:.4f},
but improves over mutation plus expression by only {expression_delta.mean():.4f}. The latter
difference is positive in {int((expression_delta>0).sum())}/5 repeat means and negative in
{int((expression_delta<0).sum())}/5. The mutation-plus-expression nnPU has higher mean top-100
and top-500 recall ({nn.loc['mutation_expression','recall100']:.1%}, {nn.loc['mutation_expression','recall500']:.1%})
than full nnPU ({nn.loc['full','recall100']:.1%}, {nn.loc['full','recall500']:.1%}). Full nnPU remains
the primary AP comparator, while mutation plus expression is a strong simpler alternative
when candidate recovery is the focus. These findings do not establish a universal winner or
that each extra regulatory assay adds value.

## 1. Separately tuned feature combinations

Each of nnPU, PU gradient boosting and PU LightGBM was evaluated with four feature groups.
Every group was tuned independently, addressing the limitation of the previous fixed-setting
ablations. Mutation uses its prevalence and missingness indicator; the CNV combination adds
amplification and deletion prevalence and their flags; the expression combination adds tumour
expression and its flag; full uses all ten biological features plus ten flags.

Five repeats of three outer folds use the same outer gene partitions as the prior robustness
study. Three inner folds choose settings within each outer training set. Imputation/scaling
are fitted only on that training partition. nnPU searches three priors and two L2 penalties,
with 600 epochs; tree models search depth 2 and 3, with six PU bags and ten sampled unlabelled
genes per positive. Exact stochastic inner-fitting seeds differ from the older run for some
tree fits; compare this addendum's feature sets with each other, not tiny differences between
snapshots. No feature group was chosen using outer-test outcomes during fitting.

| Model | Features | Mean AP | Range of repeat means | Mean top-100 recall |
|---|---|---:|---:|---:|
{table}

![Separately tuned feature comparison](figures/final_extension_feature_comparison.png)

Each mean averages three outer folds within a repeat and then five repeats. Top 100 is within
an outer test fold containing about 6,650 genes. The same genes recur across repeats, and
previous exploratory work already informed model-family choice.

### Paired differences relative to mutation-only learning

| Model | Comparison | Mean AP difference | Repeat range | Positive differences |
|---|---|---:|---:|---:|
{diffrows}

These comparisons measure the contribution of a feature group under this dataset, grid and
label definition. They do not prove general biological superiority. Mutation-based labels
can favour mutation-derived features. A benefit from CNV or expression does not establish
benefit from every regulatory assay; the full feature group should be compared to the simpler
combinations before claiming that all assays help.

## 2. Leave-one-positive-out evaluation

Each of the 17 known positives is removed completely from training, not recoded as unlabelled.
The other 16 positives and 17,950 unlabelled genes form the training set. The same 1,995
unlabelled genes are held out for every test. Each test therefore contains one known positive
and 1,995 unlabelled comparators. Three-fold inner selection is repeated using training data
only. nnPU full/mutation-only and full-feature gradient boosting/LightGBM are assessed.

| Model/features | Median positive rank | Mean single-positive AP | Top-100 credit | Top-500 credit |
|---|---:|---:|---:|---:|
{lorows}

With one positive and no ties, AP is reciprocal rank. With ties, AP credits the threshold
containing the entire tie group, while the displayed rank is the average tied rank. The saved
table also contains best/worst tied ranks. Top-k credit is fractional at a tie boundary.
Ranks refer to this fixed 1,996-gene test, not the entire genome.

| Held-out gene | nnPU mutation rank | nnPU full rank | Gradient boosting full rank | LightGBM full rank |
|---|---:|---:|---:|---:|
{generows}

![Leave-one-positive-out results](figures/final_extension_lopo.png)

These results reveal individual genes that are easy or difficult under the chosen features.
They do not independently certify novel-driver discovery. The training sets overlap heavily,
and the background is a fixed sample of unlabelled genes.

## 3. Candidate stability

For each fitted leave-one-positive-out model, we ranked the fixed held-out unlabelled set.
This isolates sensitivity to the available positive examples and the tuning they induce.
It is **not whole-genome candidate stability** and does not independently test random-seed
variability or biological cohort shifts. Training randomness is held fixed between deletions.

At the top-100 boundary, tied candidates share fractional weights summing to exactly 100.
Weighted Jaccard is the sum of pairwise minimum weights divided by the sum of pairwise maximum
weights. It equals 1 for identical selections and 0 for disjoint selections. A high overlap
indicates stability, not correctness.

| Model/features | Mean pairwise weighted Jaccard | Pairwise range | Genes fully selected in all 17 runs |
|---|---:|---:|---:|
{stabrows}

The complete candidate table records mean top-100 weight, number of full inclusions, mean rank
and rank range. A compact 30-gene shortlist per procedure is saved for review. These are
stable candidates within the held-out background, not validated biomarkers.

## 4. nnPU class-prior sensitivity

The prior is an assumption about hidden positive prevalence in PU risk estimation, not an
estimate established by these data. Five values (0.0005, 0.001, 0.003, 0.01, 0.03) span 0.05%
to 3%. For each fixed prior, L2 is tuned separately inside each outer training partition.
This diagnostic uses only the first repeat's three outer folds and does not select a new
production prior from test performance.

| Features | Fixed prior | Mean AP | Mean top-100 recall | Mean top-500 recall |
|---|---:|---:|---:|---:|
{priorrows}

![Prior sensitivity](figures/final_extension_prior.png)

Do not interpret model scores as calibrated probabilities. Stability across priors would
support rank robustness over the tested range, not prove the true biological class proportion.

## Verification and preserved scope

- 64 automated tests passed, including feature-group isolation, held-out-data isolation and tie-budget checks.
- 180 independently tuned outer model/feature evaluations and 68 leave-one-positive-out evaluations completed.
- All 248 saved model bundles reproduced predictions exactly on their evaluation inputs.
- All source/output hashes checked; each gene occurs once per outer test repeat.
- Thirty fixed-prior outer evaluations completed, corresponding to two feature sets x five priors x three folds.
- The original 17-positive labels, original models and all earlier results were preserved.
- The external MSK/autopsy and CRISPR limitations from FINAL_RESULTS.md remain unchanged.

The protocol was saved before fitting at `configs/final_extension_v1.json`. Outputs are under
`data/processed/final_extension_v1/`. Key files: `feature_summary.csv`,
`paired_feature_differences.csv`, `lopo_results.csv`, `stability_summary.csv`,
`candidate_stability.csv`, `stable_candidate_shortlists.csv`, and `prior_summary.csv`.

Reproduce from the project root:

```powershell
.\\.venv\\Scripts\\python.exe scripts/final_extension.py .
.\\.venv\\Scripts\\python.exe scripts/report_final_extension.py .
.\\.venv\\Scripts\\python.exe -m pytest -q
```

The experiment script reuses saved checkpoints only under the same source-code hash. It does
not silently replace prior analyses. The report builder verifies the completed manifest first.
Use this addendum when writing the feature-value and robustness sections of the report.
'''
(root/'docs/FINAL_EXTENSION_RESULTS.md').write_text(doc,encoding='utf-8')
readme=root/'README.md';s=readme.read_text(encoding='utf-8');marker='### 2026-11-04 — Separately tuned features and final sensitivity experiments'
if marker not in s:
 s+='\n'+marker+f'\n\n- Completed independently tuned mutation-only, mutation+CNV, mutation+expression and full-feature comparisons for nnPU, gradient boosting and LightGBM.\n- nnPU mutation-only AP {nn.loc["mutation","ap_mean"]:.4f}; full nnPU AP {nn.loc["full","ap_mean"]:.4f}. Highest mean AP in this comparison: {mlabels[best.model]} / {glabels[best.group]} ({best.ap_mean:.4f}).\n- Completed 17 leave-one-positive-out tests per four procedures, held-out candidate-stability analysis and fixed-prior sensitivity.\n- All 64 tests pass; 248 saved model round trips are exact. Original labels/models remain unchanged.\n- See `docs/FINAL_EXTENSION_RESULTS.md`; this is exploratory sensitivity analysis, not a new independent external test.\n'
 readme.write_text(s,encoding='utf-8')
newfiles=[root/'docs/FINAL_EXTENSION_RESULTS.md',root/'scripts/report_final_extension.py']+list((root/'docs/figures').glob('final_extension_*.png'))+[out/n for n in ['paired_feature_differences.csv','lopo_summary.csv','stable_candidate_shortlists.csv']]
derived={'paired_feature_differences.csv','lopo_summary.csv','stable_candidate_shortlists.csv','report_manifest.json'}
(out/'report_manifest.json').write_text(json.dumps({'verified_core_outputs':sum(n not in derived for n in manifest['outputs']),'tests_passed':64,'model_roundtrips':248,'outputs':{str(p.relative_to(root)):sha(p) for p in newfiles}},indent=2))
print(summary.to_string(index=False));print(lo_summary.to_string(index=False));print(stability.to_string(index=False));print(prior.to_string(index=False));print('Report saved and outputs verified.')
