"""Resolve author-selected gene set; quarantine broad LiFD interpretation."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,sys
import pandas as pd
import numpy as np
from sklearn.metrics import average_precision_score
from scipy.stats import rankdata
from genopred_v2.external_validation import holm
root=Path(sys.argv[1]);base=root/'data/processed/autopsy_evidence_v1';out=root/'data/processed/autopsy_evidence_reconciled_v1';out.mkdir(parents=True,exist_ok=True)
genes='KRAS TP53 CDKN2A SMAD4 TGFBR1 TGFBR2 BMPR1B ACVR1B ARID1A ARID1B ARID2 SMARCA4 SMARCB1 KDM6A KMT2A KMT2C KMT2D KMT2E DNMT3A DNMT3B RBM10 SF3B1 CDH1 FBXW7 GNAS KLF4 PCDH15 RNF43 STK11 MYC GATA6 MUC16'.split()
planpath=root/'configs/autopsy_reconciliation_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if not planpath.exists():
 planpath.write_text(json.dumps({'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'reason':'Initial broad LiFD>0 rule yielded 1098 genes, not paper-selected 32; do not use broad extraction as driver validation. Correct author-set interpretation based on Figure 1C, without choosing by performance. Broad results had already been calculated; this is a transparent post-schema-review amendment, not a pristine preregistered analysis.','author_figure1c_genes':genes,'positive_rule':'Gene in explicit Figure 1C list AND at least one LiFD>0 sequence record. Four author-set genes lack such records, so do not assume sequence-driver positives.','background':'Same eligible Table S3 mutation-ascertained gene universe excluding original 17; do not pool with MSK','reconciliation_gate':'Exactly 792 positive rows and 28 sequence-evidence genes from 32 figure genes; source 98 unique gene/protein strings differs from paper 99 variant count, retained as unresolved representation discrepancy.','models':['reference_nnpu','pu_decision_tree','pu_gradient_boosting','pu_lightgbm','mutation_prevalence_control'],'permutations':2000,'seed':20261104,'input_gene_sha256':sha(base/'gene_evidence.csv'),'input_variant_sha256':sha(base/'positive_variant_evidence.csv'),'frozen_predictions_sha256':sha(root/'data/processed/external_validation_v1/frozen_predictions.npz')},indent=2))
plan=json.loads(planpath.read_text());assert sha(base/'gene_evidence.csv')==plan['input_gene_sha256'];assert sha(base/'positive_variant_evidence.csv')==plan['input_variant_sha256']
v=pd.read_csv(base/'positive_variant_evidence.csv');v=v[v.symbol.isin(genes)];assert len(v)==792 and v.symbol.nunique()==28
g=pd.read_csv(base/'gene_evidence.csv');g['broad_lifd_positive']=g.evidence_positive;g['evidence_positive']=g.symbol.isin(v.symbol)
g.to_csv(out/'gene_evidence.csv',index=False);v.to_csv(out/'positive_variant_evidence.csv',index=False)
ev=g[g.eligible].copy();y=ev.evidence_positive.to_numpy(int);f=np.load(root/'data/processed/external_validation_v1/frozen_predictions.npz');x=f['scores'][ev.idx.to_numpy(int)]
assert sha(root/'data/processed/external_validation_v1/frozen_predictions.npz')==plan['frozen_predictions_sha256']
rr=np.column_stack([rankdata(x[:,j]) for j in range(5)]);stat=rr[y==1].mean(axis=0);count=np.zeros(5);rng=np.random.default_rng(plan['seed'])
for _ in range(plan['permutations']):count+=rr[rng.choice(len(y),int(y.sum()),replace=False)].mean(axis=0)>=stat
p=(count+1)/(plan['permutations']+1);adj=holm(p[:4]);results=[]
for j,m in enumerate(f['model_names']):
 row={'model':str(m),'observed_label_AP':float(average_precision_score(y,x[:,j])),'rank_permutation_p':float(p[j]),'holm_p':float(adj[j]) if j<4 else None}
 for k in [100,500]:
  threshold=np.sort(x[:,j])[-k];above=x[:,j]>threshold;equal=x[:,j]==threshold
  row[f'recall_at_{k}']=float((y[above].sum()+(k-above.sum())/equal.sum()*y[equal].sum())/y.sum())
 results.append(row)
pd.DataFrame(results).to_csv(out/'model_comparison.csv',index=False)
report={'status':'Secondary source-reconciled external computational evidence; schema amendment after broad extraction; not independent gold standard','source_patients':30,'source_samples':210,'author_gene_set':32,'sequence_evidence_genes':28,'positive_sequence_rows':792,'unique_gene_protein_strings':98,'paper_reported_unique_variants':99,'genes_without_positive_sequence_rows':sorted(set(genes)-set(v.symbol)),'eligible_genes':len(ev),'new_evidence_positives':int(y.sum()),'new_evidence_symbols':ev.loc[ev.evidence_positive,'symbol'].tolist(),'overlap_msk80':int((ev.evidence_positive & ev.in_msk_80).sum()),'new_outside_msk80':ev.loc[ev.evidence_positive & ~ev.in_msk_80,'symbol'].tolist(),'results':results,'plan_sha256':sha(planpath)}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
(base/'DO_NOT_USE_AS_DRIVER_VALIDATION.md').write_text('The broad LiFD>0 extraction does not reproduce the author-selected driver set. Its model metrics are retained solely for audit and must not be presented as pancreatic driver validation. Use the separately documented reconciled secondary assessment, with its explicit limitations.\n')
