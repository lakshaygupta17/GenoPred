"""Frozen secondary external assessment using published autopsy LiFD evidence."""
import json,hashlib,sys,zipfile,io
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import openpyxl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score
from genopred_v2.external_validation import holm

root=Path(sys.argv[1]);raw=root/'data/raw/validation/autopsy_tgfb_2026';out=root/'data/processed/autopsy_evidence_v1';out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
archive=raw/'tgfb_supplements.zip';frozen=root/'data/processed/external_validation_v1/frozen_predictions.npz'
planpath=root/'configs/autopsy_evidence_v1.json'
if not planpath.exists():
 planpath.write_text(json.dumps({'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(archive),'predictions_sha256':sha(frozen),'positive_rule':'At least one numeric LiFD >0 in published Table S3. Do not add thresholds or select by model outcome.','background':'Uniquely mapped genes appearing in Table S3, excluding all original 17 positives; background remains unlabelled, not proven negative. This is mutation-ascertained background, not complete exome coverage.','identity':'Use pinned approved/unique previous HGNC symbols; exact gene-universe mapping only.','models':['reference_nnpu','pu_decision_tree','pu_gradient_boosting','pu_lightgbm','mutation_prevalence_control'],'metrics':'Observed-label AP, fractional tie-adjusted top100/500 recall and one-sided positive mean-rank enrichment','permutations':2000,'seed':20261104,'multiplicity':'Holm across four primary model rank-enrichment tests; control descriptive','secondary':'Report overlap with 80 MSK evidence genes, do not pool studies.','limits':'Johns Hopkins autopsy cohort; cross-study identity not verified. LiFD computational evidence combines existing databases/predictors, not experimental proof. Source never used in fitting; original unlabelled genes were present in training.','source_url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13223546/supplementaryFiles'},indent=2))
plan=json.loads(planpath.read_text());assert sha(archive)==plan['source_sha256'] and sha(frozen)==plan['predictions_sha256']
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 name=next(n for n in z.namelist() if 'table_s3_' in n)
 content=z.read(name);(raw/name).write_bytes(content)
 w=openpyxl.load_workbook(io.BytesIO(content),read_only=True,data_only=True)
 rows=w.active.values;next(rows);next(rows);header=next(rows)
 frame=pd.DataFrame(rows,columns=header).dropna(how='all')
print('Loaded supplementary rows',len(frame),flush=True)
frame['numeric_lifd']=pd.to_numeric(frame.LiFD,errors='coerce')
h=pd.read_csv(root/'data/raw/validation/msk_pdac_2024/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('');h=h[h.status.eq('Approved')];current=set(h.symbol);previous={}
for row in h.itertuples():
 for s in row.prev_symbol.split('|'):
  if s:previous.setdefault(s,set()).add(row.symbol)
def normalize(s):return s if s in current or len(previous.get(s,set()))!=1 else next(iter(previous[s]))
frame['symbol']=frame.Gene_Name.fillna('').astype(str).map(normalize)
f=np.load(frozen);ids=pd.DataFrame({'gene_id':f['gene_ids'],'symbol':f['gene_symbols'],'label':f['labels'],'idx':np.arange(len(f['labels']))})
ids=ids[~ids.symbol.duplicated(keep=False)]
gene=frame.groupby('symbol').agg(source_rows=('Sample','size'),patients=('PAM','nunique'),max_lifd=('numeric_lifd','max')).reset_index()
gene['evidence_positive']=gene.max_lifd.gt(0)
gene=gene.merge(ids,on='symbol',how='left',validate='one_to_one')
gene['eligible']=gene.idx.notna() & gene.label.eq(-1)
msk=pd.read_csv(root/'data/processed/msk_driver_audit_v1/eligible_new_sequence_evidence.csv');gene['in_msk_80']=gene.symbol.isin(msk.symbol)
gene.to_csv(out/'gene_evidence.csv',index=False)
frame[frame.numeric_lifd.gt(0)][['Sample','PAM','Gene_Name','symbol','HGVSp_Short','LiFD','Driver','Clonality']].to_csv(out/'positive_variant_evidence.csv',index=False)
ev=gene[gene.eligible].copy();idx=ev.idx.to_numpy(int);y=ev.evidence_positive.to_numpy(int);x=f['scores'][idx]
assert y.sum()>0 and (y==0).sum()>0
rank=np.column_stack([rankdata(x[:,j]) for j in range(5)]);obs=rank[y==1].mean(axis=0);counts=np.zeros(5);rng=np.random.default_rng(plan['seed'])
for _ in range(plan['permutations']):counts+=rank[rng.choice(len(y),int(y.sum()),replace=False)].mean(axis=0)>=obs
p=(1+counts)/(1+plan['permutations']);adjusted=holm(p[:4]);results=[]
for j,m in enumerate(f['model_names']):
 row={'model':str(m),'observed_label_AP':float(average_precision_score(y,x[:,j])),'rank_permutation_p':float(p[j]),'holm_p':float(adjusted[j]) if j<4 else None}
 for k in [100,500]:
  k=min(k,len(y));cut=np.sort(x[:,j])[-k];above=x[:,j]>cut;equal=x[:,j]==cut
  row[f'recall_at_{k}']=float((y[above].sum()+(k-above.sum())/equal.sum()*y[equal].sum())/y.sum())
 results.append(row)
pd.DataFrame(results).to_csv(out/'model_comparison.csv',index=False)
report={'source_rows':len(frame),'source_patients':int(frame.PAM.nunique()),'source_samples':int(frame.Sample.nunique()),'lifd_positive_rows':int(frame.numeric_lifd.gt(0).sum()),'all_lifd_positive_genes':int(gene.evidence_positive.sum()),'paper_reported_positive_genes':32,'eligible_genes':len(ev),'new_evidence_positives':int(y.sum()),'new_evidence_symbols':ev.loc[ev.evidence_positive,'symbol'].tolist(),'overlap_with_msk_80':int((ev.evidence_positive & ev.in_msk_80).sum()),'positives_outside_both_old17_and_msk80':int((ev.evidence_positive & ~ev.in_msk_80).sum()),'unmapped_positive_symbols':gene.loc[gene.evidence_positive & gene.idx.isna(),'symbol'].tolist(),'positive_fraction':float(y.mean()),'results':results,'plan_sha256':sha(planpath),'source_sha256':sha(archive),'paper_count_reconciled':int(gene.evidence_positive.sum())==32,'status':'External computational evidence; not certified independent ground truth'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
