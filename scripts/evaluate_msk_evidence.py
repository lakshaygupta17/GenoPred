"""Provisional external evidence benchmark; observed labels, not true negatives."""
import json,hashlib,sys
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score
from genopred_v2.external_validation import holm

root=Path(sys.argv[1]);out=root/'data/processed/msk_driver_audit_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
audit_path=out/'gene_evidence_audit.csv';frozen=root/'data/processed/external_validation_v1/frozen_predictions.npz'
path=root/'configs/msk_evidence_evaluation_v1.json'
if not path.exists():
 path.write_text(json.dumps({'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'audit_sha256':sha(audit_path),'predictions_sha256':sha(frozen),'models':['reference_nnpu','pu_decision_tree','pu_gradient_boosting','pu_lightgbm','mutation_prevalence_control'],'primary':'Observed-label average precision within uniquely mapped IMPACT341 genes excluding original 17 positives. Evidence positive = author table sequence alteration in at least 3 distinct patients. Others remain unlabelled.','secondary':'Fractional tie-adjusted recall at top 20,50,100; average rank enrichment; no tuning','permutations':2000,'seed':20261104,'multiplicity':'Holm for four one-sided mean-rank permutation tests; mutation control descriptive','interpretation':'Provisional external-evidence consistency. Not a certified independent gold standard, patient diagnostic accuracy, or unseen-gene generalization test. These outcomes are consumed after this run.'},indent=2))
protocol=json.loads(path.read_text());assert sha(audit_path)==protocol['audit_sha256'] and sha(frozen)==protocol['predictions_sha256']
a=pd.read_csv(audit_path);a=a[a.in_IMPACT341 & a.unique_gene_mapping & ~a.existing_positive].copy()
y=a.primary_new_evidence.to_numpy(int);f=np.load(frozen,allow_pickle=False);assert list(f['model_names'])==protocol['models']
idx=pd.Index(f['gene_ids']).get_indexer(a.gene_id);assert (idx>=0).all() and (f['labels'][idx]==-1).all()
x=f['scores'][idx];assert len(a)==326 and y.sum()==80
ranks=np.column_stack([rankdata(x[:,j]) for j in range(x.shape[1])]);stat=ranks[y==1].mean(axis=0)
rng=np.random.default_rng(protocol['seed']);count=np.zeros(5)
for _ in range(protocol['permutations']): count+=(ranks[rng.choice(len(y),int(y.sum()),replace=False)].mean(axis=0)>=stat)
p=(count+1)/(protocol['permutations']+1);adj=holm(p[:4]);rows=[]
for j,name in enumerate(protocol['models']):
 row={'model':name,'observed_label_AP':float(average_precision_score(y,x[:,j])),'mean_positive_percentile':float((stat[j]-1)/(len(y)-1)),'rank_permutation_p':float(p[j]),'holm_p':float(adj[j]) if j<4 else None}
 for k in [20,50,100]:
  threshold=np.sort(x[:,j])[-k];above=x[:,j]>threshold;equal=x[:,j]==threshold;weight=(k-above.sum())/equal.sum()
  row[f'recall_at_{k}']=float((y[above].sum()+weight*y[equal].sum())/y.sum())
 rows.append(row)
pd.DataFrame(rows).to_csv(out/'provisional_model_comparison.csv',index=False)
result={'genes':len(y),'evidence_positives':int(y.sum()),'unlabelled_comparators':int(len(y)-y.sum()),'positive_fraction_reference':float(y.mean()),'results':rows,'protocol_sha256':sha(path),'code_sha256':sha(Path(__file__))}
(out/'provisional_evaluation.json').write_text(json.dumps(result,indent=2))
print(pd.DataFrame(rows).to_string(index=False))
