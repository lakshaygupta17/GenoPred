"""Variant provenance and model-disagreement audit; never promotes unknown evidence."""
import json,hashlib,sys
from pathlib import Path
import pandas as pd
import numpy as np
import openpyxl
from scipy.stats import rankdata,spearmanr

root=Path(sys.argv[1]);out=root/'data/processed/final_evidence_audit_v1';out.mkdir(parents=True,exist_ok=True)
raw=root/'data/raw/validation/msk_pdac_2024'
api=json.loads((raw/'msk_mutations_api.json').read_text())
hgnc=pd.read_csv(raw/'hgnc_complete_set.txt',sep='\t',dtype=str).fillna('');hgnc=hgnc[hgnc.status.eq('Approved')]
current=set(hgnc.symbol);prev={}
for row in hgnc.itertuples():
 for name in row.prev_symbol.split('|'):
  if name:prev.setdefault(name,set()).add(row.symbol)
def norm(s):
 return s if s in current or len(prev.get(s,set()))!=1 else next(iter(prev[s]))
def protein(s):return str(s).removeprefix('p.')
keys={(r['sampleId'],norm(r['gene']['hugoGeneSymbol']),protein(r.get('proteinChange',''))) for r in api}
w=openpyxl.load_workbook(raw/'msk_driver_supplement.xlsx',read_only=True,data_only=True)
rows=list(w.worksheets[1].values);v=pd.DataFrame(rows[1:],columns=rows[0]);v['symbol']=v.Hugo_Symbol.map(norm)
cohort_rows=list(w.worksheets[0].values);cohort=pd.DataFrame(cohort_rows[1:],columns=cohort_rows[0]);patients=cohort.set_index('sample_id').patient_id.to_dict()
audit=pd.read_csv(root/'data/processed/msk_driver_audit_v1/eligible_new_sequence_evidence.csv')
protocol=json.loads((root/'configs/msk_driver_audit_v1.json').read_text())
v=v[v.symbol.isin(audit.symbol)&v.Variant_Classification.isin(protocol['sequence_classes'])].copy()
v['patient_id']=v.Tumor_Sample_Barcode.map(patients)
v['variant_present_in_cbioportal']=[(r.Tumor_Sample_Barcode,r.symbol,protein(r.Alteration)) in keys for r in v.itertuples()]
v['driver_annotation_status']='Publication-level claim; no per-variant driver field in supplementary table or public API'
v.to_csv(out/'msk_80_gene_variants.csv.gz',index=False)
rows=[]
for a in audit.itertuples():
 subset=v[v.symbol.eq(a.symbol)]
 assert subset.patient_id.nunique()==a.patients_sequence_alteration
 counts=subset.groupby('Alteration').patient_id.nunique().sort_values(ascending=False)
 rows.append({'symbol':a.symbol,'gene_id':a.gene_id,'patients':int(a.patients_sequence_alteration),'distinct_alterations':len(counts),'top_alterations_with_patient_count':'; '.join(f'{name} ({n})' for name,n in counts.head(5).items()),'table_rows':len(subset),'cbioportal_variant_matches':int(subset.variant_present_in_cbioportal.sum()),'pdac_context':'Observed in published PDAC cohort; gene-specific causation not separately established','annotation_source':'Author Table 2; manuscript OncoKB v4.12 methods','cross_cohort_overlap':'Unresolved donor identity across source cohorts','decision':'Retain provisional external-evidence analysis only; do not add to training positives'})
pd.DataFrame(rows).to_csv(out/'msk_80_gene_quality_review.csv',index=False)
f=np.load(root/'data/processed/external_validation_v1/frozen_predictions.npz');x=np.load(root/'data/processed/stage5/x.npy')
names=json.loads((root/'data/processed/stage5/feature_names.json').read_text())
keep=f['labels']==-1;df=pd.DataFrame({'gene_id':f['gene_ids'][keep],'symbol':f['gene_symbols'][keep]})
scores=f['scores'][keep];models=f['model_names'].astype(str)
rank=[];top=[];corr=[]
for j,m in enumerate(models):
 rr=rankdata(-scores[:,j],method='average');rank.append(rr);df[m+'_rank']=rr
 top.append(scores[:,j]>=np.sort(scores[:,j])[-100])
 for k in [6,9]:
  valid=x[keep,k+10]==0
  corr.append({'model':m,'feature':names[k],'genes':int(valid.sum()),'spearman':float(spearmanr(scores[valid,j],x[keep][valid,k]).statistic)})
 missing=x[keep,10:].sum(axis=1)
 corr.append({'model':m,'feature':'number_of_missing_assays','genes':len(df),'spearman':float(spearmanr(scores[:,j],missing).statistic)})
rank=np.column_stack(rank);top=np.column_stack(top)
df['nnpu_minus_gb_rank']=rank[:,0]-rank[:,2]
df['in_any_shortlist_top100']=top[:,:4].any(axis=1)
for k in [6,7,8,9]:df[names[k]]=x[keep,k];df[names[k]+'_missing']=x[keep,k+10].astype(bool)
df['missing_assay_count']=missing
df=df.merge(audit[['gene_id','patients_sequence_alteration']],on='gene_id',how='left')
functional=pd.read_csv(root/'data/processed/external_validation_v1/gene_evidence.csv.gz',usecols=['gene_id','median_pancreatic_fold_change','selectivity','common_essential_or_core'])
df=df.merge(functional,on='gene_id',how='left',validate='one_to_one')
df.to_csv(out/'gene_disagreement_evidence.csv.gz',index=False)
df[df.in_any_shortlist_top100].sort_values('nnpu_minus_gb_rank').to_csv(out/'top_candidate_review.csv',index=False)
pd.DataFrame(corr).to_csv(out/'score_feature_correlations.csv',index=False)
overlap=[]
for i in range(4):
 for j in range(i+1,4):
  overlap.append({'model_a':models[i],'model_b':models[j],'top100_actual_a':int(top[:,i].sum()),'top100_actual_b':int(top[:,j].sum()),'intersection':int((top[:,i]&top[:,j]).sum()),'jaccard':float((top[:,i]&top[:,j]).sum()/(top[:,i]|top[:,j]).sum())})
pd.DataFrame(overlap).to_csv(out/'top100_overlap.csv',index=False)
summary={'reviewed_genes':len(audit),'reviewed_variant_rows':len(v),'api_matching_variant_rows':int(v.variant_present_in_cbioportal.sum()),'public_api_mutation_rows':len(api),'api_rows_with_driver_annotation':sum(any(r.get(k) for k in ['driverFilter','driverTiersFilter','driverFilterAnnotation','driverTiersFilterAnnotation']) for r in api),'promoted_to_confirmed_independent_positive':0,'source_api_sha256':hashlib.sha256((raw/'msk_mutations_api.json').read_bytes()).hexdigest(),'source_api_url':'https://www.cbioportal.org/api/molecular-profiles/pdac_msk_2024_mutations/mutations/fetch?projection=DETAILED','source_api_post_body':{'sampleListId':'pdac_msk_2024_sequenced'},'gene_specific_review':'All 80 gene records reviewed by rule-based variant/patient/provenance audit; not an expert literature adjudication of biological causality','mouse_download_status':'Original and modern PNAS supplement routes returned HTTP403; no verified mouse gene table'}
(out/'audit_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
