"""Audit published MSK evidence without changing training labels or fitting models."""
import json, hashlib, sys
from pathlib import Path
from datetime import datetime, timezone
import openpyxl
import pandas as pd
import numpy as np

root=Path(sys.argv[1]); raw=root/'data/raw/validation/msk_pdac_2024'; out=root/'data/processed/msk_driver_audit_v1'
out.mkdir(parents=True,exist_ok=True)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
protocol_path=root/'configs/msk_driver_audit_v1.json'
if not protocol_path.exists():
 protocol={'created_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Evidence eligibility audit; no model-score comparison in this audit','primary_rule':'Author-supplied Table 2 sequence mutations in at least 3 distinct patients, within IMPACT341 consistently assayed panel; uniquely mapped to frozen gene universe; exclude all 17 existing positives','sequence_classes':['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','In_Frame_Del','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation'],'secondary':'All author-supplied alteration types and one-patient sequence evidence, reported separately; no merging mouse and human labels','threshold_rationale':'Three-patient recurrence follows the study gene-level enrichment eligibility threshold; not optimized against model performance','independence':'MSK collection is institutionally distinct from consumed TCGA/CPTAC/ICGC cohorts, but cross-study donor identity is not established. Knowledge-base ascertainment is not biologically independent. Genes were previously present as unlabelled training data; this is not an unseen-gene benchmark.','model_policy':'No models retrained, no scores evaluated, no training labels changed; freeze a separate evaluation protocol after eligibility and provenance review','mouse_status':'Publisher supplements HTTP403; PMC served verification HTML; Europe PMC supplementary API reports unavailable. No validated mouse gene table or orthology mapping obtained.'}
 protocol_path.write_text(json.dumps(protocol,indent=2))
protocol=json.loads(protocol_path.read_text())
w=openpyxl.load_workbook(raw/'msk_driver_supplement.xlsx',read_only=True,data_only=True)
def sheet(n):
 rows=list(w.worksheets[n].values);return pd.DataFrame(rows[1:],columns=rows[0])
cohort=sheet(0);alts=sheet(1)
assert cohort.sample_id.is_unique and cohort.patient_id.is_unique
assert alts.Tumor_Sample_Barcode.isin(cohort.sample_id).all()
alts=alts.merge(cohort[['sample_id','patient_id']],left_on='Tumor_Sample_Barcode',right_on='sample_id',validate='many_to_one')
hgnc=pd.read_csv(raw/'hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
approved=hgnc[hgnc.status.eq('Approved')]
current=set(approved.symbol); previous={}
for row in approved.itertuples():
 for name in row.prev_symbol.split('|'):
  if name: previous.setdefault(name,set()).add(row.symbol)
def normalize(symbol):
 if symbol in current:return symbol
 candidates=previous.get(symbol,set())
 return next(iter(candidates)) if len(candidates)==1 else symbol
alts['original_symbol']=alts.Hugo_Symbol
alts['Hugo_Symbol']=alts.Hugo_Symbol.map(normalize)
alts[['original_symbol','Hugo_Symbol']].drop_duplicates().to_csv(out/'symbol_normalization.csv',index=False)
panel=json.loads((raw/'IMPACT341.json').read_text());panel_genes={g['hugoGeneSymbol'] for g in panel['genes']}
assert len(panel_genes)==341
f=np.load(root/'data/processed/external_validation_v1/frozen_predictions.npz',allow_pickle=False)
identity=pd.DataFrame({'gene_id':f['gene_ids'],'symbol':f['gene_symbols'],'existing_label':f['labels']})
unique=identity[~identity.symbol.duplicated(keep=False)].set_index('symbol')
old=set(identity.loc[identity.existing_label.eq(1),'symbol'])
assert len(old)==17
seq=alts[alts.Variant_Classification.isin(protocol['sequence_classes'])]
allcounts=alts.groupby('Hugo_Symbol').patient_id.nunique()
seqcounts=seq.groupby('Hugo_Symbol').patient_id.nunique()
rows=[]
for symbol in sorted(set(alts.Hugo_Symbol)|panel_genes):
 mapped=symbol in unique.index
 n=int(seqcounts.get(symbol,0));total=int(allcounts.get(symbol,0))
 rows.append({'symbol':symbol,'gene_id':unique.loc[symbol,'gene_id'] if mapped else None,'unique_gene_mapping':mapped,'in_IMPACT341':symbol in panel_genes,'existing_positive':symbol in old,'patients_any_author_alteration':total,'patients_sequence_alteration':n,'sequence_evidence_at_least_3':n>=3,'primary_new_evidence':mapped and symbol in panel_genes and symbol not in old and n>=3,'alteration_classes':';'.join(sorted(set(alts.loc[alts.Hugo_Symbol.eq(symbol),'Variant_Classification'])))})
audit=pd.DataFrame(rows)
audit.to_csv(out/'gene_evidence_audit.csv',index=False)
audit[audit.primary_new_evidence].to_csv(out/'eligible_new_sequence_evidence.csv',index=False)
summary={'patients':len(cohort),'alteration_rows':len(alts),'all_author_table_genes':len(allcounts),'all_table_overlap_existing':len(set(allcounts.index)&old),'all_table_genes_outside_existing':len(set(allcounts.index)-old),'sequence_genes':len(seqcounts),'recurrent_sequence_genes_all_panels':int((seqcounts>=3).sum()),'primary_new_evidence_genes':int(audit.primary_new_evidence.sum()),'primary_new_evidence_symbols':audit.loc[audit.primary_new_evidence,'symbol'].tolist(),'panel_genes':len(panel_genes),'eligible_panel_unlabelled_genes':int((audit.in_IMPACT341 & audit.unique_gene_mapping & ~audit.existing_positive).sum()),'unmapped_author_symbols':audit.loc[~audit.unique_gene_mapping & audit.patients_any_author_alteration.gt(0),'symbol'].tolist(),'panel_distribution':cohort.gene_panel.value_counts().to_dict(),'protocol_sha256':sha(protocol_path),'scores_evaluated':False,'training_labels_changed':False,'independence_status':'Provisional external evidence; donor overlap not independently resolved; author-table annotations not revalidated variant-by-variant'}
(out/'audit_summary.json').write_text(json.dumps(summary,indent=2))
sources={'msk_driver_supplement.xlsx':'https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-024-03362-3/MediaObjects/41591_2024_3362_MOESM2_ESM.xlsx','msk_article.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11835752/fullTextXML','IMPACT341.json':'https://www.cbioportal.org/api/gene-panels/IMPACT341?projection=DETAILED'}
sources['hgnc_complete_set.txt']='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
(raw/'source_manifest.json').write_text(json.dumps({name:{'url':url,'sha256':sha(raw/name),'bytes':(raw/name).stat().st_size} for name,url in sources.items()},indent=2))
print(json.dumps(summary,indent=2))
