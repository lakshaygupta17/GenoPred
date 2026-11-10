"""Separately tuned feature sets, nnPU leave-one-positive-out and prior/stability diagnostics."""
import sys,json,hashlib,itertools,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from genopred_v2.stage6 import impute_from_partition,train_nnpu,_sigmoid
from genopred_v2.stage6_benchmark import load_inputs,estimator,fit_classifier,model_score

GROUPS={'mutation':[6,16],'mutation_cnv':[6,7,8,16,17,18],'mutation_expression':[6,9,16,19],'full':list(range(20))}
MODELS=['nnpu','gradient_boosting','lightgbm']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def options(model,prior=None):
 return [dict(prior=p,l2=l) for p,l in itertools.product([prior] if prior is not None else [.001,.003,.01],[.0001,.01])] if model=='nnpu' else [dict(level=0),dict(level=1)]
def fit(x,y,mask,model,group,opt,seed):
 clean,medians=impute_from_partition(x,mask,10);cols=GROUPS[group];scaler=StandardScaler().fit(clean[mask][:,cols]);xx=scaler.transform(clean[mask][:,cols]);yy=y[mask]
 state=dict(model=model,group=group,medians=medians,scaler=scaler,opt=opt)
 if model=='nnpu':
  w,b,risk=train_nnpu(xx,yy,opt['prior'],opt['l2'],.03,600);state.update(w=w,b=b,risk=risk)
 else:
  rng=np.random.default_rng(seed);p=np.flatnonzero(yy==1);u=np.flatnonzero(yy==-1);models=[]
  for _ in range(6):
   ids=np.r_[p,rng.choice(u,min(len(u),10*len(p)),replace=False)]
   m=estimator(model,opt['level'],int(rng.integers(0,2**31-1)));models.append(fit_classifier(m,xx[ids],(yy[ids]==1).astype(int)))
  state['estimators']=models
 return state
def score(state,x):
 v=x.copy().astype(float)
 for i,median in enumerate(state['medians']):v[v[:,10+i]==1,i]=median
 v=state['scaler'].transform(v[:,GROUPS[state['group']]])
 return _sigmoid(v@state['w']+state['b']) if state['model']=='nnpu' else np.mean([model_score(m,v) for m in state['estimators']],axis=0)
def metrics(y,s):
 result={'ap':float(average_precision_score(y==1,s))}
 for k in [100,500]:
  k=min(k,len(s));t=np.sort(s)[-k];a=s>t;e=s==t
  result[f'recall{k}']=float(((y[a]==1).sum()+(k-a.sum())/e.sum()*(y[e]==1).sum())/(y==1).sum())
 return result
def select(x,y,train,model,group,seed,prior=None):
 trials=[];inner=list(StratifiedKFold(3,shuffle=True,random_state=seed+1000).split(train,y[train]))
 for opt in options(model,prior):
  vals=[]
  for k,(a,b) in enumerate(inner):
   mask=np.zeros(len(y),bool);mask[train[a]]=True
   state=fit(x,y,mask,model,group,opt,seed+k)
   vals.append(float(average_precision_score(y[train[b]]==1,score(state,x[train[b]]))))
  trials.append({'option':opt,'inner_ap':vals,'mean':float(np.mean(vals))})
 return max(trials,key=lambda r:r['mean'])['option'],trials
def top_weight(s,k):
 t=np.sort(s)[-k];a=s>t;e=s==t;return a.astype(float)+e*(k-a.sum())/e.sum()

def main(root):
 root=Path(root);out=root/'data/processed/final_extension_v1';out.mkdir(parents=True,exist_ok=True)
 x,y,names,genes,symbols,*_=load_inputs(root);p=np.flatnonzero(y==1);u=np.flatnonzero(y==-1)
 planpath=root/'configs/final_extension_v1.json'
 if not planpath.exists():
  planpath.write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'groups':GROUPS,'models':MODELS,'seeds':[113,227,349,463,587],'outer_folds':3,'inner_folds':3,'selection':'Each model/feature set independently tuned by inner AP; ties choose first candidate. Fixed 600 nnPU epochs; six PU bags.','nnpu_grid':options('nnpu'),'tree_grid':options('gradient_boosting'),'lopo':'nnPU full and mutation-only independently tuned on 16 remaining positives and fixed training-unlabelled set. Hold out fixed 10 percent of U as common test background; held positive never included as unlabelled. Compare all 17 held-positive ranks.','prior_sensitivity':'nnPU full and mutation on first repeat outer folds; fixed priors .0005,.001,.003,.01,.03; L2 tuned within each outer training set. Secondary diagnostic, not model selection.','stability':'Fit nnPU/GB/LightGBM full and nnPU mutation to same 16-positive deletion training sets; rank original U within fixed heldout U background. Top100 fractional tie weights and pairwise weighted Jaccard. Sensitivity to positive deletion, not independent biological replication.','lopo_seed':20261105,'reuse':'Reuse complete checkpoint only when source hash matches plan; no original artifacts overwritten.','interpretation':'Additional exploratory experiments on consumed data; no new independent test. No label expansion or external outcome tuning.','input_x_sha256':sha(root/'data/processed/stage5/x.npy'),'input_y_sha256':sha(root/'data/processed/stage5/y.npy')},indent=2))
 plan=json.loads(planpath.read_text());assert plan['source_sha256']==sha(__file__)
 assert sha(root/'data/processed/stage5/x.npy')==plan['input_x_sha256'] and sha(root/'data/processed/stage5/y.npy')==plan['input_y_sha256']
 rows=[];selections=[];splits=[];start=time.monotonic()
 for repeat,seed in enumerate(plan['seeds']):
  for fold,(train,test) in enumerate(StratifiedKFold(3,shuffle=True,random_state=seed).split(x,y)):
   splits.append({'repeat':repeat,'fold':fold,'train':train.tolist(),'test':test.tolist()})
   mask=np.zeros(len(y),bool);mask[train]=True
   for model in MODELS:
    for group in GROUPS:
     stem=f'cv_r{repeat}_f{fold}_{model}_{group}';checkpoint=out/(stem+'.json')
     if checkpoint.exists():record=json.loads(checkpoint.read_text())
     else:
      # Match previous evaluation's inner-fold split seed for paired comparisons.
      best,trials=select(x,y,train,model,group,seed+fold)
      state=fit(x,y,mask,model,group,best,seed+fold*100);s=score(state,x[test])
      joblib.dump(state,out/(stem+'.joblib'));np.testing.assert_array_equal(s,score(joblib.load(out/(stem+'.joblib')),x[test]))
      np.save(out/(stem+'_scores.npy'),s)
      record={'repeat':repeat,'fold':fold,'model':model,'group':group,'best':best,'trials':trials,**metrics(y[test],s)};checkpoint.write_text(json.dumps(record))
     rows.append({k:v for k,v in record.items() if k not in ['best','trials']});selections.append(record)
   print(f'Feature comparison repeat {repeat+1}/5 fold {fold+1}/3; {time.monotonic()-start:.0f}s',flush=True)
   pd.DataFrame(rows).to_csv(out/'feature_fold_metrics.csv',index=False)
 (out/'cv_splits.json').write_text(json.dumps(splits));(out/'cv_selections.json').write_text(json.dumps(selections))
 frame=pd.DataFrame(rows);rep=frame.groupby(['repeat','model','group'])[['ap','recall100','recall500']].mean().reset_index();rep.to_csv(out/'feature_repeat_metrics.csv',index=False)
 rep.groupby(['model','group']).agg(ap_mean=('ap','mean'),ap_min=('ap','min'),ap_max=('ap','max'),recall100=('recall100','mean'),recall500=('recall500','mean')).reset_index().to_csv(out/'feature_summary.csv',index=False)
 # Leave-one-positive-out with fixed strictly held-out unlabelled background.
 held_u=np.random.default_rng(plan['lopo_seed']).choice(u,int(np.ceil(.1*len(u))),replace=False)
 np.save(out/'lopo_heldout_unlabelled.npy',held_u);lopo=[];stability={}
 for i,positive in enumerate(p):
  mask=np.ones(len(y),bool);mask[held_u]=False;mask[positive]=False;train=np.flatnonzero(mask);test=np.r_[positive,held_u]
  for model,group in [('nnpu','full'),('nnpu','mutation'),('gradient_boosting','full'),('lightgbm','full')]:
   stem=f'lopo_{i}_{model}_{group}';cp=out/(stem+'.json')
   if cp.exists():record=json.loads(cp.read_text());s=np.load(out/(stem+'_scores.npy'))
   else:
    best,trials=select(x,y,train,model,group,plan['lopo_seed'])
    state=fit(x,y,mask,model,group,best,plan['lopo_seed']);s=score(state,x[test]);joblib.dump(state,out/(stem+'.joblib'))
    np.testing.assert_array_equal(s,score(joblib.load(out/(stem+'.joblib')),x[test]));np.save(out/(stem+'_scores.npy'),s)
    above=int((s>s[0]).sum());tied=int((s==s[0]).sum())
    record={'symbol':symbols[positive],'gene_id':genes[positive],'model':model,'group':group,'best_rank':above+1,'worst_rank':above+tied,'average_rank':above+(tied+1)/2,'best':best,'trials':trials,**metrics(y[test],s)};cp.write_text(json.dumps(record))
   lopo.append({k:v for k,v in record.items() if k not in ['best','trials']});stability.setdefault(model+'__'+group,[]).append(s[1:])
  print(f'Leave-one-positive-out {i+1}/17 ({symbols[positive]}); {time.monotonic()-start:.0f}s',flush=True)
  pd.DataFrame(lopo).to_csv(out/'lopo_results.csv',index=False)
 stabrows=[];candidates=[]
 for model,arrays in stability.items():
  scores=np.asarray(arrays);weights=np.array([top_weight(s,100) for s in scores]);ranks=np.array([rankdata(-s) for s in scores]);jaccard=[]
  for a,b in itertools.combinations(range(len(scores)),2):jaccard.append(float(np.minimum(weights[a],weights[b]).sum()/np.maximum(weights[a],weights[b]).sum()))
  stabrows.append({'model':model,'background_genes':len(held_u),'perturbations':17,'weighted_jaccard_mean':float(np.mean(jaccard)),'weighted_jaccard_min':float(np.min(jaccard)),'weighted_jaccard_max':float(np.max(jaccard)),'genes_full_top100_all17':int((weights==1).all(axis=0).sum())})
  for j,g in enumerate(held_u):candidates.append({'model':model,'gene_id':genes[g],'symbol':symbols[g],'mean_top100_weight':float(weights[:,j].mean()),'full_top100_runs':int((weights[:,j]==1).sum()),'mean_rank':float(ranks[:,j].mean()),'min_rank':float(ranks[:,j].min()),'max_rank':float(ranks[:,j].max())})
 pd.DataFrame(stabrows).to_csv(out/'stability_summary.csv',index=False);pd.DataFrame(candidates).to_csv(out/'candidate_stability.csv',index=False)
 # Fixed-prior sensitivity on first repeat, each L2 selected internally.
 priorrows=[]
 for fold in [s for s in splits if s['repeat']==0]:
  train=np.asarray(fold['train']);test=np.asarray(fold['test']);mask=np.zeros(len(y),bool);mask[train]=True
  for group in ['mutation','full']:
   for prior in [.0005,.001,.003,.01,.03]:
    stem=f'prior_f{fold["fold"]}_{group}_{prior}';cp=out/(stem+'.json')
    if cp.exists():record=json.loads(cp.read_text())
    else:
     best,trials=select(x,y,train,'nnpu',group,113+fold['fold'],prior)
     state=fit(x,y,mask,'nnpu',group,best,113+100*fold['fold']);s=score(state,x[test]);np.save(out/(stem+'_scores.npy'),s)
     record={'fold':fold['fold'],'group':group,'prior':prior,'best':best,'trials':trials,**metrics(y[test],s)};cp.write_text(json.dumps(record))
    priorrows.append({k:v for k,v in record.items() if k not in ['best','trials']})
  print(f'Prior sensitivity fold {fold["fold"]+1}/3 complete',flush=True)
 pd.DataFrame(priorrows).to_csv(out/'prior_fold_metrics.csv',index=False)
 pd.DataFrame(priorrows).groupby(['group','prior'])[['ap','recall100','recall500']].mean().reset_index().to_csv(out/'prior_summary.csv',index=False)
 (out/'manifest.json').write_text(json.dumps({'source_sha256':sha(__file__),'plan_sha256':sha(planpath),'model_roundtrips':len(list(out.glob('*.joblib'))),'all_exact':True,'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2))
 print('COMPLETED',time.monotonic()-start,flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=1):main(sys.argv[1])
