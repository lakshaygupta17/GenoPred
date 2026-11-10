"""Repeated nested gene CV and prespecified ablations on the original 17 labels."""
import json,sys,hashlib,itertools,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from genopred_v2.stage6_benchmark import load_inputs,fit_entry,predict,prepare,transform_saved
from genopred_v2.stage6 import train_nnpu,_sigmoid,observed_label_metrics

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def options(model):
 if model=='nnpu':return [dict(prior=p,l2=l) for p,l in itertools.product([.001,.003,.01],[.0001,.01])]
 return [dict(level=0),dict(level=1)]

def fit(x,y,mask,names,model,opt,seed,columns=None):
 entry=dict(kind='bag',family=model,level=opt.get('level',0),columns=columns)
 if model=='nnpu':
  values,prep=prepare(x,mask,entry,seed)
  w,b,risk=train_nnpu(values[mask],y[mask],opt['prior'],opt['l2'],.03,600)
  return dict(model=model,prep=prep,weights=w,bias=b)
 return dict(model=model,state=fit_entry(x,y,mask,entry,seed,names))

def score(state,x,names):
 if state['model']=='nnpu':return _sigmoid(transform_saved(x,state['prep'])@state['weights']+state['bias'])
 return predict(state['state'],x,names)

def metric(y,s):
 result=observed_label_metrics(y,s,top_k=(100,500))
 result['ap']=float(average_precision_score(y==1,s))
 for k in [100,500]:
  t=np.sort(s)[-min(k,len(s))];above=s>t;equal=s==t
  result[f'recall_{k}']=float(((y[above]==1).sum()+(min(k,len(s))-above.sum())/equal.sum()*(y[equal]==1).sum())/(y==1).sum())
 return result

def main(root):
 root=Path(root);out=root/'data/processed/final_robustness_v1';out.mkdir(parents=True,exist_ok=True)
 x,y,names,genes,symbols,_,_,_,report=load_inputs(root)
 planpath=root/'configs/final_robustness_v1.json'
 plan={'created_utc':datetime.now(timezone.utc).isoformat(),'seeds':[113,227,349,463,587],'outer_folds':3,'inner_folds':3,'models':['nnpu','decision_tree','gradient_boosting','lightgbm'],'selection':'Mean inner-fold observed-label AP. First candidate wins exact ties. No outer labels/features enter fitting or preprocessing.','grids':{m:options(m) for m in ['nnpu','decision_tree','gradient_boosting','lightgbm']},'PU_bags':6,'unlabelled_per_positive':10,'nnpu_epochs':600,'ablation_models':['nnpu','gradient_boosting'],'ablations':['no_mutation','expression_only','missing_flags_only'],'ablation_rule':'Reuse full-model inner-selected settings; refit on outer training genes only. Diagnostic, not separately optimized.','negative_control':'Shuffle outer-training labels for nnpu and gradient boosting, retaining selected settings.','limitation':'Retrospective robustness study after prior model screening. Not a fresh independent test. Repeats share genes; fold intervals are descriptive, not independent confidence intervals. Unlabelled is not negative.','source_sha256':sha(__file__),'input_files':report['output_files']}
 if planpath.exists():plan=json.loads(planpath.read_text());assert plan['source_sha256']==sha(__file__)
 else:planpath.write_text(json.dumps(plan,indent=2))
 rows=[];selection=[];all_oof=[];split_records=[];saved=[]
 start=time.monotonic()
 for repeat,seed in enumerate(plan['seeds']):
  oof={m:np.full(len(y),np.nan) for m in plan['models']+['mutation_control']}
  for m in plan['ablation_models']:
   for ab in plan['ablations']+['shuffled_labels']:oof[m+'__'+ab]=np.full(len(y),np.nan)
  splits=StratifiedKFold(3,shuffle=True,random_state=seed)
  for fold,(train_idx,test_idx) in enumerate(splits.split(x,y)):
   mask=np.zeros(len(y),bool);mask[train_idx]=True
   assert not mask[test_idx].any()
   split_records.append(dict(repeat=repeat,fold=fold,train_idx=train_idx.tolist(),test_idx=test_idx.tolist()))
   inner=list(StratifiedKFold(3,shuffle=True,random_state=seed+fold+1000).split(train_idx,y[train_idx]))
   for model in plan['models']:
    candidates=[]
    for option in options(model):
     aps=[]
     for inner_fold,(itr,iva) in enumerate(inner):
      imask=np.zeros(len(y),bool);imask[train_idx[itr]]=True
      assert not imask[test_idx].any() and not imask[train_idx[iva]].any()
      state=fit(x,y,imask,names,model,option,seed+fold*100+inner_fold)
      vals=score(state,x[train_idx[iva]],names)
      aps.append(float(average_precision_score(y[train_idx[iva]]==1,vals)))
     candidates.append(dict(option=option,inner_ap=aps,mean=float(np.mean(aps))))
    best=max(candidates,key=lambda r:r['mean'])['option']
    selection.append(dict(repeat=repeat,fold=fold,model=model,best=best,candidates=candidates))
    state=fit(x,y,mask,names,model,best,seed+fold*100)
    prediction=score(state,x[test_idx],names);oof[model][test_idx]=prediction
    path=out/f'repeat{repeat}_fold{fold}_{model}.joblib';joblib.dump(state,path)
    np.testing.assert_array_equal(prediction,score(joblib.load(path),x[test_idx],names));saved.append(path.name)
    rows.append(dict(repeat=repeat,fold=fold,model=model,**metric(y[test_idx],prediction)))
    if model in plan['ablation_models']:
     for ab in plan['ablations']+['shuffled_labels']:
      target=y.copy();group=ab
      if ab=='shuffled_labels':target[train_idx]=np.random.default_rng(seed+fold*100).permutation(target[train_idx]);group=None
      diagnostic=fit(x,target,mask,names,model,best,seed+fold*100,group)
      ss=score(diagnostic,x[test_idx],names);oof[model+'__'+ab][test_idx]=ss
      rows.append(dict(repeat=repeat,fold=fold,model=model+'__'+ab,**metric(y[test_idx],ss)))
   vals=x[test_idx,6].copy();vals[x[test_idx,16]==1]=-1
   oof['mutation_control'][test_idx]=vals
   rows.append(dict(repeat=repeat,fold=fold,model='mutation_control',**metric(y[test_idx],vals)))
   pd.DataFrame(rows).to_csv(out/'fold_metrics.csv',index=False)
   (out/'selection_log.json').write_text(json.dumps(selection,indent=2))
   print(f'Repeat {repeat+1}/5 fold {fold+1}/3 complete; {time.monotonic()-start:.1f}s',flush=True)
  assert all(np.isfinite(v).all() for v in oof.values())
  all_oof.append(oof)
 (out/'splits.json').write_text(json.dumps(split_records))
 model_names=list(all_oof[0]);np.savez_compressed(out/'oof_scores.npz',scores=np.stack([np.column_stack([d[m] for m in model_names]) for d in all_oof]),model_names=model_names,gene_ids=genes,symbols=symbols,labels=y)
 frame=pd.DataFrame(rows)
 # Each repeat averages three independent outer folds; repeat variability remains descriptive.
 repeat_scores=frame.groupby(['repeat','model'])[['ap','recall_100','recall_500']].mean().reset_index()
 repeat_scores.to_csv(out/'repeat_metrics.csv',index=False)
 summary=repeat_scores.groupby('model').agg(ap_mean=('ap','mean'),ap_min=('ap','min'),ap_max=('ap','max'),recall100_mean=('recall_100','mean'),recall500_mean=('recall_500','mean')).reset_index()
 summary.to_csv(out/'summary.csv',index=False)
 manifest={'plan_sha256':sha(planpath),'source_sha256':sha(__file__),'model_roundtrips':len(saved),'all_exact':True,'elapsed_seconds':time.monotonic()-start,'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(summary.to_string(index=False),flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=1):main(sys.argv[1])
