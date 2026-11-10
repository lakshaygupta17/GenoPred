"""Broad, bounded exploratory benchmark on the frozen Stage 6 gene split.

Run with python -m genopred_v2.stage6_benchmark ROOT.
The registry is fixed before fitting. No test scores choose settings or ensemble members.
Joblib bundles are trusted local files only. Graph bundles require the exact pinned graph.
"""
from __future__ import annotations

import csv
import gzip
import inspect
import itertools
import json
import platform
import sys
import time
import warnings
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer, PolynomialFeatures
from sklearn.decomposition import PCA
from sklearn.kernel_approximation import RBFSampler
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.svm import LinearSVC, SVC, OneClassSVM
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier,
    GradientBoostingClassifier, HistGradientBoostingClassifier, AdaBoostClassifier, IsolationForest)
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors, KernelDensity
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.utils.class_weight import compute_sample_weight
from threadpoolctl import threadpool_limits

from .stage6 import (Stage6Error, _digest, _write_tsv_gz, _sigmoid,
    deterministic_split, impute_from_partition, observed_label_metrics,
    diffuse_scores, train_nnpu, predict_saved_model, predict_saved_graph, SPLIT_NAMES)
from .stage6_trees import predict_bundle as predict_old_tree

SEEDS = [41, 137, 907]
BAGS = 6
RATIO = 10
PU_FAMILIES = ['decision_tree', 'random_forest', 'extra_trees', 'gradient_boosting',
    'hist_gradient_boosting', 'adaboost', 'xgboost', 'lightgbm', 'catboost',
    'logistic_l2', 'logistic_l1', 'ridge', 'linear_svm', 'rbf_svm', 'knn',
    'gaussian_nb', 'lda', 'qda', 'mlp', 'gaussian_process']


def registry():
    entries = []
    def add(name, kind, candidates=None, **extra):
        entries.append(dict(name=name, kind=kind, candidates=candidates or [dict(level=0)], **extra))
    for family in PU_FAMILIES:
        add('pu_'+family, 'bag', [dict(level=0),dict(level=1)], family=family)
    for family in ['logistic_l2', 'random_forest']:
        add('biased_'+family, 'biased', [dict(level=0),dict(level=1)], family=family)
    add('two_step_provisional_negatives', 'two_step', [dict(level=0),dict(level=1)])
    for family in ['centroid', 'nearest_positive', 'positive_kde', 'positive_ocsvm',
                   'positive_isolation_forest']:
        add(family, 'positive_only', [dict(level=0),dict(level=1)], family=family)
    add('rbf_features_nnpu', 'nnpu', [dict(level=0),dict(level=1)])
    add('positive_seed_graph_diffusion', 'graph', [dict(level=0),dict(level=1)])
    add('graph_degree', 'degree')
    for i in range(10):
        add(f'feature_{i:02d}', 'feature', [dict(level=0),dict(level=1)], index=i)
    for group in ['mutation_only','cnv_only','expression_only','regulation_only',
                  'no_mutation','biology_without_flags','missing_flags_only']:
        add('ablation_'+group,'bag',family='logistic_l2',columns=group)
    for prep in ['robust','quantile','pca','polynomial']:
        add('preprocess_'+prep,'bag',family='logistic_l2',preprocess=prep)
    add('random_rank_control', 'random')
    add('constant_score_control', 'constant')
    add('shuffled_training_labels_control','shuffle',family='logistic_l2')
    return entries


def columns_for(group, width):
    if width != 20:
        raise Stage6Error('Registry requires the frozen 10 biological plus 10 missingness features')
    groups = dict(mutation_only=[6,16], cnv_only=[7,8,17,18], expression_only=[9,19],
        regulation_only=list(range(6))+list(range(10,16)),
        no_mutation=[i for i in range(20) if i not in (6,16)],
        biology_without_flags=list(range(10)), missing_flags_only=list(range(10,20)))
    return groups[group] if group else list(range(width))


def prepare(x, fitting_mask, entry, seed):
    clean, medians = impute_from_partition(x, fitting_mask, 10)
    columns = columns_for(entry.get('columns'), x.shape[1])
    values = clean[:,columns]
    mode = entry.get('preprocess', 'standard')
    scaler = (RobustScaler() if mode=='robust' else
              QuantileTransformer(n_quantiles=min(100,int(fitting_mask.sum())),
                  output_distribution='normal',random_state=seed) if mode=='quantile' else StandardScaler())
    scaler.fit(values[fitting_mask]); values=scaler.transform(values)
    transform = None
    if mode=='pca':
        transform=PCA(n_components=5,svd_solver='full')
    elif mode=='polynomial':
        transform=PolynomialFeatures(degree=2,interaction_only=True,include_bias=False)
    elif entry['kind']=='nnpu':
        transform=RBFSampler(gamma=.05 if entry['level']==0 else .2,n_components=64,random_state=seed)
    if transform is not None:
        transform.fit(values[fitting_mask]); values=transform.transform(values)
    return values, dict(medians=medians,columns=columns,scaler=scaler,transform=transform)


def transform_saved(x, prep):
    values=np.asarray(x,dtype=float).copy()
    for i,median in enumerate(prep['medians']):
        values[values[:,10+i]==1,i]=median
    values=prep['scaler'].transform(values[:,prep['columns']])
    if prep['transform'] is not None:
        values=prep['transform'].transform(values)
    return values


def estimator(family, level, seed):
    depth=[2,3][level]; strength=[.1,1.][level]
    if family=='decision_tree': return DecisionTreeClassifier(max_depth=depth,min_samples_leaf=2,random_state=seed)
    if family in ('random_forest','extra_trees'):
        cls=RandomForestClassifier if family=='random_forest' else ExtraTreesClassifier
        return cls(n_estimators=40,max_depth=depth,min_samples_leaf=2,n_jobs=1,random_state=seed)
    if family=='gradient_boosting': return GradientBoostingClassifier(n_estimators=60,max_depth=depth,learning_rate=.05,random_state=seed)
    if family=='hist_gradient_boosting': return HistGradientBoostingClassifier(max_iter=60,max_leaf_nodes=2**depth,min_samples_leaf=3,learning_rate=.05,early_stopping=False,random_state=seed)
    if family=='adaboost': return AdaBoostClassifier(estimator=DecisionTreeClassifier(max_depth=1),n_estimators=60,learning_rate=[.1,.5][level],random_state=seed)
    if family=='xgboost':
        from xgboost import XGBClassifier
        return XGBClassifier(n_estimators=60,max_depth=depth,learning_rate=.05,reg_lambda=5,n_jobs=1,tree_method='hist',random_state=seed)
    if family=='lightgbm':
        from lightgbm import LGBMClassifier
        return LGBMClassifier(n_estimators=60,max_depth=depth,num_leaves=2**depth,min_child_samples=3,learning_rate=.05,reg_lambda=5,n_jobs=1,verbosity=-1,random_state=seed)
    if family=='catboost':
        from catboost import CatBoostClassifier
        return CatBoostClassifier(iterations=60,depth=depth,learning_rate=.05,l2_leaf_reg=5,thread_count=1,verbose=False,allow_writing_files=False,random_seed=seed)
    if family.startswith('logistic'): return LogisticRegression(C=strength,penalty='l1' if family.endswith('l1') else 'l2',solver='liblinear',max_iter=2000,random_state=seed)
    if family=='ridge': return RidgeClassifier(alpha=[10.,1.][level])
    if family=='linear_svm': return LinearSVC(C=strength,dual=False,max_iter=5000,random_state=seed)
    if family=='rbf_svm': return SVC(C=strength,kernel='rbf',gamma='scale',probability=False,random_state=seed)
    if family=='knn': return KNeighborsClassifier(n_neighbors=[3,7][level],weights='distance',n_jobs=1)
    if family=='gaussian_nb': return GaussianNB(var_smoothing=[1e-3,.1][level])
    if family=='lda': return LinearDiscriminantAnalysis(solver='lsqr',shrinkage=[.5,.9][level],priors=[.5,.5])
    if family=='qda': return QuadraticDiscriminantAnalysis(reg_param=[.5,.9][level],priors=[.5,.5])
    if family=='mlp': return MLPClassifier(hidden_layer_sizes=(8,),alpha=[1.,10.][level],solver='lbfgs',max_iter=800,random_state=seed)
    if family=='gaussian_process': return GaussianProcessClassifier(kernel=RBF([1.,3.][level]),optimizer=None,max_iter_predict=60,random_state=seed)
    raise Stage6Error(f'Unknown estimator {family}')


def fit_classifier(model, x, y):
    if 'sample_weight' in inspect.signature(model.fit).parameters:
        model.fit(x,y,sample_weight=compute_sample_weight('balanced',y))
    else:
        model.fit(x,y)
    return model


def model_score(model, x):
    # Chunking bounds Gaussian-process and nearest-neighbour prediction memory.
    result=[]
    for start in range(0,len(x),2048):
        part=x[start:start+2048]
        if hasattr(model,'predict_proba'):
            result.append(model.predict_proba(part)[:,list(model.classes_).index(1)])
        else:
            result.append(_sigmoid(model.decision_function(part)))
    return np.concatenate(result)


def fit_entry(x,y,mask,entry,seed,names):
    """All label access in this function is limited to mask; graph labels are masked too."""
    kind=entry['kind']; level=entry['level']
    values, prep=prepare(x,mask,entry,seed)
    xx, yy=values[mask],y[mask].copy()
    rng=np.random.default_rng(seed)
    if kind=='shuffle': yy=rng.permutation(yy)
    state=dict(entry=entry,seed=seed,prep=prep,feature_names=list(names))
    p,u=np.flatnonzero(yy==1),np.flatnonzero(yy==-1)
    if not len(p) or not len(u): raise Stage6Error('Missing training PU group')
    if kind in ('bag','shuffle'):
        models=[]
        for _ in range(BAGS):
            ids=np.r_[p,rng.choice(u,min(len(u),RATIO*len(p)),replace=False)]
            model=estimator(entry['family'],level,int(rng.integers(0,2**31-1)))
            models.append(fit_classifier(model,xx[ids],(yy[ids]==1).astype(int)))
        state['models']=models
    elif kind=='biased':
        state['models']=[fit_classifier(estimator(entry['family'],level,seed),xx,(yy==1).astype(int))]
    elif kind=='two_step':
        pilot=fit_classifier(estimator('logistic_l2',0,seed),xx,(yy==1).astype(int))
        comparator_scores=model_score(pilot,xx[u])
        cutoff=np.quantile(comparator_scores,[.1,.3][level])
        rn=u[comparator_scores<=cutoff]
        ids=np.r_[p,rn]
        state['provisional_negative_count']=len(rn)
        state['models']=[fit_classifier(estimator('random_forest',0,seed),xx[ids],(yy[ids]==1).astype(int))]
    elif kind=='positive_only':
        positive=xx[p]; family=entry['family']
        if family=='centroid': state['center']=positive.mean(axis=0)
        elif family=='nearest_positive': state['positive_model']=NearestNeighbors(n_neighbors=[1,3][level]).fit(positive)
        elif family=='positive_kde': state['positive_model']=KernelDensity(bandwidth=[1.,3.][level]).fit(positive)
        elif family=='positive_ocsvm': state['positive_model']=OneClassSVM(nu=[.1,.3][level],gamma='scale').fit(positive)
        else: state['positive_model']=IsolationForest(n_estimators=100,max_samples='auto',random_state=seed,n_jobs=1).fit(positive)
    elif kind=='nnpu':
        w,b,risk=train_nnpu(xx,yy,.001,.01,.03,300)
        state.update(weights=w,bias=b,final_risk=risk)
    elif kind=='graph':
        state['initial']=np.where(mask & (y==1),1.,0.)
    elif kind in ('degree','feature','random','constant'): pass
    else: raise Stage6Error(f'Unknown experiment {kind}')
    return state


def predict(state,x,names,edges=None,weights=None):
    if list(names)!=state['feature_names'] or x.shape[1]!=len(names):
        raise Stage6Error('Feature order or dimensions mismatch')
    entry=state['entry']; kind=entry['kind']; level=entry['level']
    values=transform_saved(x,state['prep'])
    if 'models' in state:
        scores=np.mean([model_score(m,values) for m in state['models']],axis=0)
    elif kind=='positive_only':
        family=entry['family']
        if family=='centroid':
            if level==0: scores=-np.linalg.norm(values-state['center'],axis=1)
            else:
                scores=values@state['center']/(np.maximum(np.linalg.norm(values,axis=1)*np.linalg.norm(state['center']),1e-12))
        elif family=='nearest_positive': scores=-state['positive_model'].kneighbors(values)[0].mean(axis=1)
        elif family=='positive_kde': scores=state['positive_model'].score_samples(values)
        else: scores=state['positive_model'].decision_function(values)
    elif kind=='nnpu': scores=_sigmoid(values@state['weights']+state['bias'])
    elif kind=='graph':
        scores,_=diffuse_scores(state['initial'],edges,weights,[.3,.7][level],100,1e-10,1.)
    elif kind=='degree':
        scores=np.bincount(edges[0],weights=np.asarray(weights).ravel(),minlength=len(x))
    elif kind=='feature':
        # Use raw observed values, with missing values ranked last. Missingness is not a low value.
        index=entry['index']; scores=x[:,index].astype(float).copy()*([1.,-1.][level])
        missing=x[:,index+10]==1
        scores[missing]=(scores[~missing].min()-max(1.,np.ptp(scores[~missing]))) if (~missing).any() else 0.
    elif kind=='random': scores=np.random.default_rng(state['seed']).random(len(x))
    elif kind=='constant': scores=np.zeros(len(x))
    else: raise Stage6Error('Unknown prediction mode')
    if scores.shape!=(len(x),) or not np.isfinite(scores).all():
        raise Stage6Error('Invalid predictions')
    return scores


def predict_bundle(bundle,x,names,edges=None,weights=None):
    if bundle['gene_count']!=len(x) and bundle['transductive']:
        raise Stage6Error('Transductive bundle requires the frozen gene universe')
    return np.mean([predict(s,x,names,edges,weights) for s in bundle['states']],axis=0)


def load_inputs(root):
    source=root/'data/processed/stage5'
    report=json.loads((root/'manifests/stage5_dataset_report.json').read_text())
    if report['dataset_version']!='stage5-core-v2': raise Stage6Error('Corrected data required')
    for filename,meta in report['output_files'].items():
        if _digest(source/filename)!=meta['sha256']: raise Stage6Error(f'Hash mismatch: {filename}')
    x,y=[np.load(source/f'{n}.npy',allow_pickle=False) for n in ('x','y')]
    names=json.loads((source/'feature_names.json').read_text())
    genes=(source/'gene_ids.txt').read_text().splitlines()
    with gzip.open(source/'gene_features.tsv.gz','rt',encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    if [r['ensembl_gene_id'] for r in rows]!=genes: raise Stage6Error('Gene order mismatch')
    if x.shape!=(len(genes),20) or len(y)!=len(genes) or len(names)!=20: raise Stage6Error('Input shape mismatch')
    cfg=json.loads((root/'configs/stage6.json').read_text())
    split=deterministic_split(genes,y,cfg['seed'],cfg['split']['positive_counts'],cfg['split']['unlabelled_fractions'])
    np.testing.assert_array_equal(split,np.load(root/'data/processed/stage6/split.npy',allow_pickle=False))
    edges=np.load(source/'edge_index.npy',allow_pickle=False)
    weights=np.load(source/'edge_attr.npy',allow_pickle=False)
    return x,y,names,genes,[r['gene_symbol'] for r in rows],split,edges,weights,report


def write_json(path, data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(root):
    start=time.monotonic()
    x,y,names,genes,symbols,split,edges,weights,report5=load_inputs(root)
    out=root/'data/processed/stage6_benchmark_v1'
    out.mkdir(parents=True,exist_ok=True)
    entries=registry()
    plan=dict(version='stage6-broad-benchmark-v1',created_at_utc=datetime.now(UTC).isoformat(),
        seeds=SEEDS,bags_per_seed=BAGS,unlabelled_per_positive=RATIO,experiments=entries,
        selection='Mean validation AP over fixed seeds; tie chooses first lower-complexity candidate',
        rank_ensemble='Top three bagging families by validation AP only; equal percentile-rank average',
        post_selection_diagnostics='Three development-only stratified folds for top three bagging families; fixed settings; no retuning',
        source_sha256=_digest(Path(__file__)),input_files=report5['output_files'])
    write_json(out/'plan.json',plan)
    train,val,test=split==0,split==1,split==2
    candidates,selected,failures,warning_log=[],{},[],[]
    validation_scores={}
    for entry in entries:
        print('SELECT '+entry['name'],flush=True)
        family_candidates=[]
        for option in entry['candidates']:
            spec={k:v for k,v in entry.items() if k!='candidates'}|option
            try:
                per_seed=[]; metrics=[]
                for seed in SEEDS:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        state=fit_entry(x,y,train,spec,seed,names)
                        scores=predict(state,x,names,edges,weights)
                    warning_log.extend(dict(experiment=entry['name'],phase='selection',message=str(w.message)) for w in caught)
                    per_seed.append(scores)
                    metrics.append(observed_label_metrics(y[val],scores[val]))
                row=dict(spec=spec,validation=metrics,
                    mean_validation_ap=float(np.mean([m['observed_label_auprc'] for m in metrics])))
                candidates.append(row); family_candidates.append((row,np.mean(per_seed,axis=0)))
            except Exception as exc:
                failures.append(dict(experiment=entry['name'],option=option,phase='selection',error=repr(exc)))
        if family_candidates:
            best,vs=max(family_candidates,key=lambda pair:pair[0]['mean_validation_ap'])
            selected[entry['name']]=best; validation_scores[entry['name']]=vs
        write_json(out/'selection_progress.json',dict(selected=selected,candidates=candidates,failures=failures))
    eligible=[n for n in selected if n.startswith('pu_')]
    members=sorted(eligible,key=lambda n:selected[n]['mean_validation_ap'],reverse=True)[:3]
    write_json(out/'selection_frozen.json',dict(selected=selected,ensemble_members=members,failures=failures,
        frozen_at_utc=datetime.now(UTC).isoformat()))
    # No test metrics were calculated above this line.
    results={}; final_scores={}; final_paths={}
    for name,best in selected.items():
        print('REFIT '+name,flush=True)
        try:
            states=[]; seed_scores=[]
            for seed in SEEDS:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    state=fit_entry(x,y,~test,best['spec'],seed,names)
                    scores=predict(state,x,names,edges,weights)
                warning_log.extend(dict(experiment=name,phase='final',message=str(w.message)) for w in caught)
                states.append(state); seed_scores.append(scores)
            prediction=np.mean(seed_scores,axis=0)
            bundle=dict(states=states,feature_names=names,gene_ids=genes,gene_count=len(genes),
                transductive=best['spec']['kind'] in ('graph','degree','random'),
                graph_input_hashes={k:v for k,v in report5['output_files'].items() if k in ('edge_index.npy','edge_attr.npy')})
            path=out/f'{name}.joblib'; joblib.dump(bundle,path,compress=3)
            restored=predict_bundle(joblib.load(path),x,names,edges,weights)
            error=float(np.max(np.abs(restored-prediction)))
            if error>1e-12: raise Stage6Error('Saved-model round trip failed')
            per_seed_metrics=[observed_label_metrics(y[test],s[test]) for s in seed_scores]
            correlations=[]
            for a,b in itertools.combinations(seed_scores,2):
                if np.std(a[test]) and np.std(b[test]): correlations.append(float(spearmanr(a[test],b[test]).statistic))
            results[name]=dict(validation_mean_ap=best['mean_validation_ap'],
                test=observed_label_metrics(y[test],prediction[test]),seed_test_metrics=per_seed_metrics,
                mean_test_rank_correlation=float(np.mean(correlations)) if correlations else None,
                roundtrip_max_error=error)
            final_scores[name]=prediction; final_paths[name]=path.name
        except Exception as exc:
            failures.append(dict(experiment=name,phase='final',error=repr(exc)))
        write_json(out/'results_progress.json',dict(results=results,failures=failures))
    # The ensemble membership is already frozen. Ranking uses features/scores, never test labels.
    if len(members)==3 and all(m in final_scores for m in members):
        ensemble=np.mean([rankdata(final_scores[m],method='average')/len(x) for m in members],axis=0)
        vs=np.mean([rankdata(validation_scores[m],method='average')/len(x) for m in members],axis=0)
        ensemble_spec=dict(members=members,model_paths=[final_paths[m] for m in members],
            aggregation='mean average-tie percentile rank over frozen full gene universe',gene_ids=genes)
        write_json(out/'rank_ensemble.json',ensemble_spec)
        reloaded=[]
        for path in ensemble_spec['model_paths']:
            reloaded.append(rankdata(predict_bundle(joblib.load(out/path),x,names,edges,weights))/len(x))
        np.testing.assert_array_equal(ensemble,np.mean(reloaded,axis=0))
        final_scores['validation_selected_rank_ensemble']=ensemble
        results['validation_selected_rank_ensemble']=dict(test=observed_label_metrics(y[test],ensemble[test]),
            validation_mean_ap=observed_label_metrics(y[val],vs[val])['observed_label_auprc'],
            roundtrip_max_error=0.,members=members)
    # Earlier models are reference runs, not candidates in the new selection search.
    baseline_path=root/'data/processed/stage6/model.npz'
    refs=dict(reference_nnpu=predict_saved_model(baseline_path,x,names),
              reference_nnpu_graph=predict_saved_graph(baseline_path,x,names,edges,weights))
    for name in ('pu_tree','pu_random_forest'):
        refs['reference_'+name]=predict_old_tree(joblib.load(root/f'data/processed/stage6_trees/{name}.joblib'),x,names)
    for name,scores in refs.items():
        results[name]=dict(test=observed_label_metrics(y[test],scores[test]),reference=True)
        final_scores[name]=scores
    # Fixed-setting development cross-validation: sensitivity diagnostic, not fresh validation.
    from sklearn.model_selection import StratifiedKFold
    dev=np.flatnonzero(~test); folds=StratifiedKFold(n_splits=3,shuffle=True,random_state=20261104)
    cv={}
    for name in members:
        print('DEVELOPMENT_CV '+name,flush=True)
        rows=[]
        for fit_ids,eval_ids in folds.split(dev,y[dev]):
            mask=np.zeros(len(x),dtype=bool); mask[dev[fit_ids]]=True
            with warnings.catch_warnings(record=True) as caught:
                state=fit_entry(x,y,mask,selected[name]['spec'],SEEDS[0],names)
                scores=predict(state,x,names,edges,weights)
            warning_log.extend(dict(experiment=name,phase='development_cv',message=str(w.message)) for w in caught)
            rows.append(observed_label_metrics(y[dev[eval_ids]],scores[dev[eval_ids]]))
        cv[name]=rows
    score_names=list(final_scores)
    matrix=np.column_stack([final_scores[n] for n in score_names])
    np.savez_compressed(out/'all_scores.npz',scores=matrix,model_names=np.asarray(score_names),
                        gene_ids=np.asarray(genes),gene_symbols=np.asarray(symbols),labels=y,split=split)
    rows=[]
    for i,gene in enumerate(genes):
        rows.append(dict(ensembl_gene_id=gene,gene_symbol=symbols[i],label=int(y[i]),split=SPLIT_NAMES[int(split[i])],
                         **{n:float(final_scores[n][i]) for n in score_names}))
    _write_tsv_gz(out/'all_gene_scores.tsv.gz',rows,list(rows[0]))
    positive_ranks=[]
    for name,scores in final_scores.items():
        subset=scores[test]; indices=np.flatnonzero(test)
        for i in indices[y[indices]==1]:
            value=scores[i]
            positive_ranks.append(dict(model=name,gene=symbols[i],rank_min=int(1+np.sum(subset>value)),
                rank_max=int(np.sum(subset>=value)),score=float(value)))
    write_json(out/'test_positive_ranks.json',positive_ranks)
    versions={package:version(package) for package in ('numpy','scipy','scikit-learn','joblib','xgboost','lightgbm','catboost')}
    write_json(out/'environment.json',dict(python=platform.python_version(),packages=versions))
    (out/'requirements-lock.txt').write_text(''.join(f'{p}=={v}\n' for p,v in versions.items()))
    # Retain distinct warning messages and multiplicities instead of flooding the report.
    warning_counts={}
    for item in warning_log:
        key=(item['experiment'],item['phase'],item['message'])
        warning_counts[key]=warning_counts.get(key,0)+1
    write_json(out/'warnings.json',[dict(experiment=k[0],phase=k[1],message=k[2],count=v) for k,v in warning_counts.items()])
    report=dict(experiment_version='stage6-broad-benchmark-v1',completed_at_utc=datetime.now(UTC).isoformat(),
        seconds=time.monotonic()-start,plan_sha256=_digest(out/'plan.json'),
        selection_sha256=_digest(out/'selection_frozen.json'),source_sha256=_digest(Path(__file__)),
        evaluation_status='Exploratory model search on reused validation and test split; multiple-comparison optimism',
        selected=selected,results=results,failures=failures,development_cv=cv,environment=versions,
        original_baseline_report_sha256=_digest(root/'manifests/stage6_baseline_report.json'),
        original_tree_report_sha256=_digest(root/'manifests/stage6_tree_report.json'),
        input_files=report5['output_files'],output_files={p.name:dict(bytes=p.stat().st_size,sha256=_digest(p))
             for p in out.iterdir() if p.is_file()})
    write_json(root/'manifests/stage6_benchmark_report.json',report)
    table=['# Broad Stage 6 benchmark results','',
        'All results are exploratory on the reused 3,992-gene test subset with three known positives.',
        'Rows are sorted by test average precision for inspection, not used to select a deployment model.','',
        '| Experiment | Validation AP | Test AP | Test AUROC | Recall at 100 | Recall at 500 |',
        '|---|---:|---:|---:|---:|---:|']
    for name,item in sorted(results.items(),key=lambda pair:pair[1]['test']['observed_label_auprc'],reverse=True):
        m=item['test']; va=item.get('validation_mean_ap')
        va_text=f'{va:.5f}' if va is not None else 'prior run'
        table.append(f"| {name} | {va_text} | {m['observed_label_auprc']:.5f} | {m['observed_label_auroc']:.5f} | {m['known_positive_recall_at_100']:.3f} | {m['known_positive_recall_at_500']:.3f} |")
    table+=['','## Failed attempts','',json.dumps(failures,indent=2),'',
            'See STAGE6_EXPERIMENT_CATALOG.md for scope, methods, limitations and reproduction.']
    (root/'docs/STAGE6_BENCHMARK_RESULTS.md').write_text('\n'.join(table)+'\n',encoding='utf-8')
    print(json.dumps(dict(completed=len(results),failed_attempts=len(failures),ensemble_members=members,
                         seconds=report['seconds']),indent=2),flush=True)
    return report


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main(Path(sys.argv[1]).resolve())
