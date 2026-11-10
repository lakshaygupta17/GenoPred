"""Exploratory PU tree ensembles on the preserved Stage 6 split.

Run: python -m genopred_v2.stage6_trees PROJECT_ROOT
Saved joblib files are trusted local artifacts only; never load untrusted pickle files.
"""
import itertools
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import sklearn
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

from .stage6 import (
    Stage6Error, _digest, _write_tsv_gz, deterministic_split,
    impute_from_partition, observed_label_metrics, SPLIT_NAMES,
)

CONFIG = {
    'version': 'stage6-pu-trees-v1', 'seeds': [41, 137, 907],
    'depths': [2, 3], 'unlabelled_per_positive': [5, 20],
    'bags': 40, 'forest_trees_per_bag': 10, 'min_samples_leaf': 2,
    'selection': 'mean validation observed-label average precision across all seeds',
}


def fit_bags(x, y, family, depth, ratio, seed, bags=40):
    if set(np.unique(y)) != {-1, 1}:
        raise Stage6Error('PU bags require positive and unlabelled training examples')
    rng = np.random.default_rng(seed)
    pos, un = np.flatnonzero(y == 1), np.flatnonzero(y == -1)
    models = []
    for _ in range(bags):
        ids = np.r_[pos, rng.choice(un, min(len(un), ratio * len(pos)), replace=False)]
        options = dict(max_depth=depth, min_samples_leaf=2, class_weight='balanced',
                       random_state=int(rng.integers(0, 2**31-1)))
        model = (DecisionTreeClassifier(**options) if family == 'pu_tree' else
                 RandomForestClassifier(n_estimators=10, n_jobs=1, **options))
        model.fit(x[ids], (y[ids] == 1).astype(int))
        models.append(model)
    return models


def score(models, x):
    return np.mean([m.predict_proba(x)[:, list(m.classes_).index(1)] for m in models], axis=0)


def predict_bundle(bundle, x, feature_names):
    if list(feature_names) != bundle['feature_names']:
        raise Stage6Error('Feature order mismatch')
    clean = np.asarray(x, dtype=float).copy()
    count = len(bundle['medians'])
    for i, median in enumerate(bundle['medians']):
        clean[clean[:, count+i] == 1, i] = median
    return score(bundle['models'], clean)


def run(root):
    source = root / 'data/processed/stage5'
    baseline = json.loads((root/'manifests/stage6_baseline_report.json').read_text())
    report5 = json.loads((root/'manifests/stage5_dataset_report.json').read_text())
    if report5['dataset_version'] != 'stage5-core-v2':
        raise Stage6Error('Corrected Stage 5 required')
    for name, meta in report5['output_files'].items():
        if _digest(source/name) != meta['sha256']:
            raise Stage6Error(f'Input hash mismatch: {name}')
    cfg = json.loads((root/'configs/stage6.json').read_text())
    x, y = [np.load(source/f'{name}.npy', allow_pickle=False) for name in ['x', 'y']]
    genes = (source/'gene_ids.txt').read_text().splitlines()
    names = json.loads((source/'feature_names.json').read_text())
    split = deterministic_split(genes, y, cfg['seed'], cfg['split']['positive_counts'],
                                cfg['split']['unlabelled_fractions'])
    np.testing.assert_array_equal(split, np.load(root/'data/processed/stage6/split.npy'))
    train, val, test = split == 0, split == 1, split == 2
    selection_x, _ = impute_from_partition(x, train, x.shape[1]//2)
    final_x, medians = impute_from_partition(x, train | val, x.shape[1]//2)
    out = root/'data/processed/stage6_trees'
    out.mkdir(parents=True, exist_ok=True)
    config_path = root/'configs/stage6_trees.json'
    config_path.write_text(json.dumps(CONFIG, indent=2)+'\n')
    candidates, chosen = [], {}
    for family in ['pu_tree', 'pu_random_forest']:
        for depth, ratio in itertools.product(CONFIG['depths'], CONFIG['unlabelled_per_positive']):
            metrics = []
            for seed in CONFIG['seeds']:
                models = fit_bags(selection_x[train], y[train], family, depth, ratio, seed)
                metrics.append(observed_label_metrics(y[val], score(models, selection_x[val])))
            candidates.append(dict(family=family, depth=depth, ratio=ratio, validation=metrics,
                mean_validation_ap=float(np.mean([m['observed_label_auprc'] for m in metrics]))))
        chosen[family] = max([c for c in candidates if c['family']==family],
                             key=lambda c:(c['mean_validation_ap'], -c['depth'], -c['ratio']))
    # Persist all choices before any test evaluation. Never select a seed by its test result.
    (out/'selection.json').write_text(json.dumps(chosen, indent=2)+'\n')
    results, predictions, errors = {}, {}, {}
    for family, best in chosen.items():
        models, per_seed = [], []
        for seed in CONFIG['seeds']:
            fitted = fit_bags(final_x[~test], y[~test], family, best['depth'], best['ratio'], seed)
            models.extend(fitted)
            per_seed.append(score(fitted, final_x))
        prediction = score(models, final_x)
        bundle = dict(models=models, medians=medians, feature_names=names,
                      sklearn_version=sklearn.__version__, configuration=CONFIG, selected=best)
        path = out/f'{family}.joblib'
        joblib.dump(bundle, path)
        restored = predict_bundle(joblib.load(path), x, names)
        np.testing.assert_array_equal(prediction, restored)
        errors[family] = float(np.max(np.abs(prediction-restored)))
        predictions[family] = prediction
        results[family] = dict(test=observed_label_metrics(y[test], prediction[test]),
            per_seed_test=[observed_label_metrics(y[test], p[test]) for p in per_seed],
            mean_gene_score_std_across_seeds=float(np.std(per_seed, axis=0).mean()))
    rows = [dict(ensembl_gene_id=g, label=int(y[i]), split=SPLIT_NAMES[int(split[i])],
                 **{k:float(v[i]) for k,v in predictions.items()}) for i,g in enumerate(genes)]
    _write_tsv_gz(out/'gene_scores.tsv.gz', rows, list(rows[0]))
    report = dict(configuration=CONFIG, sklearn_version=sklearn.__version__,
        numpy_version=np.__version__, input_files=report5['output_files'],
        source_sha256=_digest(Path(__file__)), config_sha256=_digest(config_path),
        baseline_report_sha256=_digest(root/'manifests/stage6_baseline_report.json'),
        validation_candidates=candidates, selected=chosen, results=results,
        baseline_test_metrics=baseline['test_metrics'], roundtrip_max_error=errors,
        evaluation_status='Exploratory comparison on previously examined test split; 3 positives',
        score_meaning='Uncalibrated PU bagging ranking score; sampled unlabelled examples are temporary comparators',
        output_files={p.name:dict(sha256=_digest(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    (root/'manifests/stage6_tree_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(results, indent=2))
    return report


if __name__ == '__main__':
    run(Path(sys.argv[1]).resolve())
