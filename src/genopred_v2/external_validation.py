"""Frozen, orthogonal Project Score assessment; no fitting or model selection."""
from pathlib import Path
import argparse
import hashlib
import json
import warnings
import zipfile
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return result


def ranks(x):
    a = rankdata(x, axis=0).astype(float)
    a -= a.mean(axis=0)
    norm = np.sqrt((a * a).sum(axis=0))
    return a / np.where(norm > 0, norm, np.nan)


def top_ties(scores, k):
    return scores >= np.sort(scores)[-min(k, len(scores))]


def collapse_donors(values, donors):
    # genes x screens -> genes x independent donors
    donors = np.asarray(donors)
    return np.column_stack([np.nanmedian(values[:, donors == d], axis=1) for d in sorted(set(donors))])


def run(root):
    root = Path(root)
    raw = root / 'data/raw/validation/project_score_2025'
    out = root / 'data/processed/external_validation_v1'
    protocol_path = root / 'configs/external_validation_v1.json'
    protocol = json.loads(protocol_path.read_text())
    schema = json.loads((root / 'configs/external_validation_v1_schema.json').read_text())
    manifest = json.loads((raw / 'source_manifest.json').read_text())
    assert sha(protocol_path) == manifest['protocol_sha256'] == schema['protocol_sha256']
    for name, item in manifest['files'].items():
        assert (raw / name).stat().st_size == item['bytes'] and sha(raw / name) == item['sha256'], name
    frozen = out / 'frozen_predictions.npz'
    assert sha(frozen) == protocol['frozen_predictions_sha256']
    pred = np.load(frozen, allow_pickle=False)
    names = pred['model_names'].astype(str).tolist()
    assert names[:4] == protocol['shortlist']
    print('Verified all source and frozen prediction hashes', flush=True)
    meta = pd.read_csv(raw / 'model_list_20261013.csv', dtype=str).fillna('')
    assert not meta.model_id.duplicated().any()
    meta = meta.set_index('model_id')
    with zipfile.ZipFile(next(raw.glob('*.zip'))) as z:
        member = next(n for n in z.namelist() if 'fold_change_values' in n)
        with z.open(member) as f:
            headers = [f.readline().decode().rstrip('\r\n').split('\t') for _ in range(5)]
        with z.open(member) as f:
            frame = pd.read_csv(f, sep='\t', skiprows=5, header=None)
    ids = frame[2].astype(str).str.split('.').str[0]
    keep = ids.str.fullmatch(r'ENSG\d+') & ~ids.duplicated(keep=False)
    matrix = frame.loc[keep].iloc[:, 3:].to_numpy(float)
    ids = ids[keep].to_numpy()
    screens = pd.DataFrame({'model_id': headers[1][3:], 'source': headers[2][3:], 'qc': headers[3][3:]})
    assert len(screens) == matrix.shape[1]
    screens['column'] = np.arange(len(screens))
    screens = screens[screens.qc.str.upper().eq('TRUE') & screens.model_id.isin(meta.index)].copy()
    # Match source documentation: retain Sanger when both sources screen one model.
    screens['priority'] = screens.source.map({'Sanger': 0, 'Broad': 1}).fillna(2)
    screens = screens.sort_values(['model_id', 'priority', 'column']).drop_duplicates('model_id')
    screens = screens.join(meta[['patient_id','model_name','synonyms','tissue','cancer_type_detail','model_type']], on='model_id')
    normal_names = (meta.model_name + '|' + meta.synonyms).str.upper().str.replace(r'[^A-Z0-9|]', '', regex=True)
    excluded = normal_names.str.contains(r'(?:^|\|)PANC1(?:$|\|)|HPNE', regex=True)
    excluded_ids = set(meta.index[excluded])
    excluded_donors = set(meta.loc[excluded, 'patient_id']) - {''}
    screens['excluded_overlap'] = screens.model_id.isin(excluded_ids) | screens.patient_id.isin(excluded_donors)
    pancreatic = screens.tissue.str.lower().eq('pancreas')
    is_adeno = screens.cancer_type_detail.str.strip().isin(schema['adenocarcinoma_details'])
    known_donor = screens.patient_id.ne('')
    screens['cohort'] = 'excluded'
    screens.loc[pancreatic & is_adeno & known_donor & ~screens.excluded_overlap, 'cohort'] = 'pancreatic'
    screens.loc[~pancreatic & ~screens.tissue.str.strip().str.lower().isin(['','unknown']) & known_donor & ~screens.excluded_overlap, 'cohort'] = 'other'
    # Ensure cross-lineage annotations cannot cause donor overlap.
    p_donors = set(screens.loc[screens.cohort.eq('pancreatic'), 'patient_id'])
    screens.loc[screens.cohort.eq('other') & screens.patient_id.isin(p_donors), 'cohort'] = 'excluded'
    screens.drop(columns=['synonyms']).to_csv(out / 'cohort_audit.csv', index=False)
    groups = {}
    for cohort in ['pancreatic', 'other']:
        rows = screens[screens.cohort.eq(cohort)]
        groups[cohort] = collapse_donors(matrix[:, rows.column.to_numpy()], rows.patient_id.to_numpy())
    p, o = groups['pancreatic'], groups['other']
    assert p.shape[1] >= protocol['minimum_pancreatic_donors']
    assert o.shape[1] >= protocol['minimum_other_donors']
    print(f'Eligible donors: pancreas={p.shape[1]}, other={o.shape[1]}', flush=True)
    annotation = pd.read_csv(raw / 'CRISPR_KO_supplemental_information_20250624.csv')
    annotation['ensembl_gene_id'] = annotation.ensembl_gene_id.astype(str).str.split('.').str[0]
    annotation = annotation[~annotation.ensembl_gene_id.duplicated(keep=False)].set_index('ensembl_gene_id')
    prediction_ids = pd.Index([str(g).split('.')[0] for g in pred['gene_ids']])
    assert not prediction_ids.duplicated().any()
    index = pd.Index(ids).get_indexer(prediction_ids)
    matched = index >= 0
    p = p[index[matched]]; o = o[index[matched]]
    scores = pred['scores'][matched].astype(float)
    genes = pd.DataFrame({'gene_id': prediction_ids[matched], 'symbol': pred['gene_symbols'][matched], 'label': pred['labels'][matched]})
    a = annotation.reindex(genes.gene_id)
    common = a.common_essential.eq(True).to_numpy() | a.core_fitness_pancan.eq(True).to_numpy()
    annotated = a.common_essential.notna().to_numpy() & a.core_fitness_pancan.notna().to_numpy()
    pcov = np.isfinite(p).mean(axis=1); ocov = np.isfinite(o).mean(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        pm = np.nanmedian(p, axis=1); om = np.nanmedian(o, axis=1)
    eligible = (pcov >= protocol['coverage_fraction']) & np.isfinite(pm) & np.isfinite(scores).all(axis=1)
    assert set(np.unique(pred['labels'])) == {-1, 1}
    primary = eligible & annotated & ~common & (genes.label.to_numpy() == -1)
    assert primary.sum() > 500, 'Too few eligible genes for the registered analysis'
    specific = primary & (ocov >= protocol['coverage_fraction']) & np.isfinite(om)
    genes['common_essential_or_core'] = common
    genes['annotation_available'] = annotated
    genes['pancreatic_coverage'] = pcov; genes['other_coverage'] = ocov
    genes['median_pancreatic_fold_change'] = pm; genes['median_other_fold_change'] = om
    genes['selectivity'] = om - pm; genes['primary_eligible'] = primary
    for j, name in enumerate(names): genes[name] = scores[:, j]
    genes.to_csv(out / 'gene_evidence.csv.gz', index=False)
    genes[genes.label.eq(1)].to_csv(out / 'known_positive_controls.csv', index=False)
    x = ranks(scores[primary]); y = ranks(-pm[primary])
    rho = x.T @ y
    rng = np.random.default_rng(protocol['seed'])
    counts = np.zeros(len(names), int)
    for _ in range(protocol['permutations']): counts += (x.T @ rng.permutation(y)) >= rho
    pv = (counts + 1) / (protocol['permutations'] + 1)
    adjusted = holm(pv[:4])
    boot = []
    for _ in range(protocol['bootstrap_donors']):
        sample = rng.integers(0, p.shape[1], p.shape[1])
        bp = np.nanmedian(p[primary][:, sample], axis=1)
        valid = np.isfinite(bp)
        boot.append(ranks(scores[primary][valid]).T @ ranks(-bp[valid]))
    ci = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    selectivity = ranks(scores[specific]).T @ ranks((om - pm)[specific])
    all_rho = ranks(scores[eligible]).T @ ranks(-pm[eligible])
    unlabelled_all = eligible & (genes.label.to_numpy() == -1)
    common_included_rho = ranks(scores[unlabelled_all]).T @ ranks(-pm[unlabelled_all])
    rows = []
    for j, name in enumerate(names):
        rows.append(dict(model=name, primary_genes=int(primary.sum()), rho=float(rho[j]), permutation_p=float(pv[j]), holm_p=float(adjusted[j]) if j < 4 else None, donor_bootstrap_low=float(ci[0,j]), donor_bootstrap_high=float(ci[1,j]), selectivity_rho=float(selectivity[j]), all_measured_rho=float(all_rho[j]), unlabelled_including_common_rho=float(common_included_rho[j])))
    metrics = pd.DataFrame(rows); metrics.to_csv(out / 'model_comparison.csv', index=False)
    top_rows = []
    candidates = []
    v = pm[primary]
    for j, name in enumerate(names):
        for k in [100,500]:
            chosen = top_ties(scores[primary,j], k)
            actual = int(chosen.sum()); observed = float(np.median(v[chosen]))
            null = np.array([np.median(v[rng.choice(len(v), actual, replace=False)]) for _ in range(protocol['permutations'])])
            top_rows.append(dict(model=name, requested_top_k=k, actual_genes_with_ties=actual, median_fold_change=observed, enrichment_p=float((1+(null<=observed).sum())/(len(null)+1)), null_median=float(np.median(null))))
            if k == 100:
                part = genes[primary].loc[chosen].copy(); part['selected_by'] = name; candidates.append(part)
    pd.DataFrame(top_rows).to_csv(out / 'top_gene_enrichment.csv', index=False)
    pd.concat(candidates).to_csv(out / 'top_candidate_evidence.csv', index=False)
    report = dict(protocol_sha256=sha(protocol_path), schema_sha256=sha(root/'configs/external_validation_v1_schema.json'), code_sha256=sha(__file__), frozen_predictions_sha256=sha(frozen), pancreatic_donors=p.shape[1], other_donors=o.shape[1], matched_genes=int(matched.sum()), primary_genes=int(primary.sum()), selectivity_genes=int(specific.sum()), excluded_overlap_models=screens.loc[screens.excluded_overlap, 'model_name'].tolist(), results=rows, top_enrichment=top_rows)
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))
    print(metrics.to_string(index=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('root', type=Path)
    with threadpool_limits(limits=1): run(parser.parse_args().root)
