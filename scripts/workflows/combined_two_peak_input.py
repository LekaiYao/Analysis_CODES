"""Versioned PbPb18/PbPb23 combined-fit input contract; no ROOT or writes."""
import json
import math
from pathlib import Path
from fit_input_contract import load_training, sha
from read_ml_punzi_optimal import read_ml_optimal, selection_signature

REVISION = 'pb18_pb23_two_peak_total_fraction_v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def ml_signature(config, metadata):
    # Sample-dependent class balancing is intentionally computed per year.
    params = {k: v for k, v in config['model_params'].items()
              if k not in ('scale_pos_weight', 'n_jobs')}
    require(metadata['class_weighting']['method'] == 'scale_pos_weight',
            'Unsupported class weighting')
    require(config['signal_weight_branch'] == 'pThatreweight', 'Expected source MC weights')
    require(config.get('signal_weight_branches') == ['pThatreweight'], 'Composite weights unsupported')
    require(config['model_params']['scale_pos_weight'] == metadata['class_weighting']['value'],
            'Class-weight metadata mismatch')
    return {**{k: config[k] for k in (
        'input_columns', 'trans_columns', 'direct_xgb_version', 'reweight_profile',
        'signal_weight_expression', 'split_policy', 'split_random_state',
        'train_fraction', 'scaler_fit_scope', 'training_signal_exclusion')},
        'model_params': params, 'class_weighting_method': metadata['class_weighting']['method'],
        'signal_selection': selection_signature(metadata['signal_selection']),
        'background_selection': selection_signature(metadata['background_selection'])}


def category(repo, tag, year, *, scope=None, xeff=None):
    require(tag.startswith('X_pb'+str(year)[2:]+'_'), 'Tag/year mismatch')
    t, points = load_training(repo, tag, year=str(year))
    require(t['weight'] == 'pThatreweight', 'Expected pThatreweight')
    require(json.loads(Path(t['thresholds_file']).read_text())['train_tag'] == tag, 'Threshold tag mismatch')
    base = Path(t['selected_dir'])
    ml = repo.parent/'XGBoost'
    model = ml/'output/models'/base.relative_to((ml/'output/selected').resolve())
    config_path, metadata_path = model/'model_config.json', model/'run_metadata.json'
    config, metadata = (json.loads(p.read_text()) for p in (config_path, metadata_path))
    require(metadata['train_tag'] == tag, 'Model tag mismatch')
    require(config['input_columns'] == metadata['input_columns'], 'Feature metadata mismatch')
    if scope:
        choice = read_ml_optimal(repo, tag, scope, year=str(year))
        score = choice['score_threshold']
    else:
        require(xeff is not None and 0 < xeff < 100, 'Specify an efficiency in (0,100)')
        matches = [p for p in points if abs(p['target_efficiency']-xeff/100) < 1.e-8]
        require(len(matches) == 1, 'Efficiency point is missing or ambiguous')
        score = matches[0]['score_threshold']
        choice = {'kind': 'weighted_efficiency', 'target_percent': xeff,
                  'point': matches[0], 'thresholds_file': t['thresholds_file'],
                  'thresholds_sha256': t['thresholds_sha256']}
    require(math.isfinite(score) and 0 <= score <= 1, 'Invalid score threshold')
    psi_path = base/'MC_psi2s_apply_manifest.json'
    psi = json.loads(psi_path.read_text())
    require(psi['status'] == 'PASS' and psi['train_tag'] == tag, 'Psi2S delivery not PASS for tag')
    require(psi['score_branch'] == 'Prediction' and psi['weight_branch'] == 'pThatreweight',
            'Psi2S score/weight contract mismatch')
    require(psi['score_cut'] is None, 'Psi2S already score-cut')
    require(Path(psi['output']).resolve().parent == base, 'Psi2S output is outside tag directory')
    require(selection_signature(psi['selection']) == selection_signature(t['mc_apply_selection']),
            'Psi2S/X APPLY selections differ')
    sources = {k: dict(t['sources'][source]) for k, source in [('data', 'data'), ('x', 'mc')]}
    sources['psi'] = {'path': str(Path(psi['output']).resolve()), 'tree': psi['tree'],
                      'expected_sha256': psi['sha256'], 'expected_entries': psi['entries']}
    for spec in sources.values():
        stat = Path(spec['path']).stat()
        spec.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    provenance = {str(p.resolve()): sha(p) for p in (
        config_path, metadata_path, psi_path, base/'batch_apply_summary.json',
        model/'xgb_model.pkl', model/'scaler.pkl', Path(t['thresholds_file']))}
    if scope:
        provenance.update({choice['scan_json']: choice['scan_sha256'],
                           choice['recipe_json']: choice['recipe_sha256']})
    return {'year': str(year), 'tag': tag, 'score': score, 'score_branch': 'Prediction',
            'comparison': '>', 'selection': t['selection'],
            'mc_selection': t['mc_apply_selection'], 'choice': choice,
            'sources': sources, 'provenance': provenance,
            'ml_signature': ml_signature(config, metadata),
            'class_weight': config['model_params']['scale_pos_weight']}


def add_upstream_bchi2prob(signature):
    """Add only the documented redundant top-level conjunction, preserving all other terms."""
    require(isinstance(signature, tuple) and signature[0] == 'And',
            'Expected a conjunctive training/APPLY selection')
    atom = selection_signature('Bchi2Prob > 0.1')
    require(atom not in signature[1], '2018 metadata already contains Bchi2Prob threshold')
    return ('And', tuple(sorted(set(signature[1]) | {atom}, key=repr)))


def verify_bchi2prob_alignment(cats):
    """Verify the sole metadata difference and all delivered entries; do not change selections."""
    import numpy as np
    import uproot
    left = dict(cats['2018']['ml_signature'])
    right = cats['2023']['ml_signature']
    for key in ('signal_selection', 'background_selection'):
        left[key] = add_upstream_bchi2prob(left[key])
    require(left == right, 'Alignment allows only the redundant Bchi2Prob > 0.1 difference')
    for key in ('selection', 'mc_selection'):
        require(add_upstream_bchi2prob(selection_signature(cats['2018'][key])) ==
                selection_signature(cats['2023'][key]), 'Other APPLY differences are not allowed')
    checks = {}
    for year, cat in cats.items():
        for key, spec in cat['sources'].items():
            with uproot.open(spec['path']) as f:
                a = f[spec['tree']]['Bchi2Prob'].array(library='np')
            require(len(a) > 0 and bool(np.isfinite(a).all()) and bool((a > .1).all()),
                    'Bchi2Prob > 0.1 is not redundant on all '+year+'/'+key+' entries')
            stat = Path(spec['path']).stat()
            require((stat.st_size, stat.st_mtime_ns) == (spec['size'], spec['mtime_ns']),
                    'Input changed during Bchi2Prob verification')
            checks[year+'/'+key] = {'entries': len(a), 'minimum': float(a.min()), 'all_above_0p1': True}
    return {'kind': 'verified_upstream_bchi2prob_gt_0p1', 'checks': checks,
            'selection_changed': False, 'metadata_difference': '2018 upstream cut; 2023 explicit cut'}


def build_contract(repo, tag18, tag23, *, scope=None, xeff18=None, xeff23=None, align_bchi2prob=False):
    require((scope is not None) != (xeff18 is not None or xeff23 is not None),
            'Choose one scope or two explicit efficiencies')
    if scope is None:
        require(xeff18 is not None and xeff23 is not None, 'Both year efficiencies are required')
    cats = {str(y): category(repo, tag, y, scope=scope, xeff=eff)
            for y, tag, eff in [(2018, tag18, xeff18), (2023, tag23, xeff23)]}
    alignment = verify_bchi2prob_alignment(cats) if align_bchi2prob else None
    if not align_bchi2prob:
        require(cats['2018']['ml_signature'] == cats['2023']['ml_signature'],
                'Independent trainings must use the same ML configuration and training selections')
        for key in ('selection', 'mc_selection'):
            require(selection_signature(cats['2018'][key]) == selection_signature(cats['2023'][key]),
                    'Year-dependent selection mismatch: '+key)
    model_hashes = [next(v for p, v in c['provenance'].items() if p.endswith('/xgb_model.pkl'))
                    for c in cats.values()]
    require(model_hashes[0] != model_hashes[1], 'Expected independently trained year models')
    return {'schema_version': 1, 'workflow_revision': REVISION, 'categories': cats,
            'mass_range': [3.62, 4.00], 'background_order': 2,
            'fraction_definition': 'f_k = N_k_2023 / (N_k_2018 + N_k_2023)',
            'shape_strategy': 'year/species independent MC double Gaussian; scale=1; independent DATA means',
            'conditional_floating': {'H1': ['x_total', 'bkg_total'], 'H0': ['bkg_total']},
            'statistical_interpretation': 'same-data fixed-parameter conditional sqrt(q0); uncalibrated',
            'presentation_label': 'Significance',
            'method_description': 'Conditional likelihood test; other parameters fixed to the initial combined fit',
            'same_ml_configuration': True, 'selection_equivalence': alignment,
            'class_weight_policy': 'same scale_pos_weight rule; sample-derived values may differ by year'}
