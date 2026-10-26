"""Shared input validation; no fitting or output side effects."""
import ast, hashlib, json, re
from pathlib import Path
MASS_RANGE = [3.75, 4.00]
POINTS = [15, 20, 25, 30, 35, 40]
REVISION = 'pb23_single_x_wide_v1'
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def selection_expression(expression):
    class BooleanArray(ast.NodeTransformer):
        def visit_BoolOp(self, node):
            values = [self.visit(v) for v in node.values]
            result = values[0]
            for v in values[1:]: result = ast.BinOp(result, ast.BitAnd() if isinstance(node.op, ast.And) else ast.BitOr(), v)
            return ast.copy_location(result, node)
    tree = ast.parse(expression.replace('&&', ' and ').replace('||', ' or '), mode='eval')
    return ast.unparse(ast.fix_missing_locations(BooleanArray().visit(tree)))
def resolve_selected_dir(repo, tag):
    """Read ML's authoritative owner registry; never pick a glob's first match.

    The 2026-09-20 single-year delivery explicitly defines this path interface.
    Read the literal registry without importing/executing upstream Python code.
    """
    if not re.fullmatch(r'[A-Za-z0-9_]+', tag):
        raise ValueError('Invalid training tag: '+tag)
    ml = repo.parent/'XGBoost'
    registry = ml/'configs/feature_selection.py'
    owner = None
    if registry.exists():
        assignments = [n for n in ast.parse(registry.read_text()).body
                       if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'FEATURE_SELECTION_EXPERIMENTS' for t in n.targets)]
        if len(assignments) != 1:
            raise ValueError('Unsupported ML owner registry format')
        experiments = ast.literal_eval(assignments[0].value)
        owners = [name for name, cfg in experiments.items()
                  if any(item['train_tag'] == tag for item in cfg['candidates'])]
        if len(owners) > 1:
            raise ValueError(f'Ambiguous ML owner for {tag}: {owners}')
        if owners:
            owner = owners[0]; cfg = experiments[owner]
            tags = [item['train_tag'] for item in cfg['candidates']]
            if len(tags) != len(set(tags)) or cfg['baseline'] in tags:
                raise ValueError('Invalid duplicate/baseline candidates in '+owner)
            varsets = re.findall(r'_(\d+v\d+)_', cfg['baseline'])
            if len(varsets) != 1 or cfg['operation'] not in ('add1', 'drop1') or owner != f"fs_{varsets[0]}_{cfg['operation']}_{cfg['baseline']}":
                raise ValueError('Invalid ML owner directory: '+owner)
    root = ml/'output/selected'
    base = root/owner/tag if owner else root/tag
    if owner and (root/tag).exists():
        raise ValueError(f'Ambiguous legacy and registered selected directories for {tag}')
    if not base.is_dir():
        raise FileNotFoundError(f'Missing or archived ML input for {tag}; expected owner path: {base}')
    return base.resolve()

def load_training(repo, tag, *, year="2023"):
    base = resolve_selected_dir(repo, tag)
    tp = base/'cut_scan/thresholds.json'
    if not tp.exists(): tp = base/'cut_scan/weighted_signal_efficiency/thresholds.json'
    ts = json.loads(tp.read_text()); weight = ts['weight_branch']
    if isinstance(weight, list):
        assert len(weight)==1; weight=weight[0]
    rec = dict(tag=tag, selected_dir=str(base), thresholds_file=str(tp), thresholds_sha256=sha(tp), weight=weight)
    mp=base/'fit_scan_manifest.pb23_pb24_simultaneous_mc_shape_nominal_v2.json'
    if mp.exists():
        if str(year) != "2023": raise ValueError("Legacy paired-year manifest is PbPb23-only in this reader")
        manifest=json.loads(mp.read_text());cat=manifest['pairing']['categories']['pb23']
        assert cat['signal_mc']['event_weight_branch']==weight
        assert cat['threshold_provenance']['sha256']==sha(tp)
        rec.update(manifest=str(mp),manifest_sha256=sha(mp),selection=cat['fiducial_selection']['expression'])
        specs={'data':cat['data'],'mc':cat['signal_mc']}
        for point in ts['thresholds']:
            matches=[p for p in manifest['working_points'] if abs(p['target_weighted_efficiency']-point['target_efficiency'])<1e-8]
            if matches: assert matches[0]['categories']['pb23']['threshold']==point['score_threshold']
    else:
        bs=base/'batch_apply_summary.json';b=json.loads(bs.read_text())
        assert str(b['input_datasets']['dataset_year'])==str(year), 'Dataset year mismatch'
        rec.update(batch_apply_summary_sha256=sha(bs),selection=b.get('apply_selection', {}).get('expression', b['draw_selection']['fiducial_cut']['expression']))
        assert rec['selection'] == b['draw_selection']['fiducial_cut']['expression'], 'Apply/draw fiducial mismatch'
        rec['mc_apply_selection'] = b.get('apply_selection', {}).get('mc_expression', rec['selection'])
        specs={'data':{'path':'DATA_with_score.root','tree':'ntmix'},'mc':{'path':'MC_with_score.root','tree':'ntmix_X3872'}}
    rec['sources']={k:{'path':str((base/v['path']).resolve()),'tree':v['tree']} for k,v in specs.items()}
    return rec,ts['thresholds']
def verify_resume(prior, current, points, width, mean):
    assert prior.get('workflow_revision')==REVISION, 'Old workflow; use a new semantic variant, not resume'
    assert prior['mass_range']==MASS_RANGE and prior['target_efficiencies_percent']==sorted(round(p['target_efficiency']*100) for p in points)
    assert prior['width_scale_range']==width and prior['mean_half_range']==mean
    previous=next(t for t in prior['trainings'] if t['tag']==current['tag'])
    for key in ['thresholds_sha256','weight','selection','sources','manifest_sha256','batch_apply_summary_sha256']:
        assert previous.get(key)==current.get(key), 'Resume input contract changed: '+key
    assert previous['points']==points
    for spec in previous['files'].values(): assert sha(spec['path'])==spec['sha256']
    for key in ['data','mc']:
        info=previous[key+'_source'];stat=Path(info['path']).stat()
        assert stat.st_size==info['size'] and stat.st_mtime_ns==info['mtime_ns'], 'Scored source changed'
    return previous
