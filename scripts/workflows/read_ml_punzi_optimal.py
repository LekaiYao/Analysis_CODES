"""Read ML-owned Punzi selections; never scan, optimize or fit."""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
from fit_input_contract import resolve_selected_dir, load_training


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expression_tree(expression):
    return ast.parse(expression.replace('&&', ' and ').replace('||', ' or ').strip(), mode='eval').body


def selection_signature(expression):
    def canonical(node):
        if isinstance(node, ast.BoolOp):
            kind=type(node.op).__name__
            def terms(n):
                if isinstance(n, ast.BoolOp) and type(n.op) is type(node.op):
                    return [v for child in n.values for v in terms(child)]
                return [canonical(n)]
            return kind, tuple(sorted(set(v for child in node.values for v in terms(child)),key=repr))
        return ast.dump(node, include_attributes=False)
    return canonical(expression_tree(expression))


def root_selection(expression):
    """Translate a deliberately restricted selection grammar to TTreeFormula."""
    def emit(n):
        if isinstance(n, ast.Name): return n.id
        if isinstance(n, ast.Constant) and type(n.value) in (int,float): return repr(n.value)
        if isinstance(n, ast.BoolOp):
            sep=' && ' if isinstance(n.op,ast.And) else ' || '
            return '('+sep.join(emit(v) for v in n.values)+')'
        if isinstance(n, ast.Compare):
            operators={ast.Gt:'>',ast.Lt:'<',ast.GtE:'>=',ast.LtE:'<=',ast.Eq:'==',ast.NotEq:'!='}
            values=[n.left]+n.comparators
            return '('+' && '.join(emit(values[i])+' '+operators[type(op)]+' '+emit(values[i+1]) for i,op in enumerate(n.ops))+')'
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='abs' and len(n.args)==1 and not n.keywords:
            return 'abs('+emit(n.args[0])+')'
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd,ast.Not)):
            return '('+{ast.USub:'-',ast.UAdd:'+',ast.Not:'!'}[type(n.op)]+emit(n.operand)+')'
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub,ast.Mult,ast.Div)):
            return '('+emit(n.left)+{ast.Add:'+',ast.Sub:'-',ast.Mult:'*',ast.Div:'/'}[type(n.op)]+emit(n.right)+')'
        raise ValueError('Unsupported selection syntax: '+ast.dump(n))
    return emit(expression_tree(expression))


def read_ml_optimal(repo, tag, scope, *, year="2023"):
    if scope not in ('training','wide'): raise ValueError('scope must explicitly be training or wide')
    repo=Path(repo).resolve();ml=repo.parent/'XGBoost';base=resolve_selected_dir(repo,tag)
    scan_path=base/'punzi/scan.json';recipe_path=base/'punzi/recipe.json'
    scan=json.loads(scan_path.read_text());recipe=json.loads(recipe_path.read_text())
    if scan.get('status')!='PASS': raise ValueError('ML Punzi is not PASS')
    if scan.get('train_tag')!=tag: raise ValueError('ML Punzi tag mismatch')
    if digest(recipe_path)!=scan.get('recipe_sha256'): raise ValueError('ML recipe hash mismatch')
    if recipe.get('recipe_id')!=scan.get('recipe_id'): raise ValueError('ML recipe ID mismatch')
    chosen=scan.get('scopes',{}).get(scope)
    if chosen is None: raise ValueError('Requested ML Punzi scope is missing: '+scope)
    optimal=chosen['optimal'];score=optimal['score']
    if type(score) not in (float,int) or not math.isfinite(score) or not 0<=score<=1: raise ValueError('Invalid ML score')
    # Match the delivered selection to its own stored row; do not recompute the maximum.
    rows=[r for r in chosen['curves']['total'] if r['score']==score]
    if len(rows)!=1: raise ValueError('Chosen score must have exactly one total-curve row')
    row=rows[0]
    if not math.isfinite(optimal['FOM']) or row['FOM']!=optimal['FOM']: raise ValueError('Chosen FOM and stored row disagree')
    epsilon=row['epsilon']
    if not math.isfinite(epsilon) or not 0<=epsilon<=1: raise ValueError('Invalid stored MC efficiency')
    rel=Path(chosen['mass_plot'])
    if rel.is_absolute(): raise ValueError('Expected mass_plot relative to ML root')
    plot=(ml/rel).resolve()
    if not plot.is_relative_to(base/'cut_scan'): raise ValueError('Mass plot is outside this tag cut_scan')
    if not plot.is_file(): raise FileNotFoundError('ML optimal mass plot not ready: '+str(plot))
    t,_=load_training(repo,tag,year=year)
    for key,expression in [('data_apply_selection',t['selection']),('mc_apply_selection',t.get('mc_apply_selection',t['selection']))]:
        if selection_signature(scan[key])!=selection_signature(expression): raise ValueError('Punzi/APPLY selection mismatch: '+key)
    model=ml/'output/models'/base.relative_to((ml/'output/selected').resolve())
    files={'apply_summary_sha256':base/'batch_apply_summary.json','model_sha256':model/'xgb_model.pkl','scaler_sha256':model/'scaler.pkl','config_sha256':model/'model_config.json','metadata_sha256':model/'run_metadata.json'}
    for key,path in files.items():
        if digest(path)!=scan.get(key): raise ValueError('Stale ML provenance: '+key)
    for item in t['sources'].values():
        if not Path(item['path']).is_file(): raise FileNotFoundError('Scored input missing: '+item['path'])
    return {'schema_version':1,'train_tag':tag,'scope':scope,'score_branch':'Prediction','comparison':'>','score_threshold':score,
            'data_apply_selection':scan['data_apply_selection'],'mc_apply_selection':scan['mc_apply_selection'],
            'background_training_selection':scan['background_training_selection'],
            'root_data_apply_selection':root_selection(scan['data_apply_selection']),
            'root_mc_apply_selection':root_selection(scan['mc_apply_selection']),
            'root_score_selection':'Prediction > '+repr(score),
            'sidebands':chosen['sidebands'],'width_ratio':chosen['width_ratio'],'optimal':optimal,
            'stored_mc_efficiency':epsilon,'efficiency_definition':'tag prompt X, 3.62<Bmass<4.00, source pThatreweight; read from ML total row',
            'mass_plot':str(plot),'plot_config':scan['plot_config'],'scan_json':str(scan_path),'scan_sha256':digest(scan_path),
            'recipe_json':str(recipe_path),'recipe_sha256':digest(recipe_path),'recipe_id':scan['recipe_id'],
            'scored_inputs':t['sources'],'source_provenance':scan['sources'],
            'verified_artifacts':{key:{'path':str(path),'sha256':scan[key]} for key,path in files.items()},
            'source_root_hashes_recomputed':False,'selection_origin':'ML scopes.'+scope+'.optimal; no downstream optimization',
            'fit_authorized_by_reader':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--scope',required=True,choices=['training','wide'])
    args=parser.parse_args()
    try: result=read_ml_optimal(Path(__file__).resolve().parents[2],args.tag,args.scope)
    except (OSError,ValueError,KeyError,TypeError,SyntaxError) as exc: parser.exit(2,'ML optimal not ready or invalid: '+str(exc)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))

if __name__=='__main__': main()
