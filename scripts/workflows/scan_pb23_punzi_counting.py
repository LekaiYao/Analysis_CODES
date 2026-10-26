#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Count-based a=5 Punzi from the delivered single-year ML input contract."""
from pathlib import Path
import argparse, csv, hashlib, json, math, re, shutil, subprocess
import numpy as np
import numexpr
import uproot
from fit_input_contract import load_training, selection_expression
from plot_pb23_punzi_a5 import plot


def fingerprint(path):
    path = Path(path); before = path.stat(); h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''): h.update(block)
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    return dict(path=str(path), sha256=h.hexdigest(), bytes=before.st_size, mtime_ns=before.st_mtime_ns)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args(); repo = Path(__file__).resolve().parents[2]
    out = args.output.resolve(); assert not out.exists(), 'Refuse to overwrite results'
    training, registered = load_training(repo, args.tag)
    assert training['weight'] == 'pThatreweight'
    base = Path(training['selected_dir'])
    batch = json.loads((base/'batch_apply_summary.json').read_text())
    header = json.loads(Path(training['thresholds_file']).read_text())
    assert header['train_tag'] == args.tag
    inputs = {k:fingerprint(base/n) for k,n in [('batch','batch_apply_summary.json'),('thresholds',str(Path(training['thresholds_file']).relative_to(base)))]}
    specs = {**training['sources'], 'reference':dict(path=str(base/'REFERENCE_MC_with_score.root'), tree='ntmix_X3872')}
    arrays = {}; counts = {}; weight_values = {}
    for key, spec in specs.items():
        inputs[key] = fingerprint(spec['path'])
        selection = training['selection'] if key == 'data' else training.get('mc_apply_selection',training['selection'])
        branches = sorted(set(['Bmass','Prediction'] + re.findall(r'\bB\w+\b',selection) + ([] if key == 'data' else ['pThatreweight'])))
        with uproot.open(spec['path']) as f: raw = f[spec['tree']].arrays(branches,library='np')
        assert all(np.isfinite(v).all() for v in raw.values())
        mask = numexpr.evaluate(selection_expression(selection),local_dict=raw)
        assert mask.all(), 'Delivered apply selection mismatch: '+key
        assert np.all((raw['Prediction'] >= 0) & (raw['Prediction'] <= 1))
        arrays[key] = {n:np.asarray(raw[n],dtype=float) for n in ['Bmass','Prediction'] + ([] if key == 'data' else ['pThatreweight'])}
        counts[key] = len(mask)
        if key != 'data':
            assert (raw['pThatreweight'] > 0).all()
            values, entries = np.unique(raw['pThatreweight'],return_counts=True)
            weight_values[key] = [dict(weight=float(w),entries=int(n)) for w,n in zip(values,entries)]
        if key in ('data','mc'): assert counts[key] == batch['apply_selection'][key+'_entries']
    mc,data,ref = (arrays[k] for k in ('mc','data','reference'))
    rs,rw = ref['Prediction'],ref['pThatreweight']; order=np.argsort(rs,kind='stable'); cum=np.cumsum(rw[order])/rw.sum()
    def threshold(e): return float(rs[order[min(int(np.searchsorted(cum,1-e,side='left')),len(rs)-1)]])
    for row in registered:
        c=threshold(row['target_efficiency']); assert c == row['score_threshold']
        assert abs(float(rw[rs>c].sum()/rw.sum())-row['achieved_efficiency']) < 1e-12
        if 'data_entries' in row: assert int(np.count_nonzero((data['Prediction']>c)&(data['Bmass']>3.62)&(data['Bmass']<4.00))) == row['data_entries']
    md=(mc['Bmass']>3.62)&(mc['Bmass']<4.00); total=float(mc['pThatreweight'][md].sum()); assert total>0
    left=(data['Bmass']>3.75)&(data['Bmass']<3.84); right=(data['Bmass']>3.90)&(data['Bmass']<3.99)
    out.mkdir(parents=True); art=out/'artifacts';art.mkdir()
    for name in ['scan_pb23_punzi_counting.py','fit_input_contract.py','plot_pb23_punzi_a5.py']: shutil.copy2(repo/'scripts/workflows'/name,art/name)
    shutil.copy2(training['thresholds_file'],art/'thresholds.json');shutil.copy2(base/'batch_apply_summary.json',art/'batch_apply_summary.json')
    contract=dict(tag=args.tag,plot_label=args.tag.replace('X_pb23_','PbPb23 ').replace('_fid20_',' ').replace('_rw0_xgb_',' ')+': counting',metric='MC weighted fraction in D / (2.5 + sqrt((N_left+N_right)/3))',a_primary=5,mass_range=[3.62,4.00],window=[3.84,3.90],signal_window_factor=1.,sideband_windows=[[3.75,3.84],[3.90,3.99]],background_width_ratio=1/3,background_method='symmetric sideband counts; no fit',comparison='strict mass boundaries; Prediction > threshold',selection=training['selection'],mc_apply_selection=training.get('mc_apply_selection',training['selection']),weight='pThatreweight',weight_values=weight_values,mc_sumw_D_no_score_cut=total,xeff_grid=list(range(15,61)),score_grid=[i/100 for i in range(5,96,5)],xeff_coordinate='achieved weighted REFERENCE MC efficiency without mass cut',training=training,inputs=inputs,input_counts=counts,error_definition='sideband Poisson only; no MC statistics or model systematic',reference_mean_used=False,reference_yield_used=False,bootstrap=0,new_data_peak_fits=0)
    dump(out/'contract.json',contract); summary={}; allrows=[]
    for grid,points in [('xeff',[(e,threshold(e/100)) for e in range(15,61)]),('score',[(None,i/100) for i in range(5,96,5)])]:
        rows=[]
        for e,c in points:
            score=data['Prediction']>c; nl=int(np.count_nonzero(left&score)); nr=int(np.count_nonzero(right&score)); n=nl+nr; assert n>0
            b=n/3; error=math.sqrt(n)/3; sw=float(mc['pThatreweight'][md&(mc['Prediction']>c)].sum()); fraction=sw/total;den=2.5+math.sqrt(b)
            row=dict(xeff=e,score_threshold=c,xeff_actual_percent=100*float(rw[rs>c].sum()/rw.sum()),mc_sumw_D=sw,mc_fraction_D=fraction,window_acceptance=1.,signal_window_fraction=fraction,left_entries=nl,right_entries=nr,sideband_entries=n,background=b,background_error=error,punzi_a5=fraction/den,punzi_a5_background_error=fraction*error/(2*math.sqrt(b)*den**2),status='PASS')
            masses=np.sort(data['Bmass'][score]); independent=[int(np.searchsorted(masses,h,side='left')-np.searchsorted(masses,l,side='right')) for l,h in [(3.75,3.84),(3.90,3.99)]]
            assert independent==[nl,nr]
            fsum=math.fsum(mc['pThatreweight'][md&(mc['Prediction']>c)])/math.fsum(mc['pThatreweight'][md]);ib=sum(independent)*(3.90-3.84)/((3.84-3.75)+(3.99-3.90))
            assert abs(fsum/(2.5+math.sqrt(ib))-row['punzi_a5'])<1e-12
            fun=lambda x:fsum/(2.5+math.sqrt(x));step=ib*1e-5
            assert abs(abs((fun(ib+step)-fun(ib-step))/(2*step))*error-row['punzi_a5_background_error'])<1e-10
            rows.append(row)
        dump(out/f'{grid}_scan.json',rows)
        with (out/f'{grid}_scan.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        best=max(rows,key=lambda r:r['punzi_a5'])
        summary[grid]=dict(points=len(rows),punzi_a5={**best,'at_grid_boundary':best is rows[0] or best is rows[-1],'above_95percent_peak':[r['xeff'] if grid=='xeff' else r['score_threshold'] for r in rows if r['punzi_a5']>=.95*best['punzi_a5']],'above_99percent_peak':[r['xeff'] if grid=='xeff' else r['score_threshold'] for r in rows if r['punzi_a5']>=.99*best['punzi_a5']]})
        allrows+=rows
    assert len(allrows)==65
    dump(out/'summary.json',summary);plot(out)
    pdfs=[out/f'{g}_scan_vs_{c}.pdf' for g in ['xeff','score'] for c in ['xeff','score']]
    for p in pdfs: assert re.search(r'^Pages:\s+1\s*$', subprocess.check_output(['pdfinfo',str(p)],text=True), re.MULTILINE)
    subprocess.run(['pdfunite',*map(str,pdfs),str(out/'punzi_all.pdf')],check=True)
    for record in inputs.values(): assert fingerprint(record['path'])==record
    dump(out/'audit.json',dict(status='PASS',points=65,registered_thresholds_exact=len(registered),apply_selection_and_entries_verified=True,independent_sorted_mass_counts=True,independent_fsum_and_width_ratio=True,numerical_derivative_error_verified=True,raw_inputs_unchanged=True,no_additional_weight_veto=True,new_fits=0))
    best=summary['xeff']['punzi_a5'];coarse=summary['score']['punzi_a5']
    (out/'README.md').write_text(f"# {args.tag}: count Punzi a=5\n\n[四页图册](punzi_all.pdf)\n\n- [细扫描 vs xeff](xeff_scan_vs_xeff.pdf) · [细扫描 vs score](xeff_scan_vs_score.pdf)\n- [粗扫描 vs xeff](score_scan_vs_xeff.pdf) · [粗扫描 vs score](score_scan_vs_score.pdf)\n\nxeff15%–60%每1个百分点（46点）；score0.05–0.95每0.05（19点）。\n\n细扫描最大：目标xeff={best['xeff']}%，实际{best['xeff_actual_percent']:.6f}%，score>{best['score_threshold']:.12g}，Punzi={best['punzi_a5']:.8g}。\n99%峰高点集：{best['above_99percent_peak']}；95%峰高点集：{best['above_95percent_peak']}。这些集合不是置信区间。\n\n粗网格最大：score>{coarse['score_threshold']}，实际xeff={coarse['xeff_actual_percent']:.6f}%。\n\n沿用计数定义：MC在3.62<m<4.00内的加权通过比例，除以2.5+sqrt((N左+N右)/3)；侧带(3.75,3.84)、(3.90,3.99)，MC窗口因子为1。阈值按REFERENCE MC无额外质量窗的pThat加权分位数，严格Prediction>cut。实际selection及权重类别见[contract](contract.json)。不增加旧v31权重剔除。\n\n误差条仅为侧带计数统计误差；各cut高度相关，未计入MC统计或背景模型系统误差。本结果是扫描范围内的优化指标，不是峰显著性或最终选点批准。未拟合DATA峰、未bootstrap。\n\n[细网格CSV](xeff_scan.csv) · [粗网格CSV](score_scan.csv) · [汇总](summary.json) · [验收](audit.json)。输入hash和选择在contract.json；冻结源码在artifacts/。\n")
    dump(out/'delivery_hashes.json',{str(p.relative_to(out)):fingerprint(p)['sha256'] for p in out.rglob('*') if p.is_file()})
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
