#!/usr/bin/env python3
"""31-point 10v1 Punzi refinement with independently reconstructed MC quantiles.

Uses the frozen coarse-study helpers without modifying that workflow or outputs.
Run in LCG_106 ROOT 6.32.02. No upstream writes or ML/model changes.
"""
from pathlib import Path
import argparse
import csv
import json
import hashlib
import numpy as np
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scan_10v1_merged_punzi as base

ROOT.gROOT.SetBatch(True)
GRID=list(range(10,41))
COARSE=base.REPO/'fitER/results/archive/earlier_workflows/10v1_merged_punzi/20260910_a3_equalmix_v1'


def quantiles(scores,weights):
    scores=np.asarray(scores,dtype=float);weights=np.asarray(weights,dtype=float)
    if len(scores)==0 or not np.isfinite(scores).all() or not np.isfinite(weights).all() or np.any(weights<=0):
        raise RuntimeError('reference MC invalid scores or weights')
    order=np.argsort(scores,kind='stable');cumulative=np.cumsum(weights[order])/weights.sum()
    output={}
    for e in GRID:
        i=min(int(np.searchsorted(cumulative,1-e/100.,side='left')),len(scores)-1)
        t=float(scores[order[i]]);eff=float(weights[scores>t].sum()/weights.sum())
        output[e]=dict(target_efficiency=e/100.,score_threshold=t,achieved_efficiency=eff)
        if abs(eff-e/100.)>1e-3:raise RuntimeError('achieved efficiency differs from target by >0.001')
    return output


def inputs():
    manifest,arrays,old_thresholds,pf=base.preflight();thresholds={};references=[];closure=[]
    algorithm=base.REPO.parent/'XGBoost/utils/score_thresholds.py'
    for y in base.YEARS:
        cat=manifest['pairing']['categories'][y]
        path=base.MANIFEST.parent.parent/cat['train_tag']/'REFERENCE_MC_with_score.root'
        before=base.fingerprint(path)
        df=ROOT.RDataFrame(cat['signal_mc']['tree'],str(path))
        cols=[str(c) for c in df.GetColumnNames()]
        score='Prediction_'+cat['train_tag'] if 'Prediction_'+cat['train_tag'] in cols else 'Prediction'
        arr=df.Filter(cat['fiducial_selection']['expression']).AsNumpy([score,'Reweight'])
        after=base.fingerprint(path)
        if before!=after:raise RuntimeError('reference changed during read')
        th=quantiles(arr[score],arr['Reweight'])
        for e in range(10,41,5):
            old=old_thresholds[y][e]
            if th[e]['score_threshold']!=old['score_threshold'] or abs(th[e]['achieved_efficiency']-old['achieved_efficiency'])>1e-12:
                raise RuntimeError(f'{y} xeff{e}: original threshold/efficiency mismatch')
            closure.append(dict(year=y,point=e,threshold_exact=True,achieved_efficiency_difference=th[e]['achieved_efficiency']-old['achieved_efficiency']))
        if min(x['score_threshold'] for x in th.values())<old_thresholds[y][40]['score_threshold']:raise RuntimeError('threshold outside cached support')
        thresholds[y]=th;references.append(dict(**before,year=y,tree=cat['signal_mc']['tree'],selection=cat['fiducial_selection']['expression'],score_branch=score,selected_entries=len(arr[score]),selected_sumw=float(np.asarray(arr['Reweight'],float).sum())))
    pf['reference_mc']=references;pf['original_threshold_closure']=closure
    pf['additional_protected']=[base.fingerprint(Path(__file__)),base.fingerprint(algorithm),base.fingerprint(COARSE/'scan.json')]
    return arrays,thresholds,pf


def scan(out,arrays,thresholds,pf,nboot):
    out.mkdir(parents=True,exist_ok=False)
    base.dump(out/'preflight.json',pf);base.dump(out/'thresholds.json',thresholds)
    contract=dict(grid_percent=GRID,a=3,window=[3.84,3.90],mass_range=[3.8,3.94],mc_mixture=[.5,.5],
                  threshold_definition='full reference MC after original fiducial cut; strict weighted quantile, no mass cut',
                  numerator='target weighted MC efficiency, same as coarse revision',bootstrap=nboot,seed=10012026,
                  plateau_fraction=.95,selection='largest valid sampled Punzi FOM; no observed signal fit used',
                  significance='uncalibrated local fit-only sqrt(q0), no toys/trials/LEE')
    base.dump(out/'contract.json',contract)
    model=base.SidebandFit([3.8,3.94],[3.84,3.90]);sb={}
    for y in base.YEARS:
        d=arrays[y]['data'];mask=(d['Bmass']<3.84)|(d['Bmass']>3.90)
        sb[y]={k:v[mask] for k,v in d.items()}
    rows=[];samples=[];masks={}
    for e in GRID:
        masks[e]={y:sb[y]['Prediction']>thresholds[y][e]['score_threshold'] for y in base.YEARS}
        masses=np.concatenate([sb[y]['Bmass'][masks[e][y]] for y in base.YEARS]);samples.append(masses)
        fit=model.fit(masses);fom=e/100./(1.5+np.sqrt(fit['background']))
        row=dict(point=e,efficiency=e/100.,**fit,fom=float(fom),
                 fom_error_background_only=float(fom*fit['background_error']/(2*np.sqrt(fit['background'])*(1.5+np.sqrt(fit['background'])))),
                 thresholds={y:thresholds[y][e]['score_threshold'] for y in base.YEARS},
                 achieved_efficiencies={y:thresholds[y][e]['achieved_efficiency'] for y in base.YEARS})
        rows.append(row)
    # Reproduce coarse sideband estimates at all seven common points before interpreting refinement.
    old={r['point']:r for r in json.loads((COARSE/'scan.json').read_text())['rows']}
    for r in rows:
        if r['point'] in old:
            for key in ['background','fom']:
                if not np.isclose(r[key],old[r['point']][key],rtol=1e-10,atol=1e-10):raise RuntimeError('coarse FOM closure failed')
    valid=[r for r in rows if r['status']=='PASS'];peak=max(valid,key=lambda r:r['fom']);idx=GRID.index(peak['point'])
    left=right=idx
    while left>0 and rows[left-1]['status']=='PASS' and rows[left-1]['fom']>=.95*peak['fom']:left-=1
    while right+1<len(rows) and rows[right+1]['status']=='PASS' and rows[right+1]['fom']>=.95*peak['fom']:right+=1
    rng=np.random.default_rng(contract['seed']);boot=[]
    for rep in range(nboot):
        counts={y:rng.poisson(1,len(sb[y]['Bmass'])) for y in base.YEARS};values=[]
        for e,mass in zip(GRID,samples):
            w=np.concatenate([counts[y][masks[e][y]] for y in base.YEARS])
            try:
                b=model.fit(mass,w);v=e/100./(1.5+np.sqrt(b['background'])) if b['status']=='PASS' else None
            except RuntimeError:v=None
            values.append(v)
        boot.append(values)
        if (rep+1)%25==0:print('BOOTSTRAP',rep+1,'/',nboot,flush=True)
    complete=[v for v in boot if all(x is not None for x in v)]
    max_points=[GRID[int(np.argmax(v))] for v in complete]
    comparisons=[]
    for e in sorted(set([10,15,25,40,max(10,peak['point']-1),min(40,peak['point']+1)])):
        if e==peak['point']:continue
        j=GRID.index(e);pairs=np.array([[v[idx],v[j]] for v in boot if v[idx] is not None and v[j] is not None]);delta=pairs[:,0]-pairs[:,1]
        comparisons.append(dict(other_point=e,valid_replicates=len(pairs),fraction_peak_greater=float(np.mean(delta>0)),delta_quantiles=np.quantile(delta,[.025,.16,.5,.84,.975]).tolist()))
    summary=dict(status='PASS',sampled_peak=peak['point'],plateau=GRID[left:right+1],boundary_limited=left==0 or right==len(GRID)-1,
                 rows=rows,bootstrap_complete=len(complete),bootstrap_requested=nboot,
                 bootstrap_invalid_by_point={str(e):sum(v[i] is None for v in boot) for i,e in enumerate(GRID)},
                 bootstrap_peak_counts={str(e):max_points.count(e) for e in GRID},bootstrap_peak_quantiles=np.quantile(max_points,[.025,.16,.5,.84,.975]).tolist(),
                 paired_comparisons=comparisons,
                 interpretation='bootstrap fractions and peak quantiles are descriptive, conditional on C2 model and valid fits; not p-values or calibrated optimal-point intervals')
    base.dump(out/'scan.json',summary);base.dump(out/'bootstrap.json',dict(grid=GRID,replicates=boot))
    base.dump(out/'selected_point.json',dict(point=peak['point'],thresholds=peak['thresholds'],scan_sha256=base.digest(out/'scan.json'),contract_sha256=base.digest(out/'contract.json'),selection='maximum Punzi FOM fixed before this point signal-window fit'))
    with (out/'scan.csv').open('w') as f:
        fields=['point','efficiency','background','background_error','fom','fom_error_background_only','status'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in rows)
    fig,axs=plt.subplots(1,2,figsize=(12,4.8));ax=axs[0]
    ax.errorbar(GRID,[r['fom']*1000 for r in rows],yerr=[r['fom_error_background_only']*1000 for r in rows],fmt='o-',ms=3,capsize=2,label='FOM ± background statistical error')
    ax.axhline(.95*peak['fom']*1000,color='gray',ls='--',label='95% of sampled peak');ax.axvspan(GRID[left],GRID[right],color='orange',alpha=.15)
    ax.axvline(peak['point'],color='red',ls=':',label=f"Peak: {peak['point']}%")
    ax.set(xlabel='Common weighted MC efficiency [%]',ylabel='Punzi FOM × 1000',title='10v1 merged Punzi, a=3');ax.legend(fontsize=8)
    axs[1].bar(GRID,[max_points.count(e) for e in GRID],width=.8)
    axs[1].set(xlabel='Efficiency at sampled maximum [%]',ylabel='Bootstrap count',title=f'Maxima: {len(complete)}/{nboot} complete replicas')
    fig.tight_layout();fig.savefig(out/'fine_punzi.pdf');fig.savefig(out/'fine_punzi.png',dpi=160);plt.close(fig)
    print('SCAN_PEAK',peak['point'],'PLATEAU',GRID[left:right+1],'BOOT_COMPLETE',len(complete),flush=True)


def fit(out,arrays,thresholds,pf):
    selected=json.loads((out/'selected_point.json').read_text())
    if base.digest(out/'scan.json')!=selected['scan_sha256'] or base.digest(out/'contract.json')!=selected['contract_sha256']:raise RuntimeError('selection receipt changed')
    if pf!=json.loads((out/'preflight.json').read_text()):raise RuntimeError('inputs changed since scan')
    fitroot=out/'merged_fits';fitroot.mkdir(exist_ok=False)
    result=base.merged_fit(fitroot,selected['point'],arrays,thresholds,[.5,.5])
    base.dump(out/'maximum_fom_fit.json',result)
    protected=pf['inputs']+pf['additional_protected']+[{k:r[k] for k in ['path','size','mtime_ns','sha256']} for r in pf['reference_mc']]
    assert [base.fingerprint(r['path']) for r in protected]==protected,'protected inputs changed'
    base.dump(out/'validation.json',dict(status='PASS' if result['quality_pass'] else 'COMPLETED_WITH_FIT_FAILURE',threshold_closure_14_of_14=True,
        common_coarse_fom_closure_7_of_7=True,protected_inputs_unchanged=True,selected_point=selected['point'],quality_pass=result['quality_pass'],
        final_working_point_selected=False,calibrated_significance=False))
    print('FIT_COMPLETE',result['point'],result['signal_yield'],result['signal_yield_error'],result['local_significance'],result['quality_pass'],flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['preflight','scan','fit','all']);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--bootstrap',type=int,default=200);args=ap.parse_args()
    arrays,thresholds,pf=inputs();out=args.output.resolve()
    if args.stage=='preflight':print(json.dumps(pf,indent=2));return
    if args.stage in ['scan','all']:scan(out,arrays,thresholds,pf,args.bootstrap)
    if args.stage in ['fit','all']:fit(out,arrays,thresholds,pf)


if __name__=='__main__':main()
