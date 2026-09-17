#!/usr/bin/env python3
"""Isolated 10v1 common-efficiency, merged-sideband Punzi and merged-fit study.

Run under LCG_106 ROOT 6.32.02. Physics parameters are explicit CLI inputs.
The scan consumes frozen historical compact caches; it does not rebuild upstream
ML artifacts or claim that original source files have been revalidated.
"""
from __future__ import annotations
import argparse
from array import array
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
TAG = "X_pb23_v3_fid3_10v1_rwr6range5v1_xgb_v1"
MANIFEST = REPO.parent / "XGBoost/output/selected" / TAG / "fit_scan_manifest.pb23_pb24_simultaneous_mc_shape_nominal_v2.json"
CACHE_ROOT = REPO / "fitER/results/archive/earlier_workflows/pbpb_x_simultaneous_year_fit" / TAG / "mc_shape_nominal_v2_fit_only_sqrtq0"
YEARS = ("pb23", "pb24")
EXPECTED_HASH = "2e78d31d598d7ade93c8094a243891f1b6c4540592c5e3e400e39b3f89b5f0c0"
GRID = [5,10,15,20,25,30,35,40]


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(p, d):
    Path(p).write_text(json.dumps(d, indent=2, allow_nan=False) + "\n")


def fingerprint(p):
    p=Path(p); before=p.stat(); h=digest(p); after=p.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
        raise RuntimeError(f"input changed while hashing: {p}")
    return dict(path=str(p),size=after.st_size,mtime_ns=after.st_mtime_ns,sha256=h)


def basis_integrals(lo,hi):
    return np.array([hi-lo,(hi*hi-lo*lo)/2,2*(hi**3-lo**3)/3-(hi-lo)])


def polynomial_min(c):
    xs=[-1.,1.]
    if c[1]>0:
        v=-c[0]/(4*c[1])
        if -1<v<1: xs.append(v)
    return min(1+c[0]*x+c[1]*(2*x*x-1) for x in xs)


class SidebandFit:
    """Conditional unbinned C2 fit; sideband count is a separate Poisson rate."""
    def __init__(self, mass_range, window):
        self.lo,self.hi=mass_range; self.sl,self.sh=window
        self.x=lambda m: 2*(np.asarray(m)-self.lo)/(self.hi-self.lo)-1
        self.sig=basis_integrals(float(self.x(self.sl)),float(self.x(self.sh)))
        self.sb=basis_integrals(-1.,float(self.x(self.sl)))+basis_integrals(float(self.x(self.sh)),1.)

    def fit(self, masses, weights=None):
        masses=np.asarray(masses)
        assert np.all((masses<self.sl)|(masses>self.sh)), "signal-window data in sideband fit"
        x=self.x(masses); u=np.column_stack([x,2*x*x-1])
        w=np.ones(len(x)) if weights is None else np.asarray(weights)
        n=float(w.sum())
        if n<10: raise RuntimeError("fewer than 10 effective sideband events")
        def fun(c):
            f=1+u@c; integ=self.sb[0]+self.sb[1:]@c
            if integ<=0 or np.any(f<=0): return 1e90,np.zeros(2)
            return n*np.log(integ)-np.dot(w,np.log(f)), n*self.sb[1:]/integ-((w/f)[:,None]*u).sum(axis=0)
        fits=[minimize(fun, start, jac=True, method="SLSQP", bounds=[(-2,2),(-2,2)],
                       constraints=[dict(type="ineq",fun=lambda c: polynomial_min(c)-1e-8)],
                       options=dict(ftol=1e-10,maxiter=300))
              for start in ([0.,0.],[-.35,-.05])]
        fits=[r for r in fits if r.success and polynomial_min(r.x)>0]
        if not fits: raise RuntimeError("positive C2 sideband fit failed")
        r=min(fits,key=lambda z:z.fun);c=r.x
        f=1+u@c; isb=self.sb[0]+self.sb[1:]@c; isig=self.sig[0]+self.sig[1:]@c
        hess=u.T@((w/f**2)[:,None]*u)-n*np.outer(self.sb[1:],self.sb[1:])/isb**2
        if np.linalg.eigvalsh(hess).min()<=0: raise RuntimeError("nonpositive sideband Hessian")
        cov=np.linalg.inv(hess); ratio=isig/isb; grad=(self.sig[1:]*isb-isig*self.sb[1:])/isb**2
        b=n*ratio; var=n*ratio**2+n*n*float(grad@cov@grad)
        boundary=bool(np.max(np.abs(c))>1.999 or polynomial_min(c)<1e-5)
        edges=np.linspace(self.lo,self.hi,29)
        hist,_=np.histogram(masses,bins=edges,weights=w)
        valid=(edges[1:]<=self.sl+1e-10)|(edges[:-1]>=self.sh-1e-10)
        expected=np.array([n*(basis_integrals(float(self.x(l)),float(self.x(h)))@np.r_[1,c])/isb for l,h in zip(edges[:-1],edges[1:])])
        pearson=float(np.sum((hist[valid]-expected[valid])**2/expected[valid]))
        ndf=int(valid.sum()-3)
        return dict(background=b,background_error=math.sqrt(var),sideband_entries=n,
                    coefficients=c.tolist(),covariance=cov.tolist(),nll=float(r.fun),
                    positive_min=polynomial_min(c),parameter_boundary=boundary,
                    pearson_chi2=pearson,pearson_ndf=ndf,pearson_pvalue_asymptotic=float(chi2.sf(pearson,ndf)),
                    minimum_expected_bin=float(expected[valid].min()),status="PASS" if not boundary else "BOUNDARY")


def read_arrays(path,tree,branches):
    import ROOT
    f=ROOT.TFile.Open(str(path)); t=f.Get(tree)
    if not t: raise RuntimeError(f"missing {path}:{tree}")
    n=t.GetEntries();f.Close()
    result=ROOT.RDataFrame(tree,str(path)).AsNumpy(branches)
    for name,a in result.items():
        if not np.all(np.isfinite(a)): raise RuntimeError(f"nonfinite {path}:{name}")
    return result,n


def preflight():
    import ROOT
    if ROOT.gROOT.GetVersion()!="6.32.02": raise RuntimeError("ROOT 6.32.02 required")
    if digest(MANIFEST)!=EXPECTED_HASH: raise RuntimeError("10v1 manifest changed")
    m=json.loads(MANIFEST.read_text());context_path=CACHE_ROOT/"run_context.json"
    ctx=json.loads(context_path.read_text())
    if ctx["input_manifest_sha256"]!=EXPECTED_HASH: raise RuntimeError("cache manifest mismatch")
    if ctx["cache_mass_range_gev"]!=[3.8,3.94]: raise RuntimeError("unexpected cache mass range")
    protected=[MANIFEST,context_path,REPO/"fitER/models/PbPbXEfficiencyFit.C",Path(__file__)]
    arrays={};thresholds={};metas={}
    for y in YEARS:
        cat=m["pairing"]["categories"][y];meta=ctx["cache_metadata"][y]
        tp=(MANIFEST.parent/cat["threshold_provenance"]["path"]).resolve()
        if digest(tp)!=cat["threshold_provenance"]["sha256"]: raise RuntimeError("threshold hash mismatch")
        protected.append(tp); th=json.loads(tp.read_text())["thresholds"]
        thresholds[y]={int(round(r["target_efficiency"]*100)):r for r in th}
        for pt in m["working_points"]:
            e=int(round(pt["target_weighted_efficiency"]*100))
            if pt["categories"][y]["threshold"]!=thresholds[y][e]["score_threshold"]: raise RuntimeError("threshold closure failed")
        arrays[y]={}
        for k,ct in [("data","data"),("mc","signal_mc")]:
            p=Path(meta[f"{k}_cache"]);protected.append(p)
            branches=["Bmass","Prediction","source_entry"]+(["Reweight"] if k=="mc" else [])
            arr,n=read_arrays(p,cat[ct]["tree"],branches)
            if n!=meta[f"{k}_entries"]: raise RuntimeError("cache entry mismatch")
            if len(np.unique(arr["source_entry"]))!=n: raise RuntimeError("duplicate cache source entries")
            if not np.all((arr["Bmass"]>3.8)&(arr["Bmass"]<3.94)): raise RuntimeError("mass cache mismatch")
            if not np.all(arr["Prediction"]>thresholds[y][40]["score_threshold"]): raise RuntimeError("score cache mismatch")
            if k=="mc" and not np.all(arr["Reweight"]>0): raise RuntimeError("nonpositive MC weight")
            arrays[y][k]=arr
        metas[y]=dict(data_entries=len(arrays[y]["data"]["Bmass"]),mc_entries=len(arrays[y]["mc"]["Bmass"]))
    fps=[fingerprint(p) for p in protected]
    return m,arrays,thresholds,dict(status="PASS",root_version=ROOT.gROOT.GetVersion(),inputs=fps,cache_counts=metas,
        provenance="Frozen historical 10v1 caches; original large scored DATA files are not reread. Manifest and threshold hashes verified; cache hashes frozen for this revision.")


def write_tree(path, tree_name, arrays):
    import ROOT
    f=ROOT.TFile(str(path),"RECREATE");t=ROOT.TTree(tree_name,tree_name)
    slots={k:array('d',[0.]) for k in arrays}
    for k,b in slots.items():t.Branch(k,b,k+"/D")
    for i in range(len(next(iter(arrays.values())))):
        for k,b in slots.items():b[0]=float(arrays[k][i])
        t.Fill()
    t.Write();f.Close()


def merged_fit(out, percent, arrays, thresholds, mixture):
    work=out/f"xeff{percent:02d}";work.mkdir()
    ds=[];ms=[];ws=[];ys=[];meta={}
    for yi,y in enumerate(YEARS):
        thr=thresholds[y][percent]["score_threshold"]
        d=arrays[y]["data"];mc=arrays[y]["mc"]
        ds.append(d["Bmass"][d["Prediction"]>thr])
        mask=mc["Prediction"]>thr;mass=mc["Bmass"][mask];weights=mc["Reweight"][mask]
        # Mixture coefficients refer to integrals in the common full fit range.
        target=mixture[yi]*10000.
        ww=weights*(target/weights.sum());ms.append(mass);ws.append(ww);ys.append(np.full(len(mass),yi))
        meta[y]=dict(data_entries=len(ds[-1]),mc_entries=len(mass),source_sumw=float(weights.sum()),
                     scaled_sumw=float(ww.sum()),target_sumw=target,mc_neff=float(weights.sum()**2/np.dot(weights,weights)))
        if not np.isclose(ww.sum(),target,rtol=1e-12):raise RuntimeError("MC mixture closure failed")
    data_path=work/"merged_data.root";mc_path=work/"merged_mc.root"
    write_tree(data_path,"ntmix",dict(Bmass=np.concatenate(ds)))
    write_tree(mc_path,"ntmix_X3872",dict(Bmass=np.concatenate(ms),MixtureWeight=np.concatenate(ws),year_index=np.concatenate(ys)))
    dump(work/"mixture_metadata.json",dict(years=meta,mixture=mixture,normalization_range=[3.8,3.94],arbitrary_common_weight_scale=10000))
    scratch=work/"build";scratch.mkdir();macro=REPO/"fitER/models/PbPbXEfficiencyFit.C";shutil.copy2(macro,scratch/macro.name)
    vals=[f"10v1_xeff{percent}_fixed_mixture",data_path,"ntmix",mc_path,"ntmix_X3872","1","1","MixtureWeight"]
    expr="PbPbXEfficiencyFit.C++("+",".join(json.dumps(str(x)) for x in vals)+",3.8,3.94,3.87169,0.005,0.9,1.5,28,"+json.dumps(str(work))+")"
    with (work/"fit.log").open('w') as log:
        rc=subprocess.run([shutil.which("root"),"-l","-b","-q",expr],cwd=scratch,stdout=log,stderr=subprocess.STDOUT).returncode
    if rc or not (work/"fit_result.json").exists():raise RuntimeError(f"merged fit failed: {work}")
    r=json.loads((work/"fit_result.json").read_text());r["point"]=percent
    r["quality_pass"]=bool(r["fit_status"]==0 and r["cov_qual"]==3 and r["edm"]<1e-3 and r["signal_mc_fit_status"]==0 and r["signal_mc_cov_qual"]==3 and r["signal_mc_edm"]<1e-3)
    # Independently check the null fit because legacy JSON does not expose its quality.
    import ROOT
    files=list(work.glob('*.root'))
    null_checks=[]
    for fpath in files:
        f=ROOT.TFile.Open(str(fpath))
        for key in f.GetListOfKeys():
            obj=f.Get(key.GetName())
            if obj.InheritsFrom('RooFitResult') and 'null' in key.GetName():
                null_checks.append(dict(status=obj.status(),covQual=obj.covQual(),edm=obj.edm()))
        f.Close()
    r["null_fit_quality"]=null_checks
    r["quality_pass"]=r["quality_pass"] and bool(null_checks) and all(x['status']==0 and x['covQual']==3 and x['edm']<1e-3 for x in null_checks)
    return r


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['preflight','scan','fit'])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--a',type=float)
    parser.add_argument('--window',type=float,nargs=2)
    parser.add_argument('--mixture',type=float,nargs=2)
    parser.add_argument('--bootstrap',type=int,default=200)
    args=parser.parse_args()
    m,arrays,thresholds,pf=preflight();out=args.output.resolve()
    if args.stage=='preflight':
        print(json.dumps(pf,indent=2));return
    if args.stage=='scan':
        if args.a is None or args.window is None or args.mixture is None:parser.error('explicit --a --window --mixture required')
        if not (args.a>0 and 3.8<args.window[0]<args.window[1]<3.94 and min(args.mixture)>0):parser.error('invalid scientific parameters')
        if not math.isclose(sum(args.mixture),1.):parser.error('mixture must sum to 1')
        out.mkdir(parents=True,exist_ok=False)
        contract=dict(created_utc=datetime.now(timezone.utc).isoformat(),a=args.a,window=args.window,mass_range=[3.8,3.94],
            mixture=args.mixture,mixture_range=[3.8,3.94],efficiency_grid_percent=GRID,bootstrap=args.bootstrap,seed=10012026,
            plateau_fraction=.95,plateau_definition='connected sampled points with FOM >= 95% of sampled peak; not a confidence interval',
            candidate_rule='plateau endpoints, sampled peak, and historical xeff25, frozen before signal fits',
            efficiency_numerator='target full-selection weighted MC efficiency, approximate; finite achieved efficiencies and mass-window acceptance recorded separately',
            observed_z_used_for_selection=False,background_model='positive Chebyshev order 2; conditional sideband unbinned likelihood + Poisson count',
            signal_fit='existing PbPbXEfficiencyFit.C; full-range fixed signal mixture, shared single-category double Gaussian, C2 background')
        dump(out/'contract.json',contract);dump(out/'preflight.json',pf)
        model=SidebandFit([3.8,3.94],args.window);rows=[];samples=[];rng=np.random.default_rng(contract['seed'])
        sb_all={}
        for y in YEARS:
            d=arrays[y]['data'];mask=(d['Bmass']<args.window[0])|(d['Bmass']>args.window[1])
            sb_all[y]={k:v[mask] for k,v in d.items()}
        for e in GRID:
            parts=[];record=dict(point=e,efficiency=e/100.,thresholds={},achieved_efficiencies={},mc_window_acceptance={},mc_neff={},year_data_entries={})
            for y in YEARS:
                th=thresholds[y][e];t=th['score_threshold'];record['thresholds'][y]=t;record['achieved_efficiencies'][y]=th['achieved_efficiency']
                if abs(th['achieved_efficiency']-e/100.)>1e-3:raise RuntimeError('efficiency mismatch')
                d=sb_all[y];parts.append(d['Bmass'][d['Prediction']>t])
                record['year_data_entries'][y]=int(np.count_nonzero(arrays[y]['data']['Prediction']>t))
                mc=arrays[y]['mc'];mask=mc['Prediction']>t;w=mc['Reweight'][mask];mass=mc['Bmass'][mask]
                record['mc_window_acceptance'][y]=float(w[(mass>args.window[0])&(mass<args.window[1])].sum()/w.sum())
                record['mc_neff'][y]=float(w.sum()**2/np.dot(w,w))
            masses=np.concatenate(parts);fit=model.fit(masses);record.update(fit)
            record['fom']=record['efficiency']/(args.a/2+math.sqrt(fit['background']))
            record['fom_error_background_only']=record['fom']*fit['background_error']/(2*math.sqrt(fit['background'])*(args.a/2+math.sqrt(fit['background'])))
            rows.append(record);samples.append(masses)
        usable=[r for r in rows if r['status']=='PASS']
        if not usable:raise RuntimeError('no interior sideband fits')
        best=max(usable,key=lambda r:r['fom']);peakidx=GRID.index(best['point']);floor=.95*best['fom']
        left=right=peakidx
        while left>0 and rows[left-1]['status']=='PASS' and rows[left-1]['fom']>=floor:left-=1
        while right+1<len(rows) and rows[right+1]['status']=='PASS' and rows[right+1]['fom']>=floor:right+=1
        plateau=GRID[left:right+1];candidates=sorted(set([plateau[0],plateau[-1],best['point'],25]))
        # Shared event-level Poisson draws retain correlations between nested cuts.
        boot=[]
        for rep in range(args.bootstrap):
            counts={y:rng.poisson(1,len(sb_all[y]['Bmass'])) for y in YEARS};values=[]
            for e in GRID:
                weights=np.concatenate([counts[y][sb_all[y]['Prediction']>thresholds[y][e]['score_threshold']] for y in YEARS])
                try:
                    fit=model.fit(samples[GRID.index(e)],weights)
                    values.append(e/100./(args.a/2+math.sqrt(fit['background'])) if fit['status']=='PASS' else None)
                except RuntimeError:values.append(None)
            boot.append(values)
        maxima={str(e):0 for e in GRID};valid=0
        for v in boot:
            if all(x is not None for x in v):maxima[str(GRID[int(np.argmax(v))])]+=1;valid+=1
        summary=dict(status='PASS',sampled_peak=best['point'],plateau=plateau,candidates=candidates,
            boundary_limited=left==0 or right==len(GRID)-1,bootstrap_complete_replicates=valid,
            bootstrap_requested=args.bootstrap,bootstrap_peak_counts=maxima,rows=rows,
            limitations=['Only frozen 5-40% grid evaluated; no global optimum claim.',
                        'Punzi uses full weighted MC efficiency; mass-window acceptance reported, not unfolded.',
                        'Sideband-only C2 extrapolation assumes negligible signal leakage; not independently validated.',
                        'Bootstrap covers background statistical variation, not transfer, MC efficiency or model systematics.',
                        '95% FOM plateau is a descriptive candidate region, not a confidence interval.',
                        'Historical signal-window scans already exist; not a blind analysis.'])
        dump(out/'scan.json',summary);dump(out/'bootstrap.json',dict(grid=GRID,replicates=boot))
        dump(out/'candidate_selection.json',dict(candidates=candidates,scan_sha256=digest(out/'scan.json'),rule=contract['candidate_rule'],contract_sha256=digest(out/'contract.json')))
        with (out/'scan.csv').open('w') as f:
            keys=['point','efficiency','background','background_error','sideband_entries','fom','fom_error_background_only','status']
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows({k:r[k] for k in keys} for r in rows)
        fig,ax=plt.subplots(figsize=(7,5));ax.errorbar(GRID,[r['fom'] for r in rows],yerr=[r['fom_error_background_only'] for r in rows],fmt='o-',capsize=3)
        ax.axhline(floor,color='gray',ls='--',label='95% of sampled peak')
        ax.set(xlabel='Common weighted MC efficiency [%]',ylabel='Merged Punzi FOM',title=f'10v1, a={args.a:g}; background statistical errors only');ax.legend();fig.tight_layout();fig.savefig(out/'punzi.pdf');fig.savefig(out/'punzi.png',dpi=150);plt.close(fig)
        fig,axs=plt.subplots(2,4,figsize=(14,7))
        for ax,r,masses in zip(axs.flat,rows,samples):
            edges=np.linspace(3.8,3.94,29);hist,_=np.histogram(masses,bins=edges);centers=(edges[:-1]+edges[1:])/2
            mask=(centers<args.window[0])|(centers>args.window[1]);ax.errorbar(centers[mask],hist[mask],np.sqrt(hist[mask]),fmt='.',color='black')
            xx=np.linspace(3.8,3.94,400);x=model.x(xx);c=np.array(r['coefficients']);isb=model.sb[0]+model.sb[1:]@c
            density=r['sideband_entries']*(1+c[0]*x+c[1]*(2*x*x-1))/isb*2/(3.94-3.8)*.005
            ax.plot(xx,density);ax.axvspan(*args.window,alpha=.15,color='gray');ax.set_title(f"xeff{r['point']:02d}, B={r['background']:.1f}");ax.set_xlabel('Bmass [GeV]')
        fig.tight_layout();fig.savefig(out/'sideband_fits.pdf');fig.savefig(out/'sideband_fits.png',dpi=150);plt.close(fig)
        print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))
    else:
        contract=json.loads((out/'contract.json').read_text());selection=json.loads((out/'candidate_selection.json').read_text())
        if digest(out/'scan.json')!=selection['scan_sha256'] or digest(out/'contract.json')!=selection['contract_sha256']:raise RuntimeError('candidate selection changed')
        oldpf=json.loads((out/'preflight.json').read_text())
        if pf['inputs']!=oldpf['inputs']:raise RuntimeError('inputs changed since scan')
        fitroot=out/'merged_fits';fitroot.mkdir(exist_ok=False)
        results=[]
        for e in selection['candidates']:
            r=merged_fit(fitroot,e,arrays,thresholds,contract['mixture']);results.append(r)
            dump(out/'merged_fit_summary.json',dict(results=results,complete=False))
            print('FIT_COMPLETE',e,r['signal_yield'],r['local_significance'],r['quality_pass'],flush=True)
        dump(out/'merged_fit_summary.json',dict(results=results,complete=True,selection='pre-frozen Punzi candidates; no observed-Z reselection'))
        current=[fingerprint(f['path']) for f in pf['inputs']]
        if current!=pf['inputs']:raise RuntimeError('protected input changed')
        dump(out/'validation.json',dict(status='PASS' if all(r['quality_pass'] for r in results) else 'COMPLETED_WITH_FIT_FAILURES',
            protected_inputs_unchanged=True,candidates_complete=True,fit_quality=[dict(point=r['point'],passed=r['quality_pass']) for r in results],
            physics_interpretation='descriptive study; no calibrated significance or final WP chosen'))


if __name__=='__main__':main()
