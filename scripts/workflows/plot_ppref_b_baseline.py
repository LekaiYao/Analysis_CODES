#!/usr/bin/env python3
"""Raw ppRef B DATA mass and scalar distributions; no selection or fit."""
import argparse, hashlib, json, math, re, shutil, subprocess
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import uproot
import ROOT as R

REPO = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True, help='New staging directory')
args = parser.parse_args()
assert not args.output.exists(), args.output
args.output.mkdir(parents=True)
assert R.gROOT.GetVersion() == '6.32.02'
R.gROOT.SetBatch(True)
R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)
R.gStyle.SetTextFont(42)
source_inventory = REPO.parent/'XGBoost/output/internal/code_validation/b_samples_request_20261013/inventory.json'
inventory = json.loads(source_inventory.read_text())
template = REPO/'plotER/validation_archive/aux.h'
configs = {}
for line in template.read_text().splitlines():
    q = re.search(r'\{"([^"]+)",\s*"([^"]+)",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)', line)
    if q:
        try: configs[q[1]] = (float(q[3]), float(q[4]), q[5]=='true', q[2].split(';')[1])
        except ValueError: pass
excluded = {'RunNo','EvtNo','LumiNo','eventIndex','candidateIndex','Bindex'}
labels = {'Bu':'B^{+}', 'Bs':'B_{s}^{0}', 'Bd':'B^{0}'}
mass_titles = {'Bu':'m(J/#psi K^{+}) [GeV/c^{2}]', 'Bs':'m(J/#psi K^{+}K^{-}) [GeV/c^{2}]', 'Bd':'m(J/#psi K#pi) [GeV/c^{2}]'}
roots = {k: args.output/k for k in ['Bmass_ppRef','Variables_ppRef']}
for p in roots.values(): (p/'metadata').mkdir(parents=True)
provenance = {'created_UTC':datetime.now(timezone.utc).isoformat(), 'selection':'1', 'weights':'unit candidate counts', 'scope':'raw ppRef DATA; no added cut, fit, sPlot or MC', 'ROOT_version':R.gROOT.GetVersion(), 'source_inventory':str(source_inventory), 'source_inventory_sha256':hashlib.sha256(source_inventory.read_bytes()).hexdigest(), 'template':str(template), 'template_sha256':hashlib.sha256(template.read_bytes()).hexdigest(), 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'source_fingerprint':'size/mtime_ns before and after; no full source SHA computed', 'channels':{}}
all_config = {}
rootfiles = {k:R.TFile(str(p/'metadata/histograms.root'),'RECREATE') for k,p in roots.items()}

def dump(path, value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def fingerprint(path):
    s=path.stat(); return {'size':s.st_size,'mtime_ns':s.st_mtime_ns}

def plot(channel, name, cfg, count, stats, destination, mass=False):
    lo,hi=cfg['range']; n=cfg['bins']
    h=R.TH1D(channel+'_'+name,'',n,lo,hi); h.SetDirectory(0); h.Sumw2()
    for i,v in enumerate(count,1): h.SetBinContent(i,int(v)); h.SetBinError(i,math.sqrt(int(v)))
    h.SetBinContent(0,stats['underflow']); h.SetBinError(0,math.sqrt(stats['underflow']))
    h.SetBinContent(n+1,stats['overflow']); h.SetBinError(n+1,math.sqrt(stats['overflow']))
    h.SetEntries(stats['finite'])
    c=R.TCanvas('canvas_'+channel+'_'+name,'',760,650)
    c.SetLeftMargin(.14); c.SetRightMargin(.04); c.SetBottomMargin(.13); c.SetTopMargin(.07)
    h.SetLineColor(R.kBlue+1); h.SetMarkerColor(R.kBlue+1); h.SetMarkerStyle(20); h.SetMarkerSize(.65); h.SetLineWidth(2)
    h.SetMinimum(0); h.SetMaximum(max(1,float(max(count)))*1.55)
    h.GetXaxis().SetTitle(cfg['title']); h.GetYaxis().SetTitle('Candidates / (10 MeV/c^{2})' if mass else 'Candidates / bin')
    h.GetYaxis().SetTitleOffset(1.5); h.Draw('E1')
    lab=R.TLatex(); lab.SetNDC(); lab.SetTextFont(42); lab.SetTextSize(.042)
    lab.DrawLatex(.18,.87,'ppRef '+labels[channel]+' DATA')
    lab.SetTextSize(.030); lab.DrawLatex(.18,.815,'Raw input; no additional cuts')
    lab.DrawLatex(.18,.77,f"N = {stats['total']:,}; in range = {int(sum(count)):,}")
    lab.SetTextSize(.025)
    fraction=(stats['underflow']+stats['overflow'])/max(1,stats['finite'])
    lab.DrawLatex(.18,.725,f"Outside range: {100*fraction:.2f}%; non-finite: {stats['nonfinite']:,}")
    if stats['min']==stats['max']: lab.DrawLatex(.18,.68,'Constant branch in this input')
    c.RedrawAxis(); c.Print(str(destination)); c.Close()
    key='Bmass_ppRef' if mass else 'Variables_ppRef'; rootfiles[key].cd(); h.Write()

for channel in ['Bu','Bs','Bd']:
    entry=next(x for x in inventory['inputs'] if x['channel']==channel and x['system']=='ppRef' and x['kind']=='DATA')
    src=Path(entry['path']); before=fingerprint(src)
    assert before=={k:entry[k] for k in ['size','mtime_ns']}
    with uproot.open(src) as f:
        tree=f[entry['tree']]; assert tree.num_entries==entry['entries']
        names=[n for n in tree.keys() if n not in excluded and n!='Bmass']
        names.sort(key=lambda n:(0 if n in ('Bpt','By') else 1 if n.startswith(('Balpha','Bdtheta','Bcos','Bnorm_sv','Bsv','BLxy','Bchi','Bvtx')) else 2 if n.startswith(('Btrk','Bnorm_trk')) else 3 if n.startswith(('Btktk','Bdoublet')) else 4 if n.startswith(('Bmu','Buj')) else 5,n))
        stats={n:{'total':0,'finite':0,'nonfinite':0,'min':None,'max':None} for n in ['Bmass']+names}
        for arrays in tree.iterate(['Bmass']+names, step_size='80 MB', library='np'):
            for n,x in arrays.items():
                assert x.ndim==1 and x.dtype.kind in 'biuf', (n,x.dtype)
                ok=np.isfinite(x); z=x[ok]; s=stats[n]; s['total']+=len(x); s['finite']+=len(z); s['nonfinite']+=int((~ok).sum())
                if len(z):
                    s['min']=float(z.min()) if s['min'] is None else min(s['min'],float(z.min()))
                    s['max']=float(z.max()) if s['max'] is None else max(s['max'],float(z.max()))
        cfgs={}
        for n in names:
            s=stats[n]; lo,hi=s['min'],s['max']; ab=False; title=n; origin='full finite support'
            if lo is None: lo,hi=0.,1.
            if n in configs and n!='Btktkmass':
                lo,hi,ab,title=configs[n]; origin='existing aux.h'
                if channel=='Bu' and n=='Btrk1Pt': hi=10.
                if channel=='Bu' and n=='Btrk1dR': hi=1.5
            elif lo==hi: lo,hi=lo-.5,hi+.5
            else:
                pad=(hi-lo)*.025; lo,hi=lo-pad,hi+pad
            # A legacy display range with no support must not hide the raw B input.
            smin,smax=s['min'],s['max']
            if smin is not None:
                dmin = (0. if smin <= 0 <= smax else min(abs(smin),abs(smax))) if ab else smin
                dmax = max(abs(smin),abs(smax)) if ab else smax
                if dmax < lo or dmin >= hi:
                    pad=.5 if dmin==dmax else (dmax-dmin)*.025
                    lo,hi=dmin-pad,dmax+pad
                    origin='full finite support: legacy range has no entries'
            cfgs[n]={'range':[lo,hi],'bins':15,'absolute_display':ab,'title':title,'range_source':origin,'log_y':False}
        cfgs['Bmass']={'range':[5.,6.],'bins':100,'absolute_display':False,'title':mass_titles[channel],'range_source':'user requested 5 < Bmass < 6 GeV','log_y':False}
        counts={n:np.zeros(c['bins'],dtype=np.int64) for n,c in cfgs.items()}
        for n,c in cfgs.items(): stats[n].update(underflow=0,overflow=0)
        for arrays in tree.iterate(['Bmass']+names,step_size='80 MB',library='np'):
            for n,x in arrays.items():
                c=cfgs[n]; x=x[np.isfinite(x)].astype('float64'); x=np.abs(x) if c['absolute_display'] else x
                lo,hi=c['range']; inside=((x>lo) if n=='Bmass' else (x>=lo))&(x<hi)
                counts[n]+=np.histogram(x[inside],bins=np.linspace(lo,hi,c['bins']+1))[0]
                stats[n]['underflow']+=int(((x<=lo) if n=='Bmass' else (x<lo)).sum()); stats[n]['overflow']+=int((x>=hi).sum())
    for n,s in stats.items():
        assert s['total']==entry['entries']
        assert int(counts[n].sum())+s['underflow']+s['overflow']==s['finite']
        assert s['finite']+s['nonfinite']==s['total']
    # Independent ROOT read and binning for the mass spectrum, including flows.
    df=R.RDataFrame(entry['tree'],str(src)); c=cfgs['Bmass']
    hp=df.Filter('std::isfinite(Bmass) && Bmass > 5 && Bmass < 6').Histo1D((channel+'_audit','',c['bins'],*c['range']),'Bmass'); rh=hp.GetValue()
    assert np.array_equal(counts['Bmass'],[rh.GetBinContent(i) for i in range(1,c['bins']+1)])
    assert df.Filter('Bmass <= 5').Count().GetValue()==stats['Bmass']['underflow'] and df.Filter('Bmass >= 6').Count().GetValue()==stats['Bmass']['overflow']
    for key,p in roots.items(): (p/channel/'pdf').mkdir(parents=True)
    plot(channel,'Bmass',cfgs['Bmass'],counts['Bmass'],stats['Bmass'],roots['Bmass_ppRef']/channel/'pdf/baseline.pdf',True)
    pdfs=[]
    for n in names:
        p=roots['Variables_ppRef']/channel/'pdf'/f'{n}.pdf'; plot(channel,n,cfgs[n],counts[n],stats[n],p); pdfs.append(str(p))
    subprocess.run(['pdfunite',*pdfs,str(roots['Variables_ppRef']/channel/'variables_all.pdf')],check=True)
    # One mass page is already directly accessible; no redundant one-page all_cuts copy.
    assert fingerprint(src)==before
    provenance['channels'][channel]={'source':str(src),'tree':entry['tree'],'fingerprint':before,'entries':entry['entries'],'excluded_identifiers':sorted(excluded),'variables':names,'statistics':stats,'mass_ROOT_numpy_bin_check':'PASS'}
    all_config[channel]=cfgs
    for key in roots:
        dump(roots[key]/'metadata'/f'{channel}_bins.json',{n:{'counts':counts[n].tolist(),'config':cfgs[n]} for n in (['Bmass'] if key=='Bmass_ppRef' else names)})
    print(channel,entry['entries'],len(names),'variables PASS',flush=True)
for f in rootfiles.values(): f.Close()
# Validate persisted bin counts/errors and PDF page counts.
checks={}
for key,p in roots.items():
    with uproot.open(p/'metadata/histograms.root') as f:
        for channel in ['Bu','Bs','Bd']:
            data=json.loads((p/'metadata'/f'{channel}_bins.json').read_text())
            for n,d in data.items():
                h=f[channel+'_'+n]; assert np.array_equal(h.values(),d['counts']); assert np.allclose(h.variances(),d['counts'],rtol=1e-14,atol=1e-12)
    for channel in ['Bu','Bs','Bd']:
        pdf=p/channel/('pdf/baseline.pdf' if key=='Bmass_ppRef' else 'variables_all.pdf')
        info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
        pages=int(re.search(r'Pages:\s+(\d+)',info)[1]); expected=1 if key=='Bmass_ppRef' else len(provenance['channels'][channel]['variables'])
        assert pages==expected; checks[str(pdf.relative_to(args.output))]=pages
    dump(p/'metadata/provenance.json',provenance)
    dump(p/'metadata/variables.json', {ch:{n:c for n,c in cs.items() if (n=='Bmass')==(key=='Bmass_ppRef')} for ch,cs in all_config.items()})
    dump(p/'metadata/validation.json',{'status':'PASS','selection':'1','candidate_count_and_flow_closure':'PASS','ROOT_numpy_mass_bins':'PASS','persisted_bins_and_sqrtN_errors':'PASS','source_stat_stable':'PASS','PDF_pages':{n:v for n,v in checks.items() if n.startswith(key)},'visual_review':'pending'})
    shutil.copy2(__file__,p/'metadata'/Path(__file__).name)
    shutil.copy2(source_inventory,p/'metadata/source_inventory.json')
    shutil.copy2(template,p/'metadata/aux.h')
    mass=key=='Bmass_ppRef'
    links='\n'.join(f'- [{ch}]({ch}/'+('pdf/baseline.pdf' if mass else 'variables_all.pdf')+')' for ch in ['Bu','Bs','Bd'])
    text='# '+key+'\n\n原始 ppRef DATA，selection=`1`，不追加 cut；单位候选计数（不是事件数），误差为 sqrt(N)，未拟合或减背景。\n\n'+links+'\n\n'
    text+='Bmass：三道均为严格 5 < Bmass < 6 GeV，10 MeV/bin；范围仅用于显示，不是 selection。\n\n' if mass else '逐变量图在各 channel 的 pdf/；15 bins，线性纵轴，未归一化。已有范围/绝对值显示沿用 aux.h；Btktkmass 和新分支使用各道完整有限取值范围并留边。图内注明范围外比例；显示范围不构成 cut。恒定分支保留并标注，不推定其物理有效性。Bu 原始 tree 中的第二径迹/双径迹字段同样仅作为 schema 诊断，不能据字段存在推断第二径迹的物理意义。未绘制运行/事件/候选标识及 Bmass（见相邻 Bmass_ppRef）。\n\n'
    text+='[变量配置](metadata/variables.json) · [输入与统计](metadata/provenance.json) · [数值验证](metadata/validation.json)。源文件按 ML inventory 的 size/mtime 核对并检查读前读后稳定；未重算源文件 SHA。\n'
    (p/'README.md').write_text(text)
print('ALL NUMERICAL CHECKS PASS',flush=True)
