#!/usr/bin/env python3
"""Rebuild the two cumulative reports from completed result JSONs; never fit data."""
from pathlib import Path
import json,os
REPO=Path(__file__).resolve().parents[2]
BASELINE='X_pb23_v19_fid13_9v9_rw0_xgb_v1'
def read(p):return json.loads(p.read_text())
def link(p,owner):return os.path.relpath(p,owner.parent)
def write(p,s):
 p.parent.mkdir(parents=True,exist_ok=True)
 import tempfile, os
 fd, name = tempfile.mkstemp(prefix=p.name+'.', suffix='.tmp', dir=p.parent)
 try:
  with os.fdopen(fd, 'w') as handle: handle.write(s)
  if p.exists(): os.chmod(name, p.stat().st_mode)
  os.replace(name, p)
 finally:
  if os.path.exists(name): os.unlink(name)
def key(x):return (x[0]!=BASELINE,x[0],str(x[1]))
COMBINED_HEADING = '## PbPb2018+2023 可选 combined fit'

def combined_section():
 doc=REPO/'docs/pb23_two_peak_results.md'
 rows=[]
 for path in sorted((REPO/'fitER/results/combined_fits').glob('**/fit_result.json')):
  result=read(path)
  if result.get('workflow_revision')!='pb18_pb23_two_peak_total_fraction_v1' or result.get('status')!='PASS':continue
  manifest=read(path.parent/'manifest.json');cats=manifest['categories'];pars=result['fits']['full']['floating']
  rows.append((path,result,manifest,cats,pars))
 if not rows:return []
 lines=['',COMBINED_HEADING,'','两年独立MC峰形、DATA均值和C2背景；各成分采用总yield与2023占比。条件检验固定首次拟合的形状、ψ总yield和全部年份fraction，仅X/背景总yield自由。','',
        '| 2018 tag / 2023 tag | score18 / score23 | 总X yield（完整拟合） | f_X,2023 | 条件显著度 | 结果 |',
        '|---|---|---:|---:|---:|---|']
 for path,result,manifest,cats,pars in rows:
  x=pars['x_total'];f=pars['x_fraction_2023'];z=result['conditional']['significance']
  lines.append(f"| `{cats['2018']['tag']}` / `{cats['2023']['tag']}` | {cats['2018']['score']:.8g} / {cats['2023']['score']:.8g} | {x['value']:.2f} ± {x['error']:.2f} | {f['value']:.4f} ± {f['error']:.4f} | {z:.4f} | [DATA图册]({link(path.parent/'all_fits.pdf',doc)}) · [全部参数/审计]({link(path,doc)}) |")
 return lines

def update_combined_only():
 """Update only the optional-method section; preserve existing single-year report text."""
 import shutil,tempfile
 doc=REPO/'docs/pb23_two_peak_results.md';old=doc.read_bytes();text=old.decode()
 if text.count(COMBINED_HEADING)>1:raise ValueError('Duplicate combined report sections')
 lines=combined_section()
 if not lines:return
 new=text.split(COMBINED_HEADING)[0].rstrip()+'\n'+'\n'.join(lines)+'\n'
 if new==text:return
 backup=Path(tempfile.mkdtemp(prefix='combined_report_backup_'))/doc.name
 shutil.copy2(doc,backup)
 if doc.read_bytes()!=old:raise RuntimeError('Report changed concurrently')
 write(doc,new)
 print('Updated combined section; previous report:',backup)

def main():
 root=REPO/'fitER/results/ml_fits';wide_doc=REPO/'docs/pb23_wide_fit_results.md';two_doc=REPO/'docs/pb23_two_peak_results.md';scans=[];peaks=[]
 roots=[root,REPO/'fitER/results/archive/ml_fits']
 for tagdir in [p for base in roots if base.exists() for p in base.iterdir()]:
  if not tagdir.is_dir():continue
  for p in tagdir.iterdir():
   if not p.is_dir() or not (p/'scan_manifest.json').exists() or not (p/'fit_summary.json').exists() or not (p/'audit.json').exists():continue
   m=read(p/'scan_manifest.json')
   if m.get('mass_range')!=[3.75,4.0]:continue
   rows=read(p/'fit_summary.json');trainings=m['trainings']
   if len(trainings)!=1 or trainings[0]['tag']!=tagdir.name:continue
   expected=sorted(round(x['target_efficiency']*100) for x in trainings[0]['points'])
   if sorted(x['xeff'] for x in rows)!=expected:continue
   scans.append((tagdir.name,p,m,rows,read(p/'audit.json')))
  for p in list((tagdir/'two_peak').glob('xeff*'))+list((tagdir/'two_peak_mc_width').glob('xeff*'))+list((tagdir/'two_peak_mc_width').glob('score*'))+list((tagdir/'checks/x_width_float').glob('xeff*'))+list((tagdir/'checks/cheb4').glob('*')):
   if (p/'fit_result.json').exists():peaks.append((tagdir.name,p,read(p/'fit_result.json')))
 scans.sort(key=key);peaks.sort(key=key)
 s=['# PbPb23 宽范围单 X 拟合：持续结果汇总','','本文件由 `scripts/workflows/update_fit_reports.py` 从结果JSON生成；新结果完成后更新本文件，不另建逐次拟合docs。基准固定列在首位。','','2023 DATA，3.75–4.00 GeV；MC双高斯定形，DATA mean/scale浮动，C2背景，Strategy2。√q₀是约定的初步指标；表中不同扫描点集及fiducial selection必须一并比较。','','## 扫描最大值','','| ML tag | 扫描xeff [%] | 最大点 [%] | √q₀ | X yield | mean [GeV] | scale | 审计 |','|---|---|---:|---:|---:|---:|---:|---|'];maxima=[]
 for tag,p,m,rows,audit in scans:
  best=max(rows,key=lambda x:x['local_significance']);label=tag+(' (baseline)' if tag==BASELINE else '');points=', '.join(str(x['xeff']) for x in rows)
  s.append(f"| `{label}` | {points} | {best['xeff']} | {best['local_significance']:.4f} | {best['signal_yield']:.2f} ± {best['signal_yield_error']:.2f} | {best['mean']:.7f} | {best['width_scale']:.6f} | {audit['status']} |")
  maxima.append(dict(best,result_directory=str(p/best['label']/f"xeff{best['xeff']}"),mass_range=m['mass_range'],target_efficiencies_percent=[x['xeff'] for x in rows],source_manifest=str(p/'scan_manifest.json')))
 # Supplement scans sharing the same input/selection are also summarized together.
 for tag in sorted({item[0] for item in scans}):
  group=[item for item in scans if item[0]==tag]
  if len(group)<2:continue
  contracts={(item[2]['trainings'][0]['selection'],item[2]['trainings'][0]['weight'],tuple(item[2]['mass_range']),tuple(item[2]['width_scale_range'])) for item in group}
  if len(contracts)!=1:continue
  combined=[(row,item[1]) for item in group for row in item[3]]
  if len({row['xeff'] for row,path in combined})!=len(combined):continue
  row,path=max(combined,key=lambda pair:pair[0]['local_significance'])
  effs=', '.join(str(x) for x in sorted(row['xeff'] for row,path in combined))
  s+=['',f"同tag补点合并：`{tag}`，选点 {effs}%；当前最高√q₀={row['local_significance']:.4f}，xeff={row['xeff']}%。各批次及边界提示见下表。"]
 for tag,p,m,rows,audit in scans:
  s+=['',f'## {tag}'+('（baseline）' if tag==BASELINE else ''),'',f"预筛选：`{m['trainings'][0]['selection']}`；MC权重：`{m['trainings'][0]['weight']}`。",'', '| xeff [%] | √q₀ | X yield | mean [GeV] | scale | 边界提示 |','|---:|---:|---:|---:|---:|---|']
  for x in rows:
   notes=[]
   if x.get('mean_at_boundary'):notes.append('mean触边')
   if min(abs(x['mean']-v) for v in m.get('mean_range',[3.86169,3.88169]))<.00001 and 'mean触边' not in notes:notes.append('mean近边界')
   lo,hi=m['width_scale_range']
   if x.get('width_scale_at_boundary'):notes.append('scale触边')
   elif min(abs(x['width_scale']-lo),abs(x['width_scale']-hi))<.001:notes.append('scale近边界')
   if x.get('parameter_boundary') and not notes:notes.append('其它参数触边')
   s.append(f"| {x['xeff']} | {x['local_significance']:.4f} | {x['signal_yield']:.2f} ± {x['signal_yield_error']:.2f} | {x['mean']:.7f} | {x['width_scale']:.6f} | {'；'.join(notes) or '无mean/scale边界提示'} |")
  booklet=p/'all_fits_manifest.json'
  if booklet.exists():
   book=read(booklet);repaired=[q['xeff'] for q in book.get('pages',[]) if q.get('repaired')]
   if repaired:s+=['','图册仅包含DATA；其中'+', '.join(f'xeff{v}%' for v in repaired)+'页采用已通过单次重启的修复结果。以下扫描数值和首轮质量记录保留原样；[页面来源与修复记录]('+link(booklet,wide_doc)+')。']
  failed=[c for c in audit['checks'] if not c['quality_pass']]
  if failed:
   warnings=[]
   for check in failed:
    reasons=[f"{k}: status={v['status']}, covQual={v['covQual']}, EDM={v['edm']:.4g}" for k,v in check['fits'].items() if v['status']!=0 or v['covQual']!=3 or v['edm']>=1e-3]
    warnings.append(f"xeff{check['xeff']}%（"+'；'.join(reasons)+'）')
   valid_eff={c['xeff'] for c in audit['checks'] if c['quality_pass']}
   valid=[x for x in rows if x['xeff'] in valid_eff]
   s+=['','质量警告：'+'；'.join(warnings)+'。上表原始最大值不因质量警告而删除，也不将失败点宣称为已验收最佳点。']
   if valid:
    accepted=max(valid,key=lambda x:x['local_significance'])
    s += [f"通过数值验收的点中最高为xeff{accepted['xeff']}%，√q₀={accepted['local_significance']:.4f}，yield={accepted['signal_yield']:.2f} ± {accepted['signal_yield_error']:.2f}；边界仍需单独查看。"]
  s+=['',f"[全部拟合图]({link(p/'all_fits.pdf',wide_doc)}) · [CSV]({link(p/'fit_summary.csv',wide_doc)}) · [配置]({link(p/'scan_manifest.json',wide_doc)}) · [审计]({link(p/'audit.json',wide_doc)})。"]
 s+=['','历史窄范围与不同Strategy扫描见 [旧报告归档](archive/pb23_legacy_fits/README.md)，不混入上述宽范围比较。']
 write(wide_doc,'\n'.join(s)+'\n');write(REPO/'fitER/results/comparisons/pb23/wide_scan_maxima.json',json.dumps(maxima,indent=2)+'\n')
 s=['# PbPb 单年份双峰及条件检验：持续结果汇总','','本文件由 `scripts/workflows/update_fit_reports.py` 从结果JSON生成；以后每个明确指定的双峰选点完成后加入这里，不另建逐次拟合docs。','','DATA范围3.62–4.00 GeV。默认先拟合X+ψ(2S)+C2背景，显式C4备选单独标记；当前默认两峰mean独立浮动、scale=1固定为MC宽度；历史mean/scale浮动策略保留并明确标记。再固定最佳形状、ψ(2S) yield及背景形状，仅X/background yields自由，移除X后仅background yield自由。两步指标分别记录，条件指标的统计可靠程度留待后续讨论。','','| ML tag | xeff [%] | 浮动拟合 X yield | 浮动 √q₀,X | 条件 X yield | 条件 √q₀,X |','|---|---:|---:|---:|---:|---:|']
 summaries=[]
 for tag,p,result in peaks:
  x=result['signals']['x'];c=read(p/'conditional/fit_result.json') if (p/'conditional/fit_result.json').exists() else None;cy=c['alt']['floating']['x_yield'] if c else None
  cs=f"{cy['value']:.2f} ± {cy['error']:.2f}" if cy else '未完成';cz=f"{c['Z_PL']:.4f}" if c else '未完成';name=tag+(' (baseline)' if tag==BASELINE else '')+(' [默认 MC width]' if result.get('shape_strategy')=='mc_width_mean_float' else ' [测试：仅X scale浮动]' if result.get('shape_strategy')=='x_scale_float_psi_mc_width' else ' [历史浮动 scale]')
  name+=(' [C4背景备选]' if result.get('background_order',2)==4 else '')
  s.append(f"| `{name}` | {result['xeff']} | {x['yield']['value']:.2f} ± {x['yield']['error']:.2f} | {x['Z_PL']:.4f} | {cs} | {cz} |")
  summaries.append({'tag':tag,'xeff':result['xeff'],'floating_Z':x['Z_PL'],'conditional_Z':c['Z_PL'] if c else None,'directory':str(p),'shape_strategy':result.get('shape_strategy','mean_scale_float'),'background_order':result.get('background_order',2)})
 for tag,p,result in peaks:
  s+=['',f"## {tag} / xeff{result['xeff']}%"+('（baseline）' if tag==BASELINE else '')+(' — 默认MC宽度，scale=1' if result.get('shape_strategy')=='mc_width_mean_float' else ' — 测试：X scale浮动，ψ scale=1' if result.get('shape_strategy')=='x_scale_float_psi_mc_width' else ' — 历史scale浮动'),'','| 峰 | yield | mean [GeV] | scale | mean/scale触边 |','|---|---:|---:|---:|---|']
  s[-4]+=('，C4背景备选' if result.get('background_order',2)==4 else '')
  for k,name in [('x','X'),('psi','ψ(2S)')]:
   x=result['signals'][k];s.append(f"| {name} | {x['yield']['value']:.2f} ± {x['yield']['error']:.2f} | {x['mean']['value']:.7f} | {x['scale']['value']:.6f} | {x['mean']['at_boundary']}/{x['scale']['at_boundary']} |")
  s+=['',f"第一步数值质量：{'PASS' if result['audit']['quality_pass'] else 'WARN'}。[双峰图]({link(p/'data_fit.pdf',two_doc)}) · [完整结果]({link(p/'fit_result.json',two_doc)}) · [输入配置]({link(p/'manifest.json',two_doc)})。"]
  if result.get('background_order',2)==4:
   s+=['','本次使用四阶Chebyshev备选；a0–a3先自由拟合，条件检验固定为第一步值，常规默认仍为二阶。背景系数：'+', '.join(f"{k}={v['value']:.6g} ± {v['error']:.3g}" for k,v in result['background'].items() if k.startswith('a'))+'。']
  for stage,hp in [('浮动双峰',p/'retry_history.json'),('条件检验',p/'conditional/retry_history.json')]:
   if hp.exists():
    retries=[h for h in read(hp) if len(h['attempts'])>1]
    for h in retries:
     first,last=h['attempts'][0],h['attempts'][-1]
     s+=['',f"{stage}触发一次同设置重启：首轮status/covQual={first['status']}/{first['covQual']}、EDM={first['edm']:.6g}；最终status/covQual={last['status']}/{last['covQual']}、EDM={last['edm']:.6g}。[重启记录]({link(hp,two_doc)})。"]
  cp=p/'conditional/fit_result.json'
  if cp.exists():
   c=read(cp);s+=['',f"条件检验：q₀={c['q0']:.6f}，√q₀={c['Z_PL']:.6f}；质量{'PASS' if c['audit']['quality_pass'] else 'WARN'}。[条件拟合图]({link(p/'conditional/data_fit.pdf',two_doc)}) · [参数固定值与验收]({link(cp,two_doc)})。"]
 req=REPO/'docs/requests/4ML/open/fid6_10v9_psi2s_scored_mc.md'
 if req.exists() and not any(t=='X_pb23_v8_fid6_10v9_rw0_xgb_v1' and x['xeff']==15 for t,p,x in peaks):s+=['','## 待输入，尚无结果','',f"`X_pb23_v8_fid6_10v9_rw0_xgb_v1` / xeff15%：缺少同模型打分的ψ(2S) MC，未执行双峰或条件检验。[输入请求]({link(req,two_doc)})。"]
 s+=combined_section()
 write(two_doc,'\n'.join(s)+'\n');write(REPO/'fitER/results/comparisons/pb23/two_peak_summary.json',json.dumps(summaries,indent=2)+'\n')
 print('Updated cumulative reports:',len(scans),'wide scans;',len(peaks),'two-peak points')
if __name__=='__main__':
 import sys
 if sys.argv[1:]==['--combined-only']:update_combined_only()
 elif len(sys.argv)==1:main()
 else:raise SystemExit('Usage: update_fit_reports.py [--combined-only]')
