"""Audit the explicitly requested six-pair study and build its single human report."""
import argparse, csv, hashlib, json, math, os, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np
import uproot


def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok, msg):
    if not ok: raise ValueError(msg)
def quality(f): return f['status'] == 0 and f['covQual'] == 3 and f['edm'] < .001


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); plan = read(args.plan)
    repo = Path(__file__).resolve().parents[2]
    doc = repo/'docs/paired_year_optimal_fit_results.md'
    out = args.output.resolve()
    require(not out.exists(), 'Comparison output exists')
    require(not doc.exists(), 'Report exists; preserve prior version before updating')
    groups = plan['groups']; singles = {}; evidence = {}; audits = {}; books = {
        'single_year_fits.pdf': [], 'single_year_conditional_fits.pdf': [],
        'combined_fits.pdf': [], 'combined_conditional_fits.pdf': [], 'mc_templates.pdf': []}
    cache = {}
    for spec in plan['singles']:
        p = Path(spec['out']); d = read(p/'fit_result.json'); c = read(p/'conditional/fit_result.json'); m = read(p/'manifest.json')
        require(d['audit']['quality_pass'] and c['audit']['quality_pass'], 'Single fit quality: '+spec['tag'])
        require(all(quality(f) for f in d['fits'].values()) and all(quality(c[k]) for k in ['alt','null']), 'Single numerical gate')
        require(d['tag'] == spec['tag'] and d['threshold'] == spec['score'], 'Tag/score mismatch')
        require(d['background_order'] == 2 and m['scale_fixed'] and d['mass_range'] == [3.62,4.], 'Single model mismatch')
        q = 2*(c['null']['nll']-c['alt']['nll'])
        require(abs(q-c['q0']) < 1e-6 and abs(math.sqrt(max(0.,q))-c['Z_PL']) < 1e-10, 'Conditional NLL audit')
        cat = spec['category']
        for k,r in m['inputs'].items():
            src=cat['sources'][k]; require((r['source'],r['size'],r['mtime_ns']) == (src['path'],src['size'],src['mtime_ns']), 'Single source mismatch')
            require(sha(r['cache']) == r['sha256'], 'Single cache hash mismatch')
            with uproot.open(r['cache']) as f: cache[spec['tag'], k] = f[r['tree']].arrays(library='np')
        row = dict(group=spec['group'],tag=spec['tag'],year=spec['year'],version=spec['version'],score=d['threshold'],
                   xeff=d['xeff'],data_entries=d['data_entries'],x_yield=d['signals']['x']['yield'],
                   conditional_q0=q,significance=c['Z_PL'],profile_significance=d['signals']['x']['Z_PL'],
                   directory=str(p),reused=spec['reuse'])
        singles[spec['tag']] = {'row':row,'full':d,'conditional':c,'manifest':m}
        for n in ['fit_result.json','conditional/fit_result.json','manifest.json']: evidence[str(p/n)] = sha(p/n)
        books['single_year_fits.pdf'].append(p/'data_fit.pdf')
        books['single_year_conditional_fits.pdf'].append(p/'conditional/data_fit.pdf')
        books['mc_templates.pdf'].extend([p/'x_mc_template_fit.pdf',p/'psi_mc_template_fit.pdf'])
    pairs = []
    for spec in plan['pairs']:
        p=Path(spec['out']);d=read(p/'fit_result.json');m=read(p/'manifest.json')
        require(d['status']=='PASS' and all(quality(f) for f in d['fits'].values()), 'Combined numerical gate')
        yr={y:singles[spec['tag'+y[2:]]] for y in ['2018','2023']}
        for y,s in yr.items():
            require(m['categories'][y]['tag']==s['row']['tag'] and m['categories'][y]['score']==s['row']['score'], 'Pair tag/cut mismatch')
            for k in ['data','x','psi']:
                r=d['inputs'][y][k];sr=s['manifest']['inputs'][k]
                require((r['path'],r['size'],r['mtime_ns'])==(sr['source'],sr['size'],sr['mtime_ns']), 'Combined source mismatch')
                fpath=p/r['cache'];require(sha(fpath)==r['cache_sha256'], 'Combined cache hash')
                with uproot.open(fpath) as f:b=f[r['tree']].arrays(library='np')
                a=cache[s['row']['tag'],k];require(set(a)==set(b) and all(np.array_equal(a[n],b[n]) for n in a), 'Single/combined selected entries differ')
                audits[spec['group']+'/'+spec['kind']+'/'+y+'/'+k]={'status':'PASS','entries':len(b['Bmass'])}
        c=d['conditional'];q=2*(c['H0']['nll']-c['H1']['nll'])
        require(abs(q-c['q0'])<1e-6 and abs(math.sqrt(max(0.,q))-c['significance'])<1e-10, 'Combined NLL/q0 audit')
        require(all(c['audit'].values()), 'Combined conditional audit')
        qi=sum(s['row']['conditional_q0'] for s in yr.values());z=math.sqrt(max(0.,qi))
        delta=d['full']['nll_audit']['sum_year_nll']-sum(s['full']['fits']['alt']['nll'] for s in yr.values())
        require(abs(delta) < .01, 'Combined and independent full-fit extrema differ: '+spec['group']+'/'+spec['kind'])
        row=dict(group=spec['group'],kind=spec['kind'],independent_q0=qi,independent_significance=z,
                 combined_q0=q,combined_significance=c['significance'],delta=c['significance']-z,
                 full_nll_difference=delta,reused=spec['reuse'],directory=str(p))
        pairs.append(row)
        for n in ['fit_result.json','manifest.json']:evidence[str(p/n)]=sha(p/n)
        for y in ['18','23']:
            books['combined_fits.pdf'].append(p/f'pb{y}_full_fit.pdf')
            books['combined_conditional_fits.pdf'].append(p/f'pb{y}_conditional_fit.pdf')
    stage=Path(tempfile.mkdtemp(prefix='paired_six_report_'));(stage/'metadata').mkdir()
    summary={'status':'PASS','single_year_results':[s['row'] for s in singles.values()], 'pairs':pairs,
             'same_selected_entries_and_weights':audits,'evidence_sha256':evidence,'calibrated':False,
             'reused':{'single_year':sum(s['row']['reused'] for s in singles.values()),'combined':sum(p['reused'] for p in pairs)}}
    (stage/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    for name,rows in [('single_year_results',summary['single_year_results']),('paired_results',pairs)]:
        with (stage/(name+'.csv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    book_sources={}
    for name,paths in books.items():
        subprocess.run(['pdfunite',*map(str,paths),str(stage/name)],check=True)
        from pypdf import PdfReader, PdfWriter
        reader=PdfReader(stage/name); writer=PdfWriter()
        for i,(page,source) in enumerate(zip(reader.pages,paths)):
            writer.add_page(page)
            single=next((s['row'] for s in singles.values() if source.is_relative_to(Path(s['row']['directory']))),None)
            pair=next((r for r in pairs if source.is_relative_to(Path(r['directory']))),None)
            label=(single['tag'] if single else pair['group']+' '+pair['kind'])+' / '+source.name
            writer.add_outline_item(label,i)
        with (stage/(name+'.tmp')).open('wb') as handle: writer.write(handle)
        os.replace(stage/(name+'.tmp'),stage/name)
        book_sources[name]=[{'page':i+1,'source':str(p),'sha256':sha(p)} for i,p in enumerate(paths)]
    (stage/'metadata/booklet_sources.json').write_text(json.dumps(book_sources,indent=2)+'\n')
    shutil.copy2(args.plan,stage/'metadata/plan.json');shutil.copy2(__file__,stage/'metadata/build_report.py')
    link=lambda p:os.path.relpath(p,doc.parent)
    lines=['# 六组2018/2023 ML新旧配对：optimal cut双峰拟合与显著度','',
      '按ML侧[配对文档](../../XGBoost/docs/reference/paired_year_ml_tags.md)执行；来源快照、SHA和实际配对保存在[运行清单]('+link(out/'metadata/plan.json')+')。全部使用各tag自己的 **training optimal**，严格`Prediction > cut`。','',
      '旧配对：2018 `v1/fid1` + 2023 `v9/fid1`；新配对：同一2018 + 2023 `v10/fid6`。新2023显式追加`Bchi2Prob > 0.1`，2018上游已包含该条件且实际输入逐条满足；没有重训或重新应用2018模型。','',
      '模型：3.62–4.00 GeV，X与ψ各自本年MC双高斯定形，scale=1，DATA均值独立浮动，背景C2。本文主表统一采用条件检验：完整拟合后固定峰形、ψ yield和背景形状，H1只浮动X/bkg yield，H0移除X只浮动bkg yield。combined另固定首次拟合得到的三个2023 fraction，H1只浮动X/bkg总yield，H0只浮动bkg总yield。','',
      '## 总显著度：新旧配对对照','',
      '直接合成为 `Z_ind = sqrt(Z_2018² + Z_2023²)`；combined为实际联合条件拟合。两个模型对背景年份分配的自由度不同，因此结果无需严格相等。','',
      '| 变量组 | 旧配对：直接合成 | 旧配对：combined | 新配对：直接合成 | 新配对：combined | 新−旧 combined |',
      '|---|---:|---:|---:|---:|---:|']
    for g in groups:
        a=next(p for p in pairs if p['group']==g and p['kind']=='old');b=next(p for p in pairs if p['group']==g and p['kind']=='new')
        lines.append(f"| {g} | {a['independent_significance']:.4f} σ | {a['combined_significance']:.4f} σ | {b['independent_significance']:.4f} σ | {b['combined_significance']:.4f} σ | {b['combined_significance']-a['combined_significance']:+.4f} σ |")
    lines+=['','数值核验说明：新3v3配对首次combined收敛于较差局部极小值（比独立拟合高0.9123 NLL）。直接复用原MC拟合与缓存，以已有独立拟合最优点作一次DATA重启后修复；表中使用修复结果，首轮证据留在该结果的artifacts/initial_local_minimum。其余配对未因本表重跑。', '', '## 单年份：每组三个tag','', '表中tag链接进入完整参数；列出的X yield来自完整双峰拟合。最后一列是另一套profile检验，仅作参考，不用于上表合成。','',
      '| 变量组 | tag | score > | DATA候选 | X yield | 条件显著度 | profile显著度 |','|---|---|---:|---:|---:|---:|---:|']
    for s in singles.values():
        r=s['row'];x=r['x_yield'];lines.append(f"| {r['group']} | [{r['tag']}]({link(Path(r['directory'])/'fit_result.json')}) | {r['score']:.8g} | {r['data_entries']} | {x['value']:.2f} ± {x['error']:.2f} | {r['significance']:.4f} σ | {r['profile_significance']:.4f} σ |")
    lines+=['','## 每个配对的似然比与图','', '| 变量组 | 配对 | 独立 q₀之和 | combined q₀ | combined−直接 [σ] | 图与完整结果 |', '|---|---|---:|---:|---:|---|']
    for r in pairs:
        p=Path(r['directory']);lines.append(f"| {r['group']} | {'新' if r['kind']=='new' else '旧'} | {r['independent_q0']:.6f} | {r['combined_q0']:.6f} | {r['delta']:+.4f} | [DATA图]({link(p/'all_fits.pdf')}) · [参数/审计]({link(p/'fit_result.json')}) |")
    lines+=['','## 图册与验收','']
    for name in books:lines.append(f'- [{name}]({link(out/name)})：{len(books[name])}页。')
    lines += ['', '图册按本文变量组顺序：单年份每组依次2018、旧2023、新2023；combined每组依次旧、新配对，每个配对依次2018、2023；MC每tag依次X、ψ。完整逐页来源及SHA在[图册来源]('+link(out/'metadata/booklet_sources.json')+')。', '',
      f"共18个单年份、12个combined结果；复用已有{summary['reused']['single_year']}个单年份和{summary['reused']['combined']}个combined拟合，其余补做。所有结果通过status=0、covQual=3、EDM<1e-3数值门槛。72份单年份/combined缓存比较均逐列一致，条件NLL/q₀独立复算通过。", '',
      'MC模板细分箱尾部残差仍需人工评估；数值收敛通过不代表峰形充分描述。各组及新旧配对使用重叠DATA，表中不同配对的显著度不能再次合成；本表不自动确定最终变量组或working point。', '',
      '[机器结果及内部审计]('+link(out/'results.json')+') · [单年份CSV]('+link(out/'single_year_results.csv')+') · [配对CSV]('+link(out/'paired_results.csv')+')']
    shutil.copytree(stage,out)
    for p in stage.rglob('*'):
        if p.is_file():require(sha(p)==sha(out/p.relative_to(stage)),'Published artifact mismatch')
    doc.write_text('\n'.join(lines)+'\n');shutil.rmtree(stage)
    print(json.dumps({'doc':str(doc),'pairs':pairs},indent=2))


if __name__=='__main__':main()
