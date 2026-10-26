"""Audit and present independent-year quadrature versus the combined conditional test."""
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import uproot


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()


def check(condition, message):
    if not condition: raise ValueError(message)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--single18', type=Path, required=True)
    ap.add_argument('--single23', type=Path, required=True)
    ap.add_argument('--combined', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    roots = {'2018': args.single18.resolve(), '2023': args.single23.resolve()}
    combined, out = args.combined.resolve(), args.output.resolve()
    check(not out.exists(), 'Refusing to overwrite comparison output')
    c = load(combined/'fit_result.json')
    cm = load(combined/'manifest.json')
    check(c['status'] == 'PASS', 'Combined fit is not PASS')
    rows, evidence, cache_checks = {}, {}, {}
    for y, root in roots.items():
        full, cond, manifest = (load(root/p) for p in ('fit_result.json', 'conditional/fit_result.json', 'manifest.json'))
        check(full['audit']['quality_pass'] and cond['audit']['quality_pass'], 'Single-year quality failure: '+y)
        check(full['tag'] == cm['categories'][y]['tag'], 'Single/combined tag mismatch')
        check(full['threshold'] == cm['categories'][y]['score'], 'Single/combined threshold mismatch')
        check(full['mass_range'] == c['inputs'][y]['data']['mass_range'] == [3.62, 4.], 'Mass-window mismatch')
        for k in ('data', 'x', 'psi'):
            a, b = manifest['inputs'][k], c['inputs'][y][k]
            check(Path(a['source']).resolve() == Path(b['path']).resolve(), 'Different source: '+y+'/'+k)
            check((a['size'], a['mtime_ns']) == (b['size'], b['mtime_ns']), 'Source version differs')
            pa, pb = Path(a['cache']), combined/b['cache']
            check(sha(pa) == a['sha256'] and sha(pb) == b['cache_sha256'], 'Cache fingerprint mismatch')
            va = uproot.open(pa)[a['tree']].arrays(library='np')
            vb = uproot.open(pb)[b['tree']].arrays(library='np')
            check(set(va) == set(vb), 'Cached columns differ')
            check(all(np.array_equal(va[n], vb[n]) for n in va), 'Selections/weights differ: '+y+'/'+k)
            cache_checks[y+'/'+k] = {'status': 'PASS', 'entries': len(va['Bmass']), 'columns': sorted(va)}
        q = 2*(cond['null']['nll']-cond['alt']['nll'])
        check(abs(q-cond['q0']) < 1.e-6 and abs(math.sqrt(max(0., q))-cond['Z_PL']) < 1.e-10,
              'Conditional NLL/q0 mismatch')
        qp = 2*(full['fits']['null_x']['nll']-full['fits']['alt']['nll'])
        check(abs(max(0., qp)-full['signals']['x']['q0']) < 1.e-6, 'Profile NLL/q0 mismatch')
        rows[y] = {'tag': full['tag'], 'score': full['threshold'], 'data_entries': full['data_entries'],
                   'x_yield_full': full['signals']['x']['yield'],
                   'psi_yield_full': full['signals']['psi']['yield'],
                   'x_yield_conditional': cond['alt']['floating']['x_yield'],
                   'conditional_q0': q, 'conditional_significance': cond['Z_PL'],
                   'profile_q0': max(0., qp), 'profile_significance': full['signals']['x']['Z_PL'],
                   'directory': str(root)}
        for p in ('fit_result.json', 'manifest.json', 'conditional/fit_result.json'):
            evidence[str(root/p)] = sha(root/p)
    cc = c['conditional']
    qcombined = 2*(cc['H0']['nll']-cc['H1']['nll'])
    check(abs(qcombined-cc['q0']) < 1.e-6, 'Combined conditional q0 mismatch')
    quadrature_q0 = sum(r['conditional_q0'] for r in rows.values())
    profile_q0 = sum(r['profile_q0'] for r in rows.values())
    full_sum = sum(load(root/'fit_result.json')['fits']['alt']['nll'] for root in roots.values())
    norm_diff = c['full']['nll_audit']['sum_year_nll']-full_sum
    result = {'status': 'PASS', 'scope': 'training optimal; conditional tests are the main comparison',
              'years': rows, 'independent_conditional': {'q0': quadrature_q0, 'significance': math.sqrt(quadrature_q0)},
              'combined_conditional': {'q0': qcombined, 'significance': math.sqrt(max(0., qcombined)),
                                       'directory': str(combined)},
              'independent_profile_reference': {'q0': profile_q0, 'significance': math.sqrt(profile_q0)},
              'combined_full_parameters': c['fits']['full']['floating'],
              'same_selected_entries_and_weights': cache_checks,
              'combined_minus_single_full_sum_nll': norm_diff,
              'calibrated': False,
              'method_difference': 'Independent nulls float each background yield; combined null fixes year fractions and floats only the background total',
              'evidence_sha256': {**evidence, str(combined/'fit_result.json'): sha(combined/'fit_result.json'),
                                  str(combined/'manifest.json'): sha(combined/'manifest.json')}}
    stage = Path(tempfile.mkdtemp(prefix='single_combined_comparison_'))
    art = stage/'metadata';art.mkdir()
    shutil.copy2(__file__, art/Path(__file__).name)
    (stage/'significance_comparison.json').write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    csv = ['method,year,score,data_entries,q0,significance']
    for y, r in rows.items(): csv.append(f"single_conditional,{y},{r['score']},{r['data_entries']},{r['conditional_q0']},{r['conditional_significance']}")
    csv += [f"independent_conditional_quadrature,combined,,,{quadrature_q0},{math.sqrt(quadrature_q0)}",
            f"combined_conditional,combined,,,{qcombined},{math.sqrt(qcombined)}",
            f"independent_profile_reference,combined,,,{profile_q0},{math.sqrt(profile_q0)}"]
    (stage/'significance_comparison.csv').write_text('\n'.join(csv)+'\n')
    pages = [roots['2018']/'data_fit.pdf', roots['2018']/'conditional/data_fit.pdf',
             roots['2023']/'data_fit.pdf', roots['2023']/'conditional/data_fit.pdf',
             combined/'pb18_full_fit.pdf', combined/'pb23_full_fit.pdf',
             combined/'pb18_conditional_fit.pdf', combined/'pb23_conditional_fit.pdf']
    subprocess.run(['pdfunite', *map(str, pages), str(stage/'all_fits.pdf')], check=True)
    # Four combined MC templates use the exact same selected MC rows as the independent fits.
    mc_pages = [combined/f'pb{y[2:]}_{k}_mc_template_fit.pdf' for y in roots for k in ('x', 'psi')]
    subprocess.run(['pdfunite', *map(str, mc_pages), str(stage/'mc_templates.pdf')], check=True)
    (art/'booklet_sources.json').write_text(json.dumps({'DATA': [{'page': i+1, 'source': str(p), 'sha256': sha(p)} for i,p in enumerate(pages)],
                                                       'MC': [{'page': i+1, 'source': str(p), 'sha256': sha(p)} for i,p in enumerate(mc_pages)]}, indent=2)+'\n')
    def link(p): return os.path.relpath(p, out)
    text = ['# PbPb2018 + PbPb2023：3v5 optimal 双峰与显著度比较', '',
            '[全部 DATA 拟合图（8页）](all_fits.pdf) · [MC 模板图（4页）](mc_templates.pdf) · [机器结果](significance_comparison.json) · [CSV](significance_comparison.csv)', '',
            '使用各年 training optimal，严格 Prediction > cut；共同 DATA 选择15<Bpt<50、|By|<1.6、BQvalue<0.15，质量3.62–4.00 GeV。X/ψ分别用本年选后MC双高斯定形，scale=1，DATA mean独立浮动，背景C2。', '',
            '| 年份 | tag | score cut | DATA候选数 | X yield（完整拟合） |', '|---|---|---:|---:|---:|']
    for y,r in rows.items():
        v=r['x_yield_full'];text.append(f"| {y} | `{r['tag']}` | {r['score']:.2f} | {r['data_entries']} | {v['value']:.2f} ± {v['error']:.2f} |")
    text += ['', '## 主对照：条件检验', '',
             '先完整拟合，再固定峰形、ψ产额和背景形状。单年份H1只浮动X/背景yield，H0移除X只浮动背景yield；combined另固定首次拟合得到的三个年份fraction，只浮动相应总yield。', '',
             '| 方法 | q₀ | 显著度 |', '|---|---:|---:|']
    for y,r in rows.items():text.append(f"| {y} 独立拟合 | {r['conditional_q0']:.6f} | {r['conditional_significance']:.4f} σ |")
    text += [f'| 两年独立结果平方和开根号 | {quadrature_q0:.6f} | {math.sqrt(quadrature_q0):.4f} σ |',
             f'| 实际 combined fit | {qcombined:.6f} | {math.sqrt(qcombined):.4f} σ |', '',
             '平方和来自独立年份条件似然比之和。combined的年份fraction固定，背景不能在两年间独立调整；因此两种方法不要求严格相等。', '',
             '## 附：完整 profile 结果', '',
             '这一组在移除X时仍对另一峰和背景做profile，与上面的条件检验分别列出，不混用。', '',
             '| 2018 | 2023 | 平方和开根号 |', '|---:|---:|---:|',
             f"| {rows['2018']['profile_significance']:.4f} σ | {rows['2023']['profile_significance']:.4f} σ | {math.sqrt(profile_q0):.4f} σ |", '',
             '## 图和证据位置', '',
             f"- 2018：[完整双峰图]({link(roots['2018']/'data_fit.pdf')}) · [条件图]({link(roots['2018']/'conditional/data_fit.pdf')}) · [全部参数]({link(roots['2018']/'fit_result.json')})。",
             f"- 2023：[完整双峰图]({link(roots['2023']/'data_fit.pdf')}) · [条件图]({link(roots['2023']/'conditional/data_fit.pdf')}) · [全部参数]({link(roots['2023']/'fit_result.json')})。",
             f"- combined：[两年份结果入口]({link(combined/'README.md')}) · [全部参数及fraction]({link(combined/'fit_result.json')})。",
             '- DATA图册页序：1–2为2018独立完整/条件，3–4为2023独立完整/条件，5–6为combined完整拟合的2018/2023图，7–8为combined条件拟合的2018/2023图。',
             '- MC图册页序：2018 X、2018 ψ、2023 X、2023 ψ；来源与SHA在metadata/booklet_sources.json。', '',
             '六份选后缓存的source_entry、质量、Prediction及MC权重均逐值一致；NLL与q₀逐项复算通过。完整MC形状、边界、首次/重启状态及内部统计解释保留在各结果JSON。', '']
    (stage/'README.md').write_text('\n'.join(text))
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(stage, out)
    for p in stage.rglob('*'):
        if p.is_file():check(sha(p) == sha(out/p.relative_to(stage)), 'Copied comparison differs')
    shutil.rmtree(stage)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__': main()
