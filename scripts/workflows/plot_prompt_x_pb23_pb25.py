#!/usr/bin/env python3
"""Compare current prompt X MC in PbPb23/25, separately for two pThat mixtures.
Run using LCG_106 x86_64-el9-gcc13-opt (ROOT 6.32.02). Never overwrite a study.
"""
import argparse
from array import array
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
import uproot
import ROOT as R

VARIABLES = ['Btrk1dR', 'Btrk2Pt', 'Btrk1Pt', 'BtrkLeadPt', 'BtrkSubPt', 'BtrkMaxdR']
MODES = ['with_pthat5', 'no_pthat5']
YEARS = ['2023', '2025']

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''):
            h.update(block)
    return h.hexdigest()

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    out = args.output.resolve()
    assert not out.exists(), out
    repo = Path(__file__).resolve().parents[2]
    manifest_path = repo / 'X3872/three_weights_pb23_pb25/metadata/mc_manifest.json'
    manifest_hash = sha(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    assert manifest['status'] == 'PASS'
    out.mkdir(parents=True)
    meta = out / 'metadata'
    meta.mkdir()
    R.gROOT.SetBatch(True)
    R.gStyle.SetOptStat(0)
    R.gStyle.SetTextFont(42)
    assert R.gROOT.GetVersion() == '6.32.02', R.gROOT.GetVersion()
    report = {'selection': 'No additional physical cut; corresponding weight > 0',
              'tree': 'ntmix_X3872', 'normalization': 'unit area in visible bins, separately by year and mode',
              'error': 'sqrt(sumw2)/visible sumw; independent-year ratio propagation, fixed normalization factors; no bin covariance',
              'ratio': '2023 / 2025', 'years': {}, 'histograms': {},
              'ranges': 'Reuse ComparePbPbYearXBaselineMC.C; derived pt/dR branches use corresponding track pt/dR ranges'}
    provenance = {'UTC': datetime.now(timezone.utc).isoformat(), 'manifest': str(manifest_path),
                  'manifest_sha256': manifest_hash, 'ROOT_version': R.gROOT.GetVersion(), 'inputs': {},
                  'template': 'plotER/validation_archive/macros/ComparePbPbYearXBaselineMC.C'}
    samples, files = {}, {}
    for year in YEARS:
        records = [r for r in manifest['files'] if r['year'] == year and r['particle'] == 'x3872' and r['kind'] == 'prompt']
        assert len(records) == 1
        record = records[0]['output']
        path = Path(record['path'])
        stat = path.stat()
        assert sha(path) == record['sha256']
        provenance['inputs'][year] = {'path': str(path), 'sha256': record['sha256'], 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
        with uproot.open(path) as f:
            a = f['ntmix_X3872'].arrays(VARIABLES + ['pThats', 'pthat'] + ['pThatreweight_' + m for m in MODES], library='np')
        assert all(np.isfinite(a[v]).all() for v in a)
        # Derived branches are consumed directly; check their definitions in this sample.
        assert np.array_equal(a['BtrkLeadPt'], np.maximum(a['Btrk1Pt'], a['Btrk2Pt']))
        assert np.array_equal(a['BtrkSubPt'], np.minimum(a['Btrk1Pt'], a['Btrk2Pt']))
        with uproot.open(path) as f:
            dr2 = f['ntmix_X3872']['Btrk2dR'].array(library='np')
        assert np.array_equal(a['BtrkMaxdR'], np.maximum(a['Btrk1dR'], dr2))
        files[year] = R.TFile.Open(str(path))
        samples[year] = a
        report['years'][year] = {'entries': len(a['pthat']), 'modes': {}}
        for mode in MODES:
            w = a['pThatreweight_' + mode]
            assert (w >= 0).all()
            assert (w[a['pthat'] < a['pThats']] == 0).all()
            if mode == 'no_pthat5':
                assert (w[a['pThats'] == 5] == 0).all()
            n = int((w > 0).sum())
            assert files[year].Get('ntmix_X3872').GetEntries('pThatreweight_' + mode + '>0') == n
            report['years'][year]['modes'][mode] = {'positive_entries': n, 'sumw': float(w.sum()), 'sumw2': float(w @ w)}
    rootout = R.TFile.Open(str(meta / 'histograms.root'), 'RECREATE')
    for mode in MODES:
        pdfdir = out / 'pdf' / mode
        pdfdir.mkdir(parents=True)
        report['histograms'][mode] = {}
        for variable in VARIABLES:
            lo, hi = (0., .450001) if 'dR' in variable else (.9, 4.500001)
            edges = np.linspace(lo, hi, 16)
            stats, histograms = {}, []
            for year, color, marker in [('2023', R.kBlue + 1, 20), ('2025', R.kRed + 1, 24)]:
                x = samples[year][variable].astype(float)
                w = samples[year]['pThatreweight_' + mode].astype(float)
                # NumPy reference includes ROOT-style underflow/overflow, right endpoint excluded.
                idx = np.searchsorted(edges, x, side='right')
                counts = np.bincount(idx, weights=w, minlength=17)
                sumw2 = np.bincount(idx, weights=w*w, minlength=17)
                name = f'{mode}_{variable}_{year}'
                R.gROOT.cd()
                h = R.TH1D(name, '', 15, array('d', edges))
                h.Sumw2()
                weight = 'pThatreweight_' + mode
                files[year].Get('ntmix_X3872').Draw(f'{variable}>>{name}', f'{weight}*({weight}>0)', 'goff')
                h.SetDirectory(0)
                np.testing.assert_allclose([h.GetBinContent(i) for i in range(17)], counts, rtol=1e-10, atol=1e-9)
                np.testing.assert_allclose([h.GetBinError(i)**2 for i in range(17)], sumw2, rtol=1e-10, atol=1e-9)
                norm = float(counts[1:16].sum())
                assert norm > 0
                stats[year] = {'edges': edges.tolist(), 'sumw_including_flow': counts.tolist(),
                               'sumw2_including_flow': sumw2.tolist(), 'visible_sumw': norm,
                               'visible_weight_fraction': norm / float(w.sum())}
                rootout.cd()
                h.Write()
                h.Scale(1/norm)
                h.SetName(name + '_unit_area')
                h.Write()
                assert abs(h.Integral() - 1) < 1e-12
                h.SetLineColor(color); h.SetMarkerColor(color); h.SetMarkerStyle(marker)
                histograms.append(h)
            h23, h25 = histograms
            ratio = h23.Clone(f'{mode}_{variable}_ratio_23_over_25')
            ratio.SetDirectory(0)
            ratio.Divide(h25)
            maxratio = 0.
            ratios = []
            for i in range(1, 16):
                a, b = h23.GetBinContent(i), h25.GetBinContent(i)
                if b <= 0:
                    ratio.SetBinContent(i, -999); ratio.SetBinError(i, 0)
                    ratios.append(None)
                    continue
                va = stats['2023']['sumw2_including_flow'][i] / stats['2023']['visible_sumw']**2
                vb = stats['2025']['sumw2_including_flow'][i] / stats['2025']['visible_sumw']**2
                error = np.sqrt(va / b**2 + a*a*vb/b**4)
                np.testing.assert_allclose([ratio.GetBinContent(i), ratio.GetBinError(i)], [a/b, error], rtol=1e-10, atol=1e-12)
                maxratio = max(maxratio, a/b + error)
                ratios.append({'value': a/b, 'error': float(error)})
            stats['ratio'] = ratios
            report['histograms'][mode][variable] = stats
            canvas = R.TCanvas('c_' + mode + variable, '', 850, 760)
            top = R.TPad('top', '', 0, .28, 1, 1)
            bottom = R.TPad('bottom', '', 0, 0, 1, .28)
            top.SetLeftMargin(.13); top.SetBottomMargin(.02)
            bottom.SetLeftMargin(.13); bottom.SetBottomMargin(.34); bottom.SetTopMargin(.02)
            top.Draw(); bottom.Draw(); top.cd()
            h23.GetYaxis().SetTitle('Unit-normalized weighted entries')
            h23.GetXaxis().SetLabelSize(0); h23.GetXaxis().SetTitle('')
            h23.SetMaximum(1.35 * max(h23.GetMaximum(), h25.GetMaximum()))
            h23.SetMinimum(0); h23.Draw('E1'); h25.Draw('E1 SAME')
            legend = R.TLegend(.57, .70, .90, .87)
            legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.032)
            legend.AddEntry(h23, 'PbPb23 prompt X MC', 'lep')
            legend.AddEntry(h25, 'PbPb25 prompt X MC', 'lep'); legend.Draw()
            note = R.TPaveText(.15, .72, .55, .88, 'NDC')
            note.SetFillStyle(0); note.SetBorderSize(0); note.SetTextAlign(12); note.SetTextSize(.027)
            note.AddText('With pThat5' if mode == 'with_pthat5' else 'Without pThat5')
            note.AddText('No additional physical cuts')
            note.AddText('Unit area in displayed range'); note.Draw()
            bottom.cd()
            ratio.GetYaxis().SetTitle('23 / 25')
            ratio.GetYaxis().SetRangeUser(0, max(2., 1.1*maxratio))
            ratio.GetYaxis().SetNdivisions(305); ratio.GetYaxis().SetTitleSize(.12)
            ratio.GetYaxis().SetLabelSize(.10); ratio.GetYaxis().SetTitleOffset(.45)
            ratio.GetXaxis().SetTitle(variable); ratio.GetXaxis().SetTitleSize(.12)
            ratio.GetXaxis().SetLabelSize(.10); ratio.Draw('E1')
            unity = R.TLine(lo, 1, hi, 1); unity.SetLineStyle(2); unity.Draw('same')
            canvas.SaveAs(str(pdfdir / (variable + '.pdf')))
            rootout.cd(); ratio.Write()
            canvas.Close()
        book = out / (mode + '_all.pdf')
        subprocess.run(['pdfunite'] + [str(pdfdir / (v + '.pdf')) for v in VARIABLES] + [str(book)], check=True)
        info = subprocess.check_output(['pdfinfo', str(book)], text=True)
        assert any(line.split() == ['Pages:', '6'] for line in info.splitlines())
    rootout.Close()
    for f in files.values():
        f.Close()
    for record in provenance['inputs'].values():
        path = Path(record['path']); stat = path.stat()
        assert (stat.st_size, stat.st_mtime_ns) == (record['size'], record['mtime_ns'])
        assert sha(path) == record['sha256']
    assert sha(manifest_path) == manifest_hash
    provenance['script_sha256'] = sha(Path(__file__))
    shutil.copy2(__file__, meta / Path(__file__).name)
    report['audit'] = {'status': 'PASS', 'checks': ['manifest input SHA256 before/after', 'ROOT selection counts',
        'ROOT TTree::Draw versus NumPy all 24 histograms including sumw2 and flows', 'unit visible area',
        'all ratio values/errors', 'derived variable definitions', 'source5 and below-threshold zeros', 'two six-page PDF books']}
    dump(meta / 'provenance.json', provenance)
    dump(meta / 'summary.json', report)
    dump(meta / 'delivery_sha256.json', {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()})
    rows = []
    for v in VARIABLES:
        cells = [f"{100*report['histograms'][m][v][y]['visible_weight_fraction']:.3f}%" for m in MODES for y in YEARS]
        rows.append('| ' + v + ' | ' + ' | '.join(cells) + ' |')
    (out / 'README.md').write_text('# PbPb23 / PbPb25 prompt X MC 变量对比\n\n'
        '[含 pThat5：六页图册](with_pthat5_all.pdf) · [不含 pThat5：六页图册](no_pthat5_all.pdf)\n\n'
        '按用户顺序绘制 Btrk1dR、Btrk2Pt、Btrk1Pt、BtrkLeadPt、BtrkSubPt、BtrkMaxdR；逐变量 PDF 在 `pdf/`。'
        '输入为 `X3872/three_weights_pb23_pb25/` 的两份 prompt X MC，读取 `ntmix_X3872`。分别使用 '
        '`pThatreweight_with_pthat5` / `pThatreweight_no_pthat5`，仅要求所用权重 > 0，不加额外物理 cut；'
        '不含 pThat5 使用已重算的混合权重，而非仅从含 pThat5 样本删除来源5。\n\n'
        '沿用历史 X MC 年份对比模板的蓝色实心点/红色空心点、上下双面板、15等宽bin；'
        'pT 范围 [0.9, 4.500001) GeV/c，dR 范围 [0, 0.450001)。新增对称变量沿用相应物理量范围。'
        '各年份在显示范围内分别归一化；下图为 pb23/pb25。误差为 sqrt(sumw2)/显示范围内sumw；'
        '比值按两年份独立传播，归一化因子固定，未计归一化带来的bin间相关；分母为0时不画点。\n\n'
        '范围外不折叠进边缘bin，也不用于图内归一化；下表为显示范围保留的加权比例。完整flow与sumw2保存在metadata。\n\n'
        '| 变量 | 含5 pb23 | 含5 pb25 | 无5 pb23 | 无5 pb25 |\n|---|---:|---:|---:|---:|\n' + '\n'.join(rows) + '\n\n'
        '[输入来源](metadata/provenance.json) · [数值与验收](metadata/summary.json) · '
        '[原始/归一化直方图和比值](metadata/histograms.root)。24组直方图独立ROOT/NumPy逐bin及sumw2核验通过，'
        '两本图册各6页。仅为候选形状比较，不是绝对产额或acceptance。\n\n'
        '复现：在仓库根目录 source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh，'
        '然后运行 `python3 scripts/workflows/plot_prompt_x_pb23_pb25.py --output <新目录>`。脚本拒绝覆盖已有结果。\n')
    print(json.dumps({'status': 'PASS', 'output': str(out), 'years': report['years']}, indent=2))

if __name__ == '__main__':
    main()
