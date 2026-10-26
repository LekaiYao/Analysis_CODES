"""Optional independent-year two-peak combined fit. Existing single-year workflows are unchanged."""
import argparse
import ast
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from combined_two_peak_input import REVISION, build_contract, require
from fit_input_contract import selection_expression


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def read_arrays(contract, artifacts):
    import numpy as np
    import numexpr
    import uproot
    arrays, records = {}, {}
    for year, cat in contract['categories'].items():
        arrays[year], records[year] = {}, {}
        for key, spec in cat['sources'].items():
            src = Path(spec['path'])
            before = src.stat()
            require((before.st_size, before.st_mtime_ns) == (spec['size'], spec['mtime_ns']),
                    'Input changed since preflight: '+str(src))
            checksum = digest(src)
            if 'expected_sha256' in spec:
                require(checksum == spec['expected_sha256'], 'Psi2S delivery hash mismatch')
            expression = cat['selection'] if key == 'data' else cat['mc_selection']
            expr = selection_expression(expression)
            columns = {n.id for n in ast.walk(ast.parse(expr, mode='eval')) if isinstance(n, ast.Name)}-{'abs'}
            columns |= {'Bmass', 'Prediction'}
            if contract.get('selection_equivalence'): columns.add('Bchi2Prob')
            if key != 'data': columns.add('pThatreweight')
            with uproot.open(src) as f:
                tree = f[spec['tree']]
                require(columns <= set(tree.keys()), 'Missing scored branches: '+str(columns-set(tree.keys())))
                a = tree.arrays(sorted(columns), library='np')
            if 'expected_entries' in spec:
                require(len(a['Bmass']) == spec['expected_entries'], 'Psi2S delivered entry mismatch')
            require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite input values')
            if contract.get('selection_equivalence'):
                require(bool((a['Bchi2Prob'] > .1).all()), 'Bchi2Prob alignment no longer holds')
            fid = numexpr.evaluate(expr, local_dict=a)
            require(bool(fid.all()), 'Scored input violates declared APPLY selection')
            if key != 'data': require(bool((a['pThatreweight'] > 0).all()), 'Nonpositive MC weights')
            mask = fid & (a['Prediction'] > cat['score']) & (a['Bmass'] > 3.62) & (a['Bmass'] < 4.)
            columns = ['Bmass', 'Prediction']+([] if key == 'data' else ['pThatreweight'])
            selected = {k: a[k][mask] for k in columns}
            selected['source_entry'] = np.flatnonzero(mask).astype('int64')
            require(len(selected['Bmass']) > 0, 'Empty selected category: '+year+'/'+key)
            cache = artifacts/(year+'_'+key+'_cache.root')
            with uproot.recreate(cache) as f:
                f.mktree(spec['tree'], {n: v.dtype for n, v in selected.items()})
                f[spec['tree']].extend(selected)
            with uproot.open(cache) as f: check = f[spec['tree']].arrays(library='np')
            require(all(np.array_equal(v, check[n]) for n, v in selected.items()), 'Cache round trip differs')
            after = src.stat()
            require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'Source changed while reading')
            arrays[year][key] = selected
            records[year][key] = {**spec, 'sha256': checksum, 'source_entries': len(mask),
                                  'selected_entries': int(mask.sum()), 'cache': 'artifacts/'+cache.name,
                                  'cache_sha256': digest(cache), 'score': cat['score'],
                                  'selection': expression, 'mass_range': [3.62, 4.0]}
            if key != 'data':
                records[year][key].update(sumw=float(selected['pThatreweight'].sum()),
                                         sumw2=float((selected['pThatreweight']**2).sum()))
    return arrays, records


def plots(model, stage, repo):
    R = model.R
    art = stage/'artifacts'
    for name in ('TwoPeakPlots.C', 'CombinedTwoPeakPlots.C'):
        shutil.copy2(repo/'fitER/models'/name, art/name)
    R.gSystem.SetBuildDir(str(art/'build'), True)
    require(R.gSystem.CompileMacro(str(art/'CombinedTwoPeakPlots.C'), 'k') > 0, 'Plot compilation failed')
    data_pages = []
    diagnostics = {}
    for kind, snap, result, z in [('full', model.full_snapshot, model.results['full'], math.nan),
                                  ('conditional', model.h1_snapshot, model.results['conditional_h1'],
                                   model.conditional['significance'])]:
        model.restore(snap)
        category_quality = {y: R.GetDataPlotQuality(model.data[y], model.models[y], model.mass, 76, 0)
                            for y in ('2018','2023')}
        total_chi2 = sum(q.chi2 for q in category_quality.values())
        total_points = sum(q.points for q in category_quality.values())
        nfloat = result.floatParsFinal().getSize()
        ndf = total_points-nfloat
        require(ndf > 0, 'Nonpositive combined chi-square ndf')
        diagnostics[kind] = {'chi2': total_chi2, 'participating_bins':total_points,
                             'floating_parameters':nfloat, 'ndf':ndf, 'chi2_ndf':total_chi2/ndf,
                             'years':{y:{'chi2':q.chi2,'participating_bins':q.points,
                                         'ndf':None,'note':'Year contribution; shared parameters deducted only once in combined ndf'}
                                      for y,q in category_quality.items()}}
        for y in ('2018', '2023'):
            path = stage/('pb'+y[2:]+'_'+kind+'_fit.pdf')
            R.PlotCombinedData(str(path), y, model.data[y], model.models[y],
                               model.shapes[y, 'x'], model.shapes[y, 'psi'], model.backgrounds[y],
                               model.mass, 76, result, z, total_chi2, ndf)
            data_pages.append(path)
    for y in ('2018', '2023'):
        for k, label in [('x', 'X'), ('psi', '#psi(2S)')]:
            shape = model.mc_shapes[y, k]
            R.PlotPeakMC(str(stage/('pb'+y[2:]+'_'+k+'_mc_template_fit.pdf')),
                         model.mc[y][k], model.mc_pdfs[y, k], model.mass,
                         model.results[k+'_'+y+'_mc'], shape['mean'], shape['sigma1'],
                         shape['sigma2'], shape['fraction'], *shape['range'], label)
    # All-fits booklet contains DATA pages only. MC templates stay separate.
    subprocess.run(['pdfunite', *map(str, data_pages), str(stage/'all_fits.pdf')], check=True)
    model.restore(model.full_snapshot)
    save_json(stage/'data_chi2.json', diagnostics)
    return diagnostics


def execute(contract, output, repo, *, no_report=False):
    import ROOT as R
    from combined_two_peak_model import CombinedModel
    R.gROOT.SetBatch(True)
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
    require(not output.exists(), 'Refusing to overwrite an existing result: '+str(output))
    require(shutil.which('pdfunite') is not None, 'pdfunite is required for the DATA booklet')
    stage = Path(tempfile.mkdtemp(prefix='pb18_pb23_combined_'))
    artifacts = stage/'artifacts'
    artifacts.mkdir()
    model = None
    try:
        save_json(stage/'manifest.json', contract)
        for name in ('run_pb18_pb23_combined_fit.py', 'combined_two_peak_model.py',
                     'combined_two_peak_input.py', 'fit_input_contract.py',
                     'read_ml_punzi_optimal.py', 'chebyshev_background.py'):
            shutil.copy2(repo/'scripts/workflows'/name, artifacts/name)
        arrays, records = read_arrays(contract, artifacts)
        model = CombinedModel(arrays)
        model.fit_mc()
        model.build_data_model()
        model.fit_full()
        model.fit_conditional()
        model.save_workspace(stage/'fit_workspace.root')
        report = model.report()
        report.update(status='PASS', workflow_revision=REVISION, inputs=records,
                      data_entries={y: model.data[y].numEntries() for y in ('2018', '2023')})
        report['data_chi2'] = plots(model, stage, repo)
        for cat in contract['categories'].values():
            for path, checksum in cat['provenance'].items():
                require(digest(path) == checksum, 'ML provenance changed during fit')
            for spec in cat['sources'].values():
                stat = Path(spec['path']).stat()
                require((stat.st_size, stat.st_mtime_ns) == (spec['size'], spec['mtime_ns']),
                        'Scored input changed during fit')
        contract['inputs'] = records
        contract['implementation'] = {p.name: digest(p) for p in artifacts.iterdir()
                                      if p.suffix in ('.py', '.C')}
        save_json(stage/'manifest.json', contract)
        save_json(stage/'fit_result.json', report)
        (stage/'README.md').write_text(
            '# PbPb2018 + PbPb2023 combined two-peak fit\n\n'
            'Independent year categories; per-species total yield and 2023 fraction.\n\n'
            '[All DATA fits](all_fits.pdf). MC template PDFs are separate.\n\n'
            'Conditional likelihood test: full-fit shapes, Psi2S total yield and all '
            'year fractions fixed; H1 floats X/background totals, H0 floats background total.\n\n'
            f"Combined significance: {model.conditional['significance']:.6f}.\n\n"
            'Input, parameter, fit-quality and statistical interpretation records: '
            '`manifest.json`, `fit_result.json`, `fit_workspace.root`.\n')
        output.parent.mkdir(parents=True, exist_ok=True)
        # copytree refuses an existing destination; compare every copied artifact.
        shutil.copytree(stage, output)
        for p in stage.rglob('*'):
            if p.is_file(): require(digest(p) == digest(output/p.relative_to(stage)), 'Output copy differs')
        shutil.rmtree(stage)
        if not no_report:
            subprocess.run([sys.executable, str(repo/'scripts/workflows/update_fit_reports.py'), '--combined-only'], check=True)
        print(json.dumps({'status': 'PASS', 'output': str(output),
                          'significance': model.conditional['significance']}, indent=2))
    except Exception:
        if not stage.exists():
            print('Fit artifacts already delivered at '+str(output)+'; report update failed', flush=True)
            raise
        if model is not None:
            # Preserve failed fit attempts in scratch, not as a permanent scientific result.
            (stage/'failure_diagnostics.json').write_text(json.dumps(model.report(), indent=2)+'\n')
        print('Diagnostics retained in scratch: '+str(stage), flush=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--align-bchi2prob', action='store_true', help='Verify upstream 2018 Bchi2Prob > 0.1 equals explicit 2023 cut; no input selection changes')
    parser.add_argument('--no-report', action='store_true', help='Batch caller updates report after all fits finish')
    parser.add_argument('--tag18', required=True)
    parser.add_argument('--tag23', required=True)
    parser.add_argument('--ml-optimal', choices=['training', 'wide'])
    parser.add_argument('--xeff18', type=float)
    parser.add_argument('--xeff23', type=float)
    parser.add_argument('--output', type=Path, help='New semantic result directory under fitER/results/combined_fits')
    parser.add_argument('--dry-run', action='store_true', help='Read and validate ML contract only; no output or fit')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    try:
        contract = build_contract(repo, args.tag18, args.tag23, scope=args.ml_optimal,
                                  xeff18=args.xeff18, xeff23=args.xeff23, align_bchi2prob=args.align_bchi2prob)
        if args.dry_run:
            print(json.dumps(contract, indent=2, ensure_ascii=False, allow_nan=False))
            return
        require(args.output is not None, '--output is required for execution')
        output = args.output.resolve()
        root = (repo/'fitER/results/combined_fits').resolve()
        require(output.is_relative_to(root) and output != root, 'Output must be a new child of '+str(root))
        execute(contract, output, repo, no_report=args.no_report)
    except (ValueError, OSError, KeyError, AssertionError) as exc:
        parser.exit(2, 'Combined fit input/error: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
