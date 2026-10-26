"""Input regression and a deterministic synthetic RooFit integration check (not calibration toys)."""
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflows'))
from combined_two_peak_input import ml_signature, build_contract
from fit_input_contract import load_training


class ContractTests(unittest.TestCase):
    def test_year_guard(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base/'cut_scan').mkdir()
            (base/'cut_scan/thresholds.json').write_text(json.dumps({'weight_branch': 'pThatreweight', 'thresholds': []}))
            (base/'batch_apply_summary.json').write_text(json.dumps({
                'input_datasets': {'dataset_year': '2018'},
                'draw_selection': {'fiducial_cut': {'expression': 'Bpt > 15'}}}))
            with patch('fit_input_contract.resolve_selected_dir', return_value=base):
                with self.assertRaisesRegex(AssertionError, 'Dataset year mismatch'):
                    load_training(base, 'tag')
                record, _ = load_training(base, 'tag', year='2018')
                self.assertEqual(record['selection'], 'Bpt > 15')

    def test_choice_must_be_explicit(self):
        with self.assertRaises(ValueError):
            build_contract(Path('/tmp'), 'a', 'b')
        with self.assertRaises(ValueError):
            build_contract(Path('/tmp'), 'a', 'b', scope='training', xeff18=40)
        with self.assertRaises(ValueError):
            build_contract(Path('/tmp'), 'a', 'b', xeff18=40)

    def test_same_configuration(self):
        c = {k: None for k in ('input_columns', 'trans_columns', 'direct_xgb_version',
             'reweight_profile', 'signal_weight_expression', 'split_policy', 'split_random_state',
             'train_fraction', 'scaler_fit_scope', 'training_signal_exclusion')}
        c.update(model_params={'max_depth': 4, 'scale_pos_weight': 2., 'n_jobs': 4},
                 signal_weight_branch='pThatreweight', signal_weight_branches=['pThatreweight'])
        m = {'class_weighting': {'method': 'scale_pos_weight', 'value': 2.},
             'signal_selection': 'Bpt > 15', 'background_selection': 'Bpt > 15'}
        d, n = copy.deepcopy(c), copy.deepcopy(m)
        d['model_params']['scale_pos_weight'] = n['class_weighting']['value'] = 3.
        self.assertEqual(ml_signature(c, m), ml_signature(d, n))
        d['model_params']['max_depth'] = 5
        self.assertNotEqual(ml_signature(c, m), ml_signature(d, n))



class AlignmentTests(unittest.TestCase):
    def test_only_documented_conjunction_is_accepted(self):
        from combined_two_peak_input import add_upstream_bchi2prob
        from read_ml_punzi_optimal import selection_signature as sig
        source = sig('Bpt > 15 and Bpt < 50')
        self.assertEqual(add_upstream_bchi2prob(source), sig('Bchi2Prob > 0.1 and Bpt > 15 and Bpt < 50'))
        self.assertNotEqual(add_upstream_bchi2prob(source), sig('Bchi2Prob > 0.1 and Bpt > 10 and Bpt < 50'))
        with self.assertRaises(ValueError): add_upstream_bchi2prob(sig('Bpt > 15 or Bpt < 50'))
        with self.assertRaises(ValueError): add_upstream_bchi2prob(sig('Bchi2Prob > 0.1 and Bpt > 15'))

    def test_all_entries_and_other_configuration_are_checked(self):
        import numpy as np
        import uproot
        from combined_two_peak_input import verify_bchi2prob_alignment
        from read_ml_punzi_optimal import selection_signature as sig
        with tempfile.TemporaryDirectory() as tmp:
            cats = {}
            for year in ('2018', '2023'):
                src = Path(tmp)/(year+'.root')
                with uproot.recreate(src) as f:
                    f.mktree('t', {'Bchi2Prob': 'float64'})
                    f['t'].extend({'Bchi2Prob': np.array([.1001, .3, .8])})
                st = src.stat();sel = 'Bpt > 15 and Bpt < 50'+(' and Bchi2Prob > 0.1' if year == '2023' else '')
                cats[year] = {'selection': sel, 'mc_selection': sel+' and pThatreweight > 0',
                    'ml_signature': {'features': ['a'], 'signal_selection': sig(sel), 'background_selection': sig(sel)},
                    'sources': {'data': {'path': str(src), 'tree': 't', 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}}}
            self.assertFalse(verify_bchi2prob_alignment(cats)['selection_changed'])
            bad = copy.deepcopy(cats);bad['2023']['ml_signature']['features'] = ['b']
            with self.assertRaisesRegex(ValueError, 'only the redundant'): verify_bchi2prob_alignment(bad)
            src = Path(cats['2018']['sources']['data']['path'])
            with uproot.recreate(src) as f:
                f.mktree('t', {'Bchi2Prob': 'float64'});f['t'].extend({'Bchi2Prob': np.array([.1, .3, .8])})
            st=src.stat();cats['2018']['sources']['data'].update(size=st.st_size, mtime_ns=st.st_mtime_ns)
            with self.assertRaisesRegex(ValueError, 'not redundant'): verify_bchi2prob_alignment(cats)


class PrecisionTests(unittest.TestCase):
    def test_float_input_uses_double_density_arithmetic(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest('Run this numerical regression in LCG_106')
        from types import SimpleNamespace
        from combined_two_peak_model import CombinedModel
        def value(v): return SimpleNamespace(getVal=lambda: v)
        m = CombinedModel.__new__(CombinedModel)
        m.mc_shapes = {('2018', k): {'sigma1': .008, 'sigma2': .004, 'fraction': .25}
                       for k in ('x', 'psi')}
        m.means = {('2018', 'x'): value(3.87169), ('2018', 'psi'): value(3.68610)}
        m.yields = {('2018', k): value(v) for k, v in [('x', 100), ('psi', 200), ('bkg', 1000)]}
        m.coefficients = {'2018': [value(.4), value(.01)]}
        mass = np.linspace(3.63, 3.99, 111, dtype=np.float32)
        m.arrays = {'2018': {'data': {'Bmass': mass}}}
        n32 = m.array_nll('2018')
        m.arrays['2018']['data']['Bmass'] = mass.astype(np.float64)
        self.assertEqual(n32, m.array_nll('2018'))


class ReportTests(unittest.TestCase):
    def test_combined_section_preserves_single_year(self):
        import update_fit_reports as reports
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            doc = repo/'docs/pb23_two_peak_results.md'
            doc.parent.mkdir()
            original = '# Existing report\n\nSingle-year science remains unchanged.\n'
            doc.write_text(original)
            result_dir = repo/'fitER/results/combined_fits/example'
            result_dir.mkdir(parents=True)
            (result_dir/'fit_result.json').write_text(json.dumps({
                'workflow_revision': 'pb18_pb23_two_peak_total_fraction_v1', 'status': 'PASS',
                'fits': {'full': {'floating': {'x_total': {'value': 50., 'error': 8.},
                                              'x_fraction_2023': {'value': .4, 'error': .1}}}},
                'conditional': {'significance': 3.}}))
            (result_dir/'manifest.json').write_text(json.dumps({'categories': {
                '2018': {'tag': 'year18', 'score': .8}, '2023': {'tag': 'year23', 'score': .7}}}))
            with patch.object(reports, 'REPO', repo):
                reports.update_combined_only()
                first = doc.read_text()
                self.assertTrue(first.startswith(original))
                self.assertIn('year18', first)
                reports.update_combined_only()
                self.assertEqual(first, doc.read_text())


def cache_integration():
    import numpy as np
    import uproot
    from run_pb18_pb23_combined_fit import read_arrays, digest
    with tempfile.TemporaryDirectory(prefix='combined_cache_test_') as tmp:
        root = Path(tmp)
        cats = {}
        for y in ('2018', '2023'):
            specs = {}
            for k in ('data', 'x', 'psi'):
                src = root/(y+'_'+k+'.root')
                a = {'Bmass': np.array([3.62, 3.63, 3.7, 3.85, 4.]),
                     'Prediction': np.array([.9, .5, .6, .7, .9]), 'Bpt': np.full(5, 20.)}
                if k != 'data': a['pThatreweight'] = np.arange(1., 6.)
                with uproot.recreate(src) as f:
                    f.mktree('ntmix', {n: v.dtype for n, v in a.items()})
                    f['ntmix'].extend(a)
                stat = src.stat()
                specs[k] = {'path': str(src), 'tree': 'ntmix', 'size': stat.st_size,
                            'mtime_ns': stat.st_mtime_ns}
                if k == 'psi': specs[k].update(expected_entries=5, expected_sha256=digest(src))
            cats[y] = {'sources': specs, 'selection': 'Bpt > 15',
                       'mc_selection': 'Bpt > 15 and pThatreweight > 0', 'score': .5}
        out = root/'artifacts'
        out.mkdir()
        arrays, records = read_arrays({'categories': cats}, out)
        for y in cats:
            for k in ('data', 'x', 'psi'):
                assert arrays[y][k]['source_entry'].tolist() == [2, 3]
                assert records[y][k]['selected_entries'] == 2
        bad = copy.deepcopy(cats)
        bad['2018']['sources']['psi']['expected_sha256'] = 'stale'
        out2 = root/'rejected'
        out2.mkdir()
        try:
            read_arrays({'categories': bad}, out2)
        except ValueError as exc:
            assert 'hash mismatch' in str(exc)
        else:
            raise AssertionError('Stale Psi2S hash was accepted')
    print('CACHE INTEGRATION PASS: strict score/mass boundaries, provenance, weights, round trip')


def synthetic_integration(output):
    import numpy as np
    import ROOT as R
    from combined_two_peak_model import CombinedModel, quality
    from run_pb18_pb23_combined_fit import plots, save_json
    R.gROOT.SetBatch(True)
    R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
    rng = np.random.default_rng(30092026)
    def peak(n, mu, factor):
        sigma = np.where(rng.random(n) < .58, .004*factor, .010*factor)
        return rng.normal(mu, sigma)
    def bg(n):
        out = []
        while len(out) < n:
            z = rng.uniform(-1, 1, 2*n)
            keep = rng.random(2*n)*1.4 < 1-.25*z+.05*(2*z*z-1)
            out.extend((3.81+.19*z[keep]).tolist())
        return np.array(out[:n])
    arrays = {}
    for year, nx, np_, nb, scale in [('2018', 550, 1000, 8000, .9), ('2023', 950, 1800, 12000, 1.1)]:
        arrays[year] = {'data': {'Bmass': np.r_[peak(nx, 3.87169, scale), peak(np_, 3.68610, scale), bg(nb)]}}
        for k, mu in [('x', 3.87169), ('psi', 3.68610)]:
            mass = peak(16000, mu, scale)
            arrays[year][k] = {'Bmass': mass, 'pThatreweight': rng.uniform(.8, 1.2, len(mass))}
    m = CombinedModel(arrays)
    m.fit_mc()
    m.build_data_model()
    m.fit_full()
    m.fit_conditional()
    assert m.conditional['q0'] > 0
    for k in ('x', 'psi', 'bkg'):
        m.restore(m.full_snapshot)
        assert abs(sum(m.yields[y, k].getVal() for y in ('2018', '2023'))-m.totals[k].getVal()) < 1.e-9
    # Independent year-yield parameterization must reach the same full-fit optimum.
    direct_nll = 0.
    m.restore(m.full_snapshot)
    for y in ('2018', '2023'):
        yy = [m.var('direct_'+k+'_'+y, m.yields[y, k].getVal(), 0., 2*m.data[y].numEntries())
              for k in ('x', 'psi', 'bkg')]
        pdf = R.RooAddPdf('direct_'+y, '', R.RooArgList(m.shapes[y, 'x'], m.shapes[y, 'psi'], m.backgrounds[y]),
                         R.RooArgList(*yy))
        m.keep.append(pdf)
        result = m.fit('direct_'+y, pdf, m.data[y])
        assert quality(result)
        direct_nll += result.minNll()
    assert abs(direct_nll-m.full_nll['sum_year_nll']) < .02, (direct_nll, m.full_nll)
    m.restore(m.full_snapshot)
    output.mkdir(parents=True, exist_ok=False)
    (output/'artifacts').mkdir()
    m.save_workspace(output/'fit_workspace.root')
    plots(m, output, Path(__file__).resolve().parents[2])
    with R.TFile.Open(str(output/'fit_workspace.root')) as f:
        w = f.Get('ws_combined')
        for name, expected in [('full', 14), ('conditional_h1', 2)]:
            w.loadSnapshot(name)
            pars = w.pdf('combined_model').getParameters(w.data('combined_data'))
            assert sum(not p.isConstant() for p in pars) == expected
        w.loadSnapshot('conditional_h0')
        pars = w.pdf('combined_null').getParameters(w.data('combined_data'))
        assert {p.GetName() for p in pars if not p.isConstant()} == {'bkg_total'}
    save_json(output/'validation.json', {'status': 'PASS', 'fixture': 'deterministic synthetic, not physics result or calibration',
              'direct_year_nll_difference': direct_nll-m.full_nll['sum_year_nll'],
              'result': m.report()})
    print('ROOT INTEGRATION PASS', output)


if __name__ == '__main__':
    if sys.argv[1:] == ['--cache-check']:
        cache_integration()
    elif len(sys.argv) == 3 and sys.argv[1] == '--root-output':
        synthetic_integration(Path(sys.argv[2]))
    else:
        unittest.main()
