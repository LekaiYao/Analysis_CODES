"""Two-year extended double-peak model and fixed-parameter conditional test."""
import math
import numpy as np
from chebyshev_background import chebyshev_minimum

YEARS = ('2018', '2023')
SPECIES = ('x', 'psi')
NOMINAL = {'x': 3.87169, 'psi': 3.68610}
MC_RANGES = {'x': (3.84, 3.90), 'psi': (3.65610, 3.71610)}


def quality(result):
    return (result.status() == 0 and result.covQual() == 3
            and math.isfinite(result.minNll()) and math.isfinite(result.edm())
            and result.edm() < 1.e-3)


def fit_info(result):
    def info(p):
        return {'value': p.getVal(), 'error': p.getError(),
                'range': [p.getMin(), p.getMax()],
                'at_boundary': min(p.getVal()-p.getMin(), p.getMax()-p.getVal())
                               <= 1.e-4*(p.getMax()-p.getMin())}
    return {'status': result.status(), 'covQual': result.covQual(),
            'edm': result.edm(), 'nll': result.minNll(), 'quality_pass': quality(result),
            'floating': {p.GetName(): info(p) for p in result.floatParsFinal()}}


class CombinedModel:
    def __init__(self, arrays):
        import ROOT as R
        self.R = R
        self.keep, self.results, self.attempts = [], {}, {}
        self.mass = R.RooRealVar('Bmass', 'Bmass', 3.62, 4.00)
        self.mass.setRange('all', 3.62, 4.00)
        self.category = R.RooCategory('year', 'year')
        self.data, self.mc = {}, {}
        self.arrays = arrays
        for year in YEARS:
            self.category.defineType(year)
            self.data[year] = R.RooDataSet.from_numpy(
                {'Bmass': arrays[year]['data']['Bmass']}, [self.mass], name='data_'+year)
            self.mc[year] = {}
            for k in SPECIES:
                a = arrays[year][k]
                w = R.RooRealVar('pThatreweight', 'pThatreweight', 0., max(1., float(a['pThatreweight'].max())*2))
                self.keep.append(w)
                self.mc[year][k] = R.RooDataSet.from_numpy(
                    {'Bmass': a['Bmass'], 'pThatreweight': a['pThatreweight']},
                    [self.mass, w], name='mc_'+k+'_'+year, weight_name='pThatreweight')
                if self.mc[year][k].numEntries() != len(a['Bmass']):
                    raise ValueError('MC dataset dropped rows')
                if not math.isclose(self.mc[year][k].sumEntries(), float(a['pThatreweight'].sum()), rel_tol=1.e-10):
                    raise ValueError('MC RooDataSet sumw mismatch')
        self.combined = R.RooDataSet('combined_data', 'combined_data', R.RooArgSet(self.mass),
                                    R.RooFit.Index(self.category),
                                    R.RooFit.Import('2018', self.data['2018']),
                                    R.RooFit.Import('2023', self.data['2023']))
        if self.combined.numEntries() != sum(len(arrays[y]['data']['Bmass']) for y in YEARS):
            raise ValueError('Combined dataset dropped rows')
        self.shapes, self.mc_pdfs, self.mc_shapes = {}, {}, {}
        self.means, self.scales, self.backgrounds, self.coefficients = {}, {}, {}, {}
        self.totals, self.fractions, self.yields, self.models = {}, {}, {}, {}

    def var(self, name, value, lo, hi, constant=False):
        p = self.R.RooRealVar(name, name, value, lo, hi)
        p.setConstant(constant)
        self.keep.append(p)
        return p

    def fit(self, name, pdf, data, window='all', weighted=False):
        R = self.R
        opts = [R.RooFit.Save(), R.RooFit.Range(window), R.RooFit.Strategy(2),
                R.RooFit.Hesse(True), R.RooFit.PrintLevel(-1), R.RooFit.Warnings(False),
                R.RooFit.Verbose(False), R.RooFit.Offset(False), R.RooFit.EvalBackend("legacy")]
        opts.append(R.RooFit.SumW2Error(False) if weighted else R.RooFit.Extended(True))
        r = pdf.fitTo(data, *opts)
        self.keep.append(r)
        self.attempts[name] = [r]
        if math.isfinite(r.minNll()) and not quality(r):
            pdf.getParameters(data).assignValueOnly(r.floatParsFinal())
            r = pdf.fitTo(data, *opts)
            self.keep.append(r)
            self.attempts[name].append(r)
        self.results[name] = r
        return r

    def double_gaussian(self, prefix, mean, s1, s2, frac):
        R = self.R
        g1 = R.RooGaussian(prefix+'_g1', '', self.mass, mean, s1)
        g2 = R.RooGaussian(prefix+'_g2', '', self.mass, mean, s2)
        pdf = R.RooAddPdf(prefix, '', R.RooArgList(g1, g2), R.RooArgList(frac))
        self.keep.extend([g1, g2, pdf])
        return pdf

    def fit_mc(self):
        for y in YEARS:
            for k in SPECIES:
                prefix = k+'_'+y
                nom = NOMINAL[k]
                mu = self.var(prefix+'_mc_mean', nom, nom-.01, nom+.01)
                s1 = self.var(prefix+'_sigma1', .01, .001, .1)
                s2 = self.var(prefix+'_sigma2', .005, .001, .1)
                frac = self.var(prefix+'_gauss_fraction', .5, .01, 1.)
                pdf = self.double_gaussian(prefix+'_mc', mu, s1, s2, frac)
                window = prefix+'_mc_range'
                self.mass.setRange(window, *MC_RANGES[k])
                result = self.fit(prefix+'_mc', pdf, self.mc[y][k], window, weighted=True)
                self.mc_pdfs[y, k] = pdf
                self.mc_shapes[y, k] = {'mean': mu.getVal(), 'sigma1': s1.getVal(),
                                        'sigma2': s2.getVal(), 'fraction': frac.getVal(),
                                        'range': list(MC_RANGES[k])}
                if not quality(result):
                    raise RuntimeError('MC quality gate failed: '+prefix)
                for p in (mu, s1, s2, frac):
                    p.setConstant(True)
                mean = self.var(prefix+'_mean', mu.getVal(), nom-.01, nom+.01)
                scale = self.var(prefix+'_scale', 1., .5, 1.5, constant=True)
                ws1 = self.R.RooProduct(prefix+'_scaled_sigma1', '', self.R.RooArgList(scale, s1))
                ws2 = self.R.RooProduct(prefix+'_scaled_sigma2', '', self.R.RooArgList(scale, s2))
                self.keep.extend([ws1, ws2])
                self.means[y, k], self.scales[y, k] = mean, scale
                self.shapes[y, k] = self.double_gaussian(prefix+'_signal', mean, ws1, ws2, frac)

    def build_data_model(self):
        R = self.R
        counts = {y: self.data[y].numEntries() for y in YEARS}
        n = sum(counts.values())
        for k in (*SPECIES, 'bkg'):
            if k == 'bkg':
                init = {y: .7*counts[y] for y in YEARS}
                lo, hi = .1*n, n
            else:
                near = {y: int(np.sum(abs(self.arrays[y]['data']['Bmass']-NOMINAL[k]) < .005))
                        for y in YEARS}
                init = {y: max(1., .4*near[y]) for y in YEARS}
                lo, hi = 0., sum(max(10., 2.*near[y]) for y in YEARS)
            self.totals[k] = self.var(k+'_total', sum(init.values()), lo, hi)
            self.fractions[k] = self.var(k+'_fraction_2023', init['2023']/sum(init.values()), 0., 1.)
            for y in YEARS:
                formula = '@0*@1' if y == '2023' else '@0*(1-@1)'
                v = R.RooFormulaVar(k+'_yield_'+y, formula,
                                    R.RooArgList(self.totals[k], self.fractions[k]))
                self.yields[y, k] = v
                self.keep.append(v)
        self.sim = R.RooSimultaneous('combined_model', '', self.category)
        self.null = R.RooSimultaneous('combined_null', '', self.category)
        self.null_models = {}
        for y in YEARS:
            a0 = self.var('a0_'+y, -.35, -2., 2.)
            a1 = self.var('a1_'+y, -.05, -2., 2.)
            self.coefficients[y] = [a0, a1]
            bg = R.RooChebychev('background_'+y, '', self.mass, R.RooArgList(a0, a1))
            model = R.RooAddPdf('model_'+y, '',
                               R.RooArgList(self.shapes[y, 'x'], self.shapes[y, 'psi'], bg),
                               R.RooArgList(*(self.yields[y, k] for k in (*SPECIES, 'bkg'))))
            null = R.RooAddPdf('null_'+y, '', R.RooArgList(self.shapes[y, 'psi'], bg),
                              R.RooArgList(self.yields[y, 'psi'], self.yields[y, 'bkg']))
            self.backgrounds[y], self.models[y], self.null_models[y] = bg, model, null
            self.keep.extend([bg, model, null])
            self.sim.addPdf(model, y)
            self.null.addPdf(null, y)
        self.keep.extend([self.sim, self.null])
        self.params = self.sim.getParameters(self.combined)

    def snapshot(self):
        return {p.GetName(): {'value': p.getVal(), 'constant': p.isConstant()}
                for p in self.params}

    def restore(self, snapshot):
        for p in self.params:
            p.setVal(snapshot[p.GetName()]['value'])
            p.setConstant(snapshot[p.GetName()]['constant'])

    def nll(self, model, data):
        R = self.R
        obj = model.createNLL(data, R.RooFit.Extended(True), R.RooFit.Range('all'), R.RooFit.Offset(False), R.RooFit.EvalBackend("legacy"))
        return obj.getVal()

    def array_nll(self, year, include_x=True):
        """Independent normalized density calculation, without RooFit evaluators."""
        # ROOT evaluates RooRealVar densities in double precision, even for Float_t input.
        masses = np.asarray(self.arrays[year]['data']['Bmass'], dtype=np.float64)
        density = np.zeros(len(masses))
        expected = 0.
        for k in SPECIES:
            if k == 'x' and not include_x:
                continue
            shape = self.mc_shapes[year, k]
            mu = self.means[year, k].getVal()
            components = []
            for sigma in (shape['sigma1'], shape['sigma2']):
                norm = .5*(math.erf((4.-mu)/(math.sqrt(2)*sigma))
                           - math.erf((3.62-mu)/(math.sqrt(2)*sigma)))
                components.append(np.exp(-.5*((masses-mu)/sigma)**2)/(math.sqrt(2*math.pi)*sigma*norm))
            frac = shape['fraction']
            n = self.yields[year, k].getVal()
            density += n*(frac*components[0]+(1-frac)*components[1])
            expected += n
        a0, a1 = (p.getVal() for p in self.coefficients[year])
        z = (masses-3.81)/.19
        bg = (1+a0*z+a1*(2*z*z-1))/(.38*(1-a1/3))
        n = self.yields[year, 'bkg'].getVal()
        density += n*bg
        expected += n
        if not np.isfinite(density).all() or not (density > 0).all():
            raise RuntimeError('Invalid independent intensity')
        return float(expected-np.log(density).sum())

    def check_nll(self, model, result, models):
        joint = self.nll(model, self.combined)
        separate = sum(self.nll(models[y], self.data[y]) for y in YEARS)
        independent = sum(self.array_nll(y, include_x=(model is not self.null)) for y in YEARS)
        if abs(independent-separate) > 1.e-5:
            raise RuntimeError('NumPy/analytic density NLL differs from RooFit')
        # RooSimultaneous may include a constant category term. Delta NLL must agree.
        if not math.isclose(joint, result.minNll(), abs_tol=1.e-5, rel_tol=0):
            raise RuntimeError('Fitted/recomputed NLL mismatch')
        return {'joint': joint, 'sum_year_nll': separate, 'independent_array_nll': independent,
                'category_constant': joint-separate}

    def derived_yields(self, result):
        return {y: {k: {'value': self.yields[y, k].getVal(),
                        'error': self.yields[y, k].getPropagatedError(result)}
                    for k in (*SPECIES, 'bkg')} for y in YEARS}

    def check_background(self):
        minima = {y: chebyshev_minimum([p.getVal() for p in self.coefficients[y]]) for y in YEARS}
        if any(v <= 0 for v in minima.values()):
            raise RuntimeError('Nonpositive background over full mass range')
        return minima

    def fit_full(self):
        r = self.fit('full', self.sim, self.combined)
        expected = {k+'_total' for k in (*SPECIES, 'bkg')} | {
            k+'_fraction_2023' for k in (*SPECIES, 'bkg')} | {
            k+'_'+y+'_mean' for y in YEARS for k in SPECIES} | {
            a+'_'+y for y in YEARS for a in ('a0', 'a1')}
        if {p.GetName() for p in r.floatParsFinal()} != expected:
            raise RuntimeError('Unexpected full-fit free parameters')
        if not quality(r):
            raise RuntimeError('Full fit failed quality gate; conditional test not run')
        self.full_snapshot = self.snapshot()
        self.full_yields = self.derived_yields(r)
        self.full_nll = self.check_nll(self.sim, r, self.models)
        self.background_minima = self.check_background()
        return r

    def fit_conditional(self):
        self.restore(self.full_snapshot)
        for p in self.params:
            p.setConstant(p.GetName() not in ('x_total', 'bkg_total'))
        fixed = {p.GetName(): p.getVal() for p in self.params if p.isConstant()}
        h1 = self.fit('conditional_h1', self.sim, self.combined)
        if not quality(h1):
            raise RuntimeError('Conditional H1 quality gate failed')
        self.h1_snapshot = self.snapshot()
        self.h1_yields = self.derived_yields(h1)
        nll1 = self.check_nll(self.sim, h1, self.models)
        h0 = self.fit('conditional_h0', self.null, self.combined)
        if not quality(h0):
            raise RuntimeError('Conditional H0 quality gate failed')
        nll0 = self.check_nll(self.null, h0, self.null_models)
        for name, value in fixed.items():
            if self.params.find(name).getVal() != value:
                raise RuntimeError('Conditional fixed parameter changed: '+name)
        for r, names in [(h1, {'x_total', 'bkg_total'}), (h0, {'bkg_total'})]:
            if {p.GetName() for p in r.floatParsFinal()} != names:
                raise RuntimeError('Unexpected conditional free parameters')
        self.totals['x'].setVal(0.)
        self.totals['x'].setConstant(True)
        if abs(self.nll(self.sim, self.combined)-nll0['joint']) > 1.e-5:
            raise RuntimeError('Removed X and zero X yield NLL differ')
        self.h0_snapshot = self.snapshot()
        q0 = 2*(nll0['joint']-nll1['joint'])
        if q0 < -1.e-5:
            raise RuntimeError('H0 NLL below H1: numerical failure')
        if abs(q0-2*(nll0['sum_year_nll']-nll1['sum_year_nll'])) > 1.e-5:
            raise RuntimeError('Combined q0 does not equal sum of year likelihood differences')
        if abs(nll1['joint']-self.full_nll['joint']) > .01:
            raise RuntimeError('Conditional H1 does not reproduce full-fit maximum')
        self.restore(self.h1_snapshot)
        self.conditional = {'fixed_parameters': fixed, 'H1': fit_info(h1), 'H0': fit_info(h0),
                            'nll_audit_H1': nll1, 'nll_audit_H0': nll0,
                            'q0': max(0., q0), 'significance': math.sqrt(max(0., q0)),
                            'calibrated': False,
                            'method': 'conditional likelihood; other parameters fixed from same-data full fit',
                            'audit': {'fixed_values_unchanged': True, 'removed_x_equals_zero_yield': True,
                                      'q0_equals_sum_year_differences': True}}
        return self.conditional

    def report(self):
        return {'fits': {k: fit_info(v) for k, v in self.results.items()},
                'attempts': {k: [fit_info(v) for v in rs] for k, rs in self.attempts.items()},
                'mc_shapes': {y: {k: self.mc_shapes.get((y, k)) for k in SPECIES} for y in YEARS},
                'full': {'parameters': getattr(self, 'full_snapshot', None),
                         'year_yields': getattr(self, 'full_yields', None),
                         'nll_audit': getattr(self, 'full_nll', None)},
                'conditional': getattr(self, 'conditional', None),
                'background_minima': getattr(self, 'background_minima', None)}

    def save_workspace(self, path):
        R = self.R
        self.restore(self.full_snapshot)
        w = R.RooWorkspace('ws_combined')
        imp = getattr(w, 'import')
        imp(self.combined)
        imp(self.sim)
        imp(self.null, R.RooFit.RecycleConflictNodes())
        for y in YEARS:
            imp(self.data[y])
            for k in SPECIES:
                imp(self.mc[y][k])
                imp(self.mc_pdfs[y, k], R.RooFit.RecycleConflictNodes())
        for name, snap in [('full', self.full_snapshot), ('conditional_h1', self.h1_snapshot),
                           ('conditional_h0', self.h0_snapshot)]:
            for n, v in snap.items():
                w.var(n).setVal(v['value'])
                w.var(n).setConstant(v['constant'])
            w.saveSnapshot(name, w.allVars())
        w.loadSnapshot('full')
        f = R.TFile.Open(str(path), 'RECREATE')
        if not f or f.IsZombie():
            raise OSError('Cannot create workspace')
        w.Write()
        for name, attempts in self.attempts.items():
            for i, result in enumerate(attempts):
                result.Write('fit_result_'+name+'_attempt'+str(i+1))
            attempts[-1].Write('fit_result_'+name)
        f.Close()
        self.keep.append(w)
