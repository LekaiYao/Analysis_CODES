#!/usr/bin/env python3
"""Produce three-weight MC files using explicit, recorded normalization authority.

Inputs are frozen in metadata/input_inventory.json. Keep original ROOT objects
by copying the source, then append only two Double_t branches. Validate original
columns bitwise and new weights against a separate Decimal/mask implementation.
"""
from pathlib import Path
from decimal import Decimal, localcontext
from contextlib import nullcontext
import argparse, csv, datetime, hashlib, json, os, shutil, tempfile
import numpy as np
import uproot
import ROOT

THRESHOLDS = np.array([5, 10, 15, 30, 50], dtype=np.int32)
NAMES = ('pThatreweight_no_pthat5', 'pThatreweight_with_pthat5')
RTOL = 1e-12

def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def save(p, obj):
    tmp = p.with_suffix(p.suffix + '.tmp')
    with tmp.open('w') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')
    os.replace(tmp, p)

def check_input(rec, full=False):
    p = Path(rec['source']['path']); st = p.stat()
    assert (st.st_size, st.st_mtime_ns) == (rec['source']['size'], rec['source']['mtime_ns']), p
    if full:
        assert digest(p) == rec['source']['sha256'], p

def weights(actual, source, luminosities, include5):
    if not np.isfinite(actual).all() or (actual < 0).any():
        raise ValueError('Illegal/nonfinite pthat')
    if not np.isin(source, THRESHOLDS).all():
        raise ValueError('Unknown source')
    if len(luminosities) != 5 or not np.isfinite(luminosities).all() or (luminosities <= 0).any():
        raise ValueError('Invalid normalization')
    retained = np.ones(5, dtype=bool)
    if not include5:
        retained[0] = False
    denominator = np.cumsum(np.where(retained, luminosities, 0.))
    valid = (actual >= source) & (include5 | (source != 5))
    bins = np.searchsorted(THRESHOLDS, actual, side='right') - 1
    if (bins[valid] < 0).any() or (denominator[bins[valid]] <= 0).any():
        raise ValueError('Invalid denominator')
    answer = np.zeros(len(actual), dtype=np.float64)
    answer[valid] = 1. / denominator[bins[valid]]
    assert np.isfinite(answer).all()
    return answer

def decimal_intervals(rows, include5):
    # Independent decimal arithmetic, no cumulative NumPy implementation.
    with localcontext() as ctx:
        ctx.prec = 50
        total = Decimal(0); result = []
        for row in rows:
            if include5 or int(row['pthat']) != 5:
                total += Decimal(row['n_gen']) / (Decimal(row['xsec_pb']) * Decimal(row['filter_eff']))
            result.append(float(1 / total) if total else 0.)
        return result

def independent(actual, source, interval_values, include5):
    result = np.zeros(len(actual), dtype=np.float64)
    for i, lo in enumerate([5, 10, 15, 30, 50]):
        hi = [10, 15, 30, 50, np.inf][i]
        mask = (actual >= lo) & (actual < hi) & (actual >= source)
        if not include5:
            mask &= source != 5
        result[mask] = interval_values[i]
    return result

def fixtures():
    rows = [{'pthat': str(t), 'n_gen': str(i+1), 'xsec_pb': '1', 'filter_eff': '1'} for i,t in enumerate(THRESHOLDS)]
    values = [0., 1., 100., 20.]
    for t in THRESHOLDS:
        x = np.float32(t)
        values.extend([float(np.nextafter(x, np.float32(-np.inf))), float(x), float(np.nextafter(x, np.float32(np.inf)))])
        values.extend([np.nextafter(float(t), -np.inf), float(t), np.nextafter(float(t), np.inf)])
    actual = np.tile(values, 5); source = np.repeat(THRESHOLDS, len(values))
    for include in [False, True]:
        got = weights(actual, source, np.arange(1., 6.), include)
        want = independent(actual, source, decimal_intervals(rows, include), include)
        assert np.allclose(got, want, rtol=RTOL, atol=0)
        assert np.array_equal(got == 0, (actual < source) | ((source == 5) & (not include)))
    invalid = [([np.nan], [5], np.ones(5)), ([np.inf], [5], np.ones(5)), ([-1.], [5], np.ones(5)), ([20.], [7], np.ones(5)), ([20.], [5], np.zeros(5)), ([20.], [5], np.array([1., 1., np.nan, 1., 1.]))]
    for p,s,l in invalid:
        try:
            weights(np.array(p), np.array(s), l, True)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid fixture not rejected')
    return {'status': 'PASS', 'boundary_entries_per_mode': len(actual), 'invalid_cases_rejected': len(invalid), 'rtol': RTOL, 'zero_atol': 0}

def summary(w, actual, source):
    finite = np.isfinite(w); safe = w[finite]
    sw = float(safe.sum()); sw2 = float(np.square(safe).sum())
    bins = [-np.inf, 5, 10, 15, 30, 50, np.inf]
    cells = {}
    for s in THRESHOLDS:
        cells[str(s)] = []
        for lo,hi in zip(bins[:-1], bins[1:]):
            mask = (source == s) & (actual >= lo) & (actual < hi)
            v = w[mask]
            cells[str(s)].append({'lower': str(lo), 'upper': str(hi), 'entries': int(mask.sum()), 'zero': int((v == 0).sum()), 'sumw': float(v.sum()), 'sumw2': float(np.square(v).sum())})
    return {'zero': int((w == 0).sum()), 'positive': int((w > 0).sum()), 'negative': int((w < 0).sum()), 'nonfinite': int((~finite).sum()), 'sumw': sw, 'sumw2': sw2, 'Neff': sw * sw / sw2 if sw2 else None, 'source_interval_statistics': cells}

def difference(new, old, actual, source, include5):
    changed = ~np.isclose(new, old, rtol=RTOL, atol=0)
    reasons = {}
    masks = [('excluded_source5', (source == 5) & (not include5)), ('below_own_threshold', actual < source), ('exact_threshold_endpoint', np.isin(actual, THRESHOLDS))]
    if not include5:
        masks.append(('removed_source5_luminosity_from_denominator', np.ones(len(actual), dtype=bool)))
    left = changed.copy()
    for name,mask in masks:
        reasons[name] = int((left & mask).sum()); left &= ~mask
    reasons['unexplained_or_normalization_difference'] = int(left.sum())
    assert not left.any(), reasons
    return {'changed_at_rtol_1e12': int(changed.sum()), 'exclusive_reasons_in_priority_order': reasons}

def audit(rec, output, rows=None, expected_by_mode=None):
    results = {}; oldpath = Path(__file__).resolve().parent / 'PbPb23' / Path(rec['source']['path']).name
    legacy_context = uproot.open(oldpath) if rec['year'] == '2023' else nullcontext(None)
    with uproot.open(rec['source']['path']) as original, uproot.open(output) as out, legacy_context as legacy:
        assert original.classnames() == out.classnames(), 'Original ROOT objects or cycles changed'
        for name,contract in rec['trees'].items():
            a,b = original[name],out[name]
            assert b.num_entries == a.num_entries == contract['entries'] and b.title == a.title
            types = dict(contract['schema']); types.update({n: 'double' for n in NAMES})
            assert b.typenames() == types
            for start in range(0,a.num_entries,100000):
                stop = min(start+100000,a.num_entries)
                x = a.arrays(entry_start=start,entry_stop=stop,library='np')
                y = b.arrays(list(contract['schema']),entry_start=start,entry_stop=stop,library='np')
                for branch,v in x.items():
                    assert v.dtype == y[branch].dtype and v.shape == y[branch].shape and v.tobytes() == y[branch].tobytes(), (name,branch,start)
                if legacy is not None and name in legacy:
                    old = legacy[name]
                    assert old.num_entries == a.num_entries
                    shared = [k for k in old.keys() if k in x and k != 'pThatreweight']
                    z = old.arrays(shared,entry_start=start,entry_stop=stop,library='np')
                    for branch,v in z.items():
                        assert v.dtype == x[branch].dtype and v.tobytes() == x[branch].tobytes(), ('legacy identity', name,branch,start)
            data = b.arrays(['pthat','pThats','pThatreweight',*NAMES],library='np')
            actual,source = data['pthat'],data['pThats']
            result = {'status':'PASS','entries':b.num_entries,'original_columns_bitwise':'PASS','schema_title_order':'PASS','branches':{},'new_weights':{}}
            for branch in ['pThatreweight',*NAMES]:
                result['branches'][branch] = summary(data[branch],actual,source)
            for branch,include in zip(NAMES,[False,True]):
                intervals = expected_by_mode[include] if expected_by_mode is not None else decimal_intervals(rows,include)
                expected = independent(actual,source,intervals,include)
                got = data[branch]; zero = (actual < source) | ((source == 5) & (not include))
                assert np.isfinite(got).all() and (got >= 0).all()
                assert np.array_equal(got == 0,zero)
                assert np.allclose(got,expected,rtol=RTOL,atol=0)
                pos = expected > 0
                result['new_weights'][branch] = {'decimal_recompute':'PASS','max_relative_error':float((np.abs(got[pos]-expected[pos])/expected[pos]).max()) if pos.any() else 0.,'zero_set_exact':'PASS','difference_from_Henrique':difference(got,data['pThatreweight'],actual,source,include)}
            if legacy is not None and name in legacy:
                oldw = legacy[name]['pThatreweight'].array(library='np')
                assert np.allclose(oldw,data[NAMES[0]],rtol=RTOL,atol=0)
                assert np.array_equal(oldw == 0,data[NAMES[0]] == 0)
                positive = oldw > 0
                result['legacy_no5_weight_bitwise'] = 'PASS' if np.array_equal(oldw,data[NAMES[0]]) else 'NOT_BITWISE_EQUAL: floating-point arithmetic only; numerical closure PASS'
                result['legacy_no5_weight_relative_error'] = float(np.max(np.abs(oldw[positive]-data[NAMES[0]][positive])/oldw[positive])) if positive.any() else 0.
                result['legacy_no5_weight_numerical'] = 'PASS'
                result['legacy_shared_nonweight_columns_bitwise'] = 'PASS'
            else:
                result['legacy_comparison'] = 'N/A: no old PbPb25 file' if legacy is None else 'N/A: ntEvent absent from old local files'
            results[name] = result
    rf = ROOT.TFile.Open(str(output),'READ')
    for name in rec['trees']:
        t=rf.Get(name)
        for branch,include in zip(NAMES,[False,True]):
            expression = '(pthat < pThats)' if include else '(pThats == 5 || pthat < pThats)'
            assert t.GetEntries(f'(({branch} == 0) != ({expression}))') == 0
            assert t.GetLeaf(branch).GetTypeName() == 'Double_t'
        results[name]['ROOT_zero_sets_and_types'] = 'PASS'
    rf.Close()
    return results

ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare(r"""
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <vector>
#include <stdexcept>
void AppendThreeWeights(TTree* tree, const std::vector<double>& no5, const std::vector<double>& with5) {
 if (!tree || no5.size()!=static_cast<size_t>(tree->GetEntries()) || with5.size()!=no5.size()) throw std::runtime_error("Length mismatch");
 if(tree->GetBranch("pThatreweight_no_pthat5") || tree->GetBranch("pThatreweight_with_pthat5")) throw std::runtime_error("Branch collision");
 tree->SetAutoSave(0);
 double a=0.,b=0.;
 auto ba=tree->Branch("pThatreweight_no_pthat5",&a,"pThatreweight_no_pthat5/D");
 auto bb=tree->Branch("pThatreweight_with_pthat5",&b,"pThatreweight_with_pthat5/D");
 for(size_t i=0;i<no5.size();++i){a=no5[i];b=with5[i]; if(ba->Fill()<0 || bb->Fill()<0) throw std::runtime_error("Branch fill failed");}
 if(tree->Write("",TObject::kOverwrite)<=0) throw std::runtime_error("Tree write failed");
 ba->ResetAddress(); bb->ResetAddress();
}
""")


def inverse_intervals(bands):
    """Independent high precision reference directly from original plateau weights."""
    with localcontext() as ctx:
        ctx.prec = 50
        b = [Decimal.from_float(float(w)) for w in bands]
        no5 = [0.] + [float(1 / (1 / w - 1 / b[0])) for w in b[1:]]
        return {False: no5, True: [float(w) for w in b]}

def infer_normalization(rec):
    reference = None; checks = {}
    with uproot.open(rec['source']['path']) as f:
        for name in rec['trees']:
            t=f[name];a=t.arrays(['pthat','pThats','pThatreweight'],library='np')
            p,s,w=a['pthat'],a['pThats'],a['pThatreweight']
            assert np.isfinite(p).all() and (p>=0).all() and np.isin(s,THRESHOLDS).all()
            assert np.isfinite(w).all() and (w>0).all()
            bands=[];interior_counts=[]
            for i,lo in enumerate(THRESHOLDS):
                hi=THRESHOLDS[i+1] if i<4 else np.inf
                mask=(p>lo)&(p<hi)
                assert mask.any(), ('Missing interior plateau',name,int(lo))
                candidate=float(w[mask][0]);assert np.allclose(w[mask],candidate,rtol=RTOL,atol=0)
                bands.append(candidate);interior_counts.append(int(mask.sum()))
            if reference is None:reference=np.array(bands,dtype=np.float64)
            assert np.allclose(bands,reference,rtol=RTOL,atol=0), 'Inconsistent normalization between trees'
            # Henrique currently uses right-closed boundaries and first plateau below 5.
            bins=np.maximum(0,np.searchsorted(THRESHOLDS,p,side='left')-1)
            original_prediction=reference[bins]
            err=np.abs(w-original_prediction)/w
            assert (err<=RTOL).all()
            checks[name]={'entries':t.num_entries,'original_weight_closure':'PASS','plateau_interior_counts':interior_counts,'max_relative_error':float(err.max()),'source_event_counts':{str(x):int((s==x).sum()) for x in THRESHOLDS}}
    cumulative=1/reference;luminosities=np.diff(np.r_[0.,cumulative])
    assert np.isfinite(luminosities).all() and (luminosities>0).all()
    assert np.allclose(1/np.cumsum(luminosities),reference,rtol=RTOL,atol=0)
    normalization={'method':'user-approved inversion of current Henrique five-band weights','authority':'2026-10-22 user: 用五档权重反演L,np MC只要这个版本自洽即可；据此完成pb23,pb25 MC的weight计算','original_band_weights':reference.tolist(),'thresholds':THRESHOLDS.tolist(),'luminosities':luminosities.tolist(),'L_definition':'L5=1/w5; Lj=1/wj-1/w_previous; effective contributions, not independently measured production parameters','production_Ngen_cross_section_filter_efficiency':'not inferred individually; not independently verified','three_tree_self_consistency':checks}
    if rec['year']=='2023' and rec['particle']=='psi2s' and rec['kind']=='nonprompt':
        normalization['historical_source5_discrepancy']={'old_table_Ngen':489681,'current_ntEvent_source5':checks['ntEvent']['source_event_counts']['5'],'Ngen_inferred_using_historical_sigma_eff':float(luminosities[0]*2.759e7*0.61),'resolution':'use current version self-consistent original weights as explicitly approved; historical count does not define new L'}
    return normalization,inverse_intervals(reference)

def revalidate_existing(rec, item, expectation):
    """Reuse already bitwise-audited files; check only new normalization concern."""
    path=Path(item['output']['path']);assert digest(path)==item['output']['sha256']
    checks={}
    with uproot.open(path) as f:
        for name in rec['trees']:
            a=f[name].arrays(['pthat','pThats',*NAMES],library='np');checks[name]={}
            for branch,include in zip(NAMES,[False,True]):
                want=independent(a['pthat'],a['pThats'],expectation[include],include);got=a[branch]
                assert np.allclose(got,want,rtol=RTOL,atol=0)
                assert np.array_equal(got==0,want==0)
                good=want>0;error=float(np.max(np.abs(got[good]-want[good])/want[good])) if good.any() else 0.
                checks[name][branch]={'status':'PASS','max_relative_error':error,'zero_set_exact':True}
    return checks

def resume_from_henrique_bands(output):
    meta=output/'metadata';manifest_path=meta/'mc_manifest.json'
    manifest=json.loads(manifest_path.read_text());assert manifest['status']=='PARTIAL_BLOCKED_NORMALIZATION'
    inventory=json.loads((meta/'input_inventory.json').read_text())
    saved={(f['year'],f['particle'],f['kind']):f for f in manifest['files']}
    manifest['initial_producer_and_auditor_sha256']=manifest['producer_and_auditor_sha256']
    manifest['producer_and_auditor_sha256']=digest(Path(__file__))
    manifest['schema_version']=2
    manifest['normalization_policy']='User-approved current-version five-band inversion for every file; already produced files reused only after all-row numerical equivalence verification'
    manifest['authorization_override']='2026-10-22 explicit user approval supersedes the original request restriction on inverse-derived luminosities'
    manifest['historical_normalization_table']=manifest.pop('normalization_table')
    manifest['algorithm']['formula']='1 / sum(inferred_effective_L_j for retained j with threshold_j <= pthat)'
    manifest['algorithm']['inversion']='Lj=1/wj-1/w_previous; L5=1/w5, from original weight plateaus, common to all three trees'
    manifest['status']='BUILDING_FROM_APPROVED_INVERSION';manifest['files']=[];manifest.pop('completed_UTC',None)
    manifest['boundary_tests']=fixtures()
    # Exercise inverse normalization with synthetic known L, including source5 tail.
    synthetic=np.arange(1.,6.);bands=1/np.cumsum(synthetic);expected=inverse_intervals(bands)
    for include in [False,True]:
        actual=np.tile(np.array([4.9,5.,9.9,10.,14.9,15.,29.9,30.,49.9,50.,100.]),5);source=np.repeat(THRESHOLDS,11)
        got=weights(actual,source,np.diff(np.r_[0.,1/bands]),include)
        assert np.allclose(got,independent(actual,source,expected[include],include),rtol=RTOL,atol=0)
    manifest['inverse_formula_fixture']='PASS';save(manifest_path,manifest)
    for rec in inventory['files']:
        check_input(rec,True);normalization,expected=infer_normalization(rec)
        key=(rec['year'],rec['particle'],rec['kind']);previous=saved[key]
        if previous['status']=='PASS':
            checks=revalidate_existing(rec,previous,expected)
            item=previous.copy();item['original_production_normalization']=item['normalization'];item['normalization']=normalization
            item['normalization']['existing_ROOT_reuse']='Original file unchanged; verified numerically equivalent to current inverse reference within rtol=1e-12'
            item['inversion_revalidation']=checks
            check_input(rec)
            manifest['files'].append(item);save(manifest_path,manifest)
            print('REUSED AND INVERSE-VALIDATED',Path(item['output']['path']).name,flush=True)
            continue
        assert previous['status']=='BLOCKED_NORMALIZATION' and previous['output'] is None
        luminosities=np.array(normalization['luminosities']);destination=Path(rec['output']);assert not destination.exists()
        destination.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='three_weights_inverse_',dir='/tmp/leyao') as temp:
            temporary=Path(temp)/destination.name;shutil.copy2(rec['source']['path'],temporary)
            assert digest(temporary)==rec['source']['sha256']
            target=ROOT.TFile.Open(str(temporary),'UPDATE');assert target and not target.IsZombie()
            with uproot.open(rec['source']['path']) as src:
                for name in rec['trees']:
                    a=src[name].arrays(['pthat','pThats'],library='np')
                    pair=[weights(a['pthat'],a['pThats'],luminosities,include) for include in [False,True]]
                    ROOT.AppendThreeWeights(target.Get(name),ROOT.std.vector('double')(pair[0]),ROOT.std.vector('double')(pair[1]))
            target.Close()
            verification=audit(rec,temporary,expected_by_mode=expected)
            check_input(rec,True)
            artifact_hash=digest(temporary);partial=destination.with_suffix('.root.partial');assert not partial.exists()
            shutil.copy2(temporary,partial);assert digest(partial)==artifact_hash;os.replace(partial,destination)
        item={k:rec[k] for k in ['year','particle','kind','source']}
        item.update(status='PASS',output={'path':str(destination),'size':destination.stat().st_size,'mtime_ns':destination.stat().st_mtime_ns,'sha256':artifact_hash},normalization=normalization,trees={name:{'entries':v['entries'],'schema':{**v['schema'],**{n:'double' for n in NAMES}},'title':v['title']} for name,v in rec['trees'].items()},acceptance=verification)
        manifest['files'].append(item);save(manifest_path,manifest)
        print('INVERSE PRODUCTION AND INDEPENDENT AUDIT PASS',destination.name,flush=True)
    backup=json.loads(Path(inventory['backup']).read_text())
    for b in backup['files']:
        assert digest(Path(b['source']['path']))==b['source']['sha256']
        assert digest(Path(b['backup']))==b['source']['sha256']
    manifest['old_inputs_and_backups_unchanged']='PASS';manifest['passed_files']=len(manifest['files']);manifest['blocked_files']=0
    assert manifest['passed_files']==8 and all(f['status']=='PASS' for f in manifest['files'])
    manifest['status']='PASS';manifest['completed_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();save(manifest_path,manifest)
    print('ALL EIGHT FILES PASS; 24 TREES; ORIGINALS AND BACKUPS UNCHANGED',flush=True)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True); parser.add_argument('--resume-from-henrique-bands',action='store_true'); args=parser.parse_args()
    if args.resume_from_henrique_bands:
        return resume_from_henrique_bands(args.output.resolve())
    output=args.output.resolve(); meta=output/'metadata'; inventory=json.loads((meta/'input_inventory.json').read_text())
    table=meta/'Henrique_mc_normalization_4855f01.csv'
    rows=list(csv.DictReader(line for line in table.read_text().splitlines() if not line.startswith('#')))
    manifest_path=meta/'mc_manifest.json'
    assert not manifest_path.exists(), 'Do not overwrite an existing delivery'
    manifest={'schema_version':1,'status':'BUILDING','UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'request':inventory['request'],'input_inventory_sha256':digest(meta/'input_inventory.json'),'backup_manifest':inventory['backup'],'normalization_table':{'path':str(table),'sha256':digest(table),'upstream_commit':'4855f0156052a15096e04c07acb466e658940179','validation':'Historical production table, independently checked against current original five-band weights; no inversion used'},'branches':{'pThatreweight':{'origin':'Henrique','role':'original, preserved bitwise','type':'Double_t'},NAMES[0]:{'origin':'Analysis_CODES','retained_sources':[10,15,30,50],'type':'Double_t'},NAMES[1]:{'origin':'Analysis_CODES','retained_sources':[5,10,15,30,50],'type':'Double_t'}},'algorithm':{'formula':'1 / sum(N_gen_j/(cross_section_pb_j*filter_efficiency_j) for retained j with threshold_j <= pthat)','zero_rule':'source excluded OR pthat < pThats','intervals':'[5,10), [10,15), [15,30), [30,50), [50,+inf); no5 starts at 10','equality':'retained','errors':'unknown source; negative/nonfinite pthat; missing/nonpositive/nonfinite normalization','no_rows_removed':True,'no_class_balance_or_shape_reweighting':True},'boundary_tests':fixtures(),'producer_and_auditor_sha256':digest(Path(__file__)),'files':[]}
    save(manifest_path,manifest)
    for rec in inventory['files']:
        item={'year':rec['year'],'particle':rec['particle'],'kind':rec['kind'],'source':rec['source']}
        if rec['normalization_status']=='BLOCKED':
            item.update(status='BLOCKED_NORMALIZATION',output=None,missing_evidence=rec['blocker']);manifest['files'].append(item);save(manifest_path,manifest);continue
        group=sorted([r for r in rows if (r['system'],r['particle'],r['promptness'])==('PbPb'+rec['year'][2:],rec['particle'],rec['kind'])],key=lambda r:int(r['pthat']))
        assert group==rec['normalization_rows'] and [int(r['pthat']) for r in group]==THRESHOLDS.tolist()
        luminosities=np.array([int(r['n_gen'])/(float(r['xsec_pb'])*float(r['filter_eff'])) for r in group])
        destination=Path(rec['output']);assert not destination.exists();destination.parent.mkdir(parents=True,exist_ok=True)
        check_input(rec,True)
        with tempfile.TemporaryDirectory(prefix='three_weights_',dir='/tmp/leyao') as temp:
            temporary=Path(temp)/destination.name;shutil.copy2(rec['source']['path'],temporary)
            assert digest(temporary)==rec['source']['sha256']
            target=ROOT.TFile.Open(str(temporary),'UPDATE');assert target and not target.IsZombie()
            with uproot.open(rec['source']['path']) as src:
                for name in rec['trees']:
                    a=src[name].arrays(['pthat','pThats'],library='np')
                    pair=[weights(a['pthat'],a['pThats'],luminosities,include) for include in [False,True]]
                    ROOT.AppendThreeWeights(target.Get(name),ROOT.std.vector('double')(pair[0]),ROOT.std.vector('double')(pair[1]))
            target.Close()
            verification=audit(rec,temporary,group)
            check_input(rec,True)
            artifact_hash=digest(temporary);partial=destination.with_suffix('.root.partial');assert not partial.exists()
            shutil.copy2(temporary,partial);assert digest(partial)==artifact_hash;os.replace(partial,destination)
        item.update(status='PASS',output={'path':str(destination),'size':destination.stat().st_size,'mtime_ns':destination.stat().st_mtime_ns,'sha256':artifact_hash},normalization={'rows':group,'luminosities':luminosities.tolist(),'source':'historical production CSV; current upstream five-band closure within 1e-12'},trees={name:{'entries':v['entries'],'schema':{**v['schema'],**{n:'double' for n in NAMES}},'title':v['title']} for name,v in rec['trees'].items()},acceptance=verification)
        manifest['files'].append(item);save(manifest_path,manifest)
        print('PRODUCTION AND INDEPENDENT AUDIT PASS',destination.name,flush=True)
    # Old training inputs and explicit backups remain byte-identical.
    backup=json.loads(Path(inventory['backup']).read_text())
    for b in backup['files']:
        assert digest(Path(b['source']['path']))==b['source']['sha256']
        assert digest(Path(b['backup']))==b['source']['sha256']
    manifest['old_inputs_and_backups_unchanged']='PASS'
    manifest['passed_files']=sum(r['status']=='PASS' for r in manifest['files'])
    manifest['blocked_files']=sum(r['status']=='BLOCKED_NORMALIZATION' for r in manifest['files'])
    manifest['status']='PARTIAL_BLOCKED_NORMALIZATION' if manifest['blocked_files'] else 'PASS'
    manifest['completed_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();save(manifest_path,manifest)
    print(manifest['status'],'passed',manifest['passed_files'],'blocked',manifest['blocked_files'],flush=True)

if __name__=='__main__':
    main()
