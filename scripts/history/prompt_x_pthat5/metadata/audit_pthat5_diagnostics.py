from pathlib import Path
import ROOT as R,json,numpy as np,hashlib
R.gROOT.SetBatch(True)
out=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/plotER/results/prompt_x_pthat5');report=json.loads((out/'metadata/summary.json').read_text());provenance=json.loads((out/'metadata/provenance.json').read_text());result={}
cut='Bpt>15 && Bpt<50 && abs(By)<1.6 && BQvalue<0.15'
for year,source in [(key,value) for key,value in provenance.items() if key in ['2018','2023']]:
 file=R.TFile.Open(source['path']);tree=file.Get('ntmix_X3872');expected=report['years'][year]['groups']['all'];assert tree.GetEntries(cut)==expected['candidates']
 histogram=R.TH1D('audit_'+year,'',1,0,1);histogram.Sumw2();tree.Draw('0.5>>audit_'+year,'pThatreweight*('+cut+')','goff');assert np.isclose(histogram.GetBinContent(1),expected['sumw'],rtol=1e-10);assert np.isclose(histogram.GetBinError(1)**2,expected['sumw2'],rtol=1e-10)
 if year=='2023':
  for name,extra in [('source5','pThats==5'),('source5_actual_pthat_le10','pThats==5 && pthat<=10')]:assert tree.GetEntries(cut+' && '+extra)==report['years'][year]['groups'][name]['candidates']
 result[year]={'ROOT_selection_count_sumw_sumw2':'PASS'};file.Close()
identity=json.loads((out/'metadata/split_identity.json').read_text())
for source in identity['sources']:
 stat=Path(source['path']).stat();assert (stat.st_size,stat.st_mtime_ns)==(source['size'],source['mtime_ns'])
old=json.loads(Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/docs/reference/upstream_review_evidence/2026-09-29/pthat_layers/pthat_layers_audit.json').read_text())['raw2018']['prompt_x3872.root'];stat=Path(old['path']).stat();assert stat.st_size==old['size'] and stat.st_mtime_ns==old['mtime_ns']
result['source_unchanged']='PASS';result['2018_mapping_input']={key:old[key] for key in ['path','size','mtime_ns','uuid','branch_sha256']};result['status']='PASS'
(out/'metadata/independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS')
