"""Tests of the ML-to-fit selection boundary; no ROOT or fitting."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from read_ml_punzi_optimal import read_ml_optimal, digest, root_selection, selection_signature

class ReaderTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
  self.parent=Path(self.temp.name);self.repo=self.parent/'Analysis_CODES';self.repo.mkdir()
  self.tag='X_pb23_v1_fid1_2v1_rw0_xgb_v1';self.ml=self.parent/'XGBoost';self.base=self.ml/'output/selected'/self.tag
  (self.base/'punzi').mkdir(parents=True);(self.base/'cut_scan').mkdir()
  self.data='(Bpt > 15) and (abs(By) < 1.6)';self.mc=self.data+' and (pThats != 5)'
  self.write('cut_scan/thresholds.json',{'weight_branch':'pThatreweight','thresholds':[]})
  self.write('batch_apply_summary.json',{'input_datasets':{'dataset_year':'2023'},'apply_selection':{'expression':self.data,'mc_expression':self.mc},'draw_selection':{'fiducial_cut':{'expression':self.data}}})
  for name in ['DATA_with_score.root','MC_with_score.root']:(self.base/name).write_bytes(b'fixture')
  model=self.ml/'output/models'/self.tag;model.mkdir(parents=True)
  self.scan={'status':'PASS','train_tag':self.tag,'recipe_id':'fixture','recipe_sha256':'','data_apply_selection':self.data+' and '+self.data,'mc_apply_selection':self.mc,'background_training_selection':self.data,'plot_config':{},'sources':[],'scopes':{}}
  self.write('punzi/recipe.json',{'recipe_id':'fixture'});self.scan['recipe_sha256']=digest(self.base/'punzi/recipe.json')
  self.scan['apply_summary_sha256']=digest(self.base/'batch_apply_summary.json')
  for name,key in [('xgb_model.pkl','model_sha256'),('scaler.pkl','scaler_sha256'),('model_config.json','config_sha256'),('run_metadata.json','metadata_sha256')]:
   (model/name).write_bytes(b'fixture');self.scan[key]=digest(model/name)
  for scope,score,suffix in [('training',.67,'optimal'),('wide',.68,'wide_optimal')]:
   plot=self.base/'cut_scan'/('score_'+suffix+'.pdf');plot.write_bytes(b'fixture')
   self.scan['scopes'][scope]={'optimal':{'score':score,'FOM':.01},'curves':{'total':[{'score':score,'FOM':.01,'epsilon':.4}]},'mass_plot':str(plot.relative_to(self.ml)),'sidebands':[[3.8,3.85],[3.89,3.94]],'width_ratio':.6}
  self.flush()
 def write(self,name,value): (self.base/name).write_text(json.dumps(value))
 def flush(self): self.write('punzi/scan.json',self.scan)
 def read(self,scope='training'): return read_ml_optimal(self.repo,self.tag,scope)
 def test_explicit_scopes(self):
  self.assertEqual(self.read()['score_threshold'],.67);self.assertEqual(self.read('wide')['score_threshold'],.68)
  self.assertEqual(self.read()['scan_sha256'],digest(self.base/'punzi/scan.json'))
 def test_missing_scan(self):
  (self.base/'punzi/scan.json').unlink()
  with self.assertRaises(FileNotFoundError):self.read()
 def test_failed_status(self):
  self.scan['status']='FAIL';self.flush()
  with self.assertRaisesRegex(ValueError,'not PASS'):self.read()
 def test_wrong_tag(self):
  self.scan['train_tag']='another';self.flush()
  with self.assertRaisesRegex(ValueError,'tag mismatch'):self.read()
 def test_absent_scope(self):
  del self.scan['scopes']['wide'];self.flush()
  with self.assertRaisesRegex(ValueError,'scope is missing'):self.read('wide')
 def test_invalid_scope(self):
  with self.assertRaisesRegex(ValueError,'explicitly'):self.read('auto')
 def test_recipe_drift(self):
  self.write('punzi/recipe.json',{'recipe_id':'changed'})
  with self.assertRaisesRegex(ValueError,'recipe hash'):self.read()
 def test_model_drift(self):
  (self.ml/'output/models'/self.tag/'xgb_model.pkl').write_bytes(b'changed')
  with self.assertRaisesRegex(ValueError,'Stale ML provenance'):self.read()
 def test_missing_plot(self):
  (self.ml/self.scan['scopes']['training']['mass_plot']).unlink()
  with self.assertRaises(FileNotFoundError):self.read()
 def test_plot_cannot_cross_tag(self):
  self.scan['scopes']['training']['mass_plot']='output/selected/other/cut_scan/plot.pdf';self.flush()
  with self.assertRaisesRegex(ValueError,'outside'):self.read()
 def test_selection_drift(self):
  self.scan['mc_apply_selection']=self.data+' and (pthat >= 10)';self.flush()
  with self.assertRaisesRegex(ValueError,'selection mismatch'):self.read()
 def test_stored_row_consistency(self):
  self.scan['scopes']['training']['optimal']['FOM']=.02;self.flush()
  with self.assertRaisesRegex(ValueError,'stored row'):self.read()
 def test_root_translation(self):
  translated=root_selection('(Bpt > 15 and abs(By) < 1.6) or (pThats != 5)')
  self.assertIn('&&',translated);self.assertIn('||',translated);self.assertNotIn(' and ',translated)
  self.assertEqual(selection_signature(self.data),selection_signature(self.data+' and '+self.data))
  with self.assertRaises(ValueError):root_selection('arbitrary(Bpt)')

if __name__=='__main__':unittest.main()
