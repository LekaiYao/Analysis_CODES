"""Regression checks for deterministic ML add/drop ownership and legacy inputs."""
from pathlib import Path
import json,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows'))
from fit_input_contract import resolve_selected_dir,load_training,sha

class InputPaths(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.repo=Path(self.tmp.name)/'Analysis_CODES';self.repo.mkdir()
        self.ml=self.repo.parent/'XGBoost';(self.ml/'configs').mkdir(parents=True)
        self.root=self.ml/'output/selected';self.root.mkdir(parents=True)
        self.tag='X_pb23_v31_fid20_7v17_rw0_xgb_v1'
    def registry(self,operation='add1'):
        name=f'fs_6v13_{operation}_X_pb23_v31_fid20_6v13_rw0_xgb_v1'
        cfg={name:{'baseline':'X_pb23_v31_fid20_6v13_rw0_xgb_v1','operation':operation,'candidates':[{'train_tag':self.tag}]}}
        return name,cfg
    def save(self,cfg):
        (self.ml/'configs/feature_selection.py').write_text('FEATURE_SELECTION_EXPERIMENTS = '+repr(cfg))
    def test_add_and_drop_owners(self):
        for operation in ['add1','drop1']:
            name,cfg=self.registry(operation);self.save(cfg);p=self.root/name/self.tag;p.mkdir(parents=True)
            self.assertEqual(resolve_selected_dir(self.repo,self.tag),p)
    def test_missing_owner_does_not_use_legacy(self):
        name,cfg=self.registry();self.save(cfg)
        with self.assertRaises(FileNotFoundError):resolve_selected_dir(self.repo,self.tag)
        (self.root/self.tag).mkdir()
        with self.assertRaises(ValueError):resolve_selected_dir(self.repo,self.tag)
    def test_duplicate_owners_rejected(self):
        name,cfg=self.registry();cfg['another_owner']=cfg[name];self.save(cfg)
        with self.assertRaises(ValueError):resolve_selected_dir(self.repo,self.tag)
    def test_legacy_schema_v2(self):
        p=self.root/self.tag;(p/'cut_scan').mkdir(parents=True);tp=p/'cut_scan/thresholds.json'
        points=[{'target_efficiency':.4,'score_threshold':.5}]
        tp.write_text(json.dumps({'weight_branch':['Reweight'],'thresholds':points}))
        category={'signal_mc':{'path':'mc.root','tree':'ntmix_X3872','event_weight_branch':'Reweight'},'data':{'path':'data.root','tree':'ntmix'},'threshold_provenance':{'sha256':sha(tp)},'fiducial_selection':{'expression':'Bpt > 10'}}
        (p/'fit_scan_manifest.pb23_pb24_simultaneous_mc_shape_nominal_v2.json').write_text(json.dumps({'pairing':{'categories':{'pb23':category}},'working_points':[{'target_weighted_efficiency':.4,'categories':{'pb23':{'threshold':.5}}}]}))
        t,actual=load_training(self.repo,self.tag);self.assertEqual(actual,points);self.assertEqual(t['weight'],'Reweight');self.assertEqual(t['selection'],'Bpt > 10');self.assertEqual(t['sources']['mc']['path'],str(p/'mc.root'))
    def test_apply_not_training_selection(self):
        p=self.root/self.tag;(p/'cut_scan').mkdir(parents=True)
        (p/'cut_scan/thresholds.json').write_text(json.dumps({'weight_branch':'pThatreweight','thresholds':[]}))
        (p/'batch_apply_summary.json').write_text(json.dumps({'input_datasets':{'dataset_year':'2023'},'input_selection':{'signal_selection':'Bpt > 7.5'},'apply_selection':{'expression':'Bpt > 15','mc_expression':'Bpt > 15 and pThatreweight < 1'},'draw_selection':{'fiducial_cut':{'expression':'Bpt > 15'}}}))
        t,_=load_training(self.repo,self.tag);self.assertEqual(t['selection'],'Bpt > 15');self.assertEqual(t['mc_apply_selection'],'Bpt > 15 and pThatreweight < 1')
if __name__=='__main__':unittest.main()
