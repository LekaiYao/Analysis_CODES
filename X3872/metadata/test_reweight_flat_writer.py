from pathlib import Path
import ROOT as R,numpy as np,uproot
from array import array
R.gROOT.SetBatch(True);R.gInterpreter.Declare((Path(__file__).resolve().parents[1]/'ReweightFlatTree.C').read_text())
f=R.TFile.Open('/tmp/pthat_writer_input.root','RECREATE');t=R.TTree('ntmix','kept title');value=array('f',[0]);weight=array('d',[0]);flag=array('b',[0]);t.Branch('value',value,'value/F');t.Branch('pThatreweight',weight,'pThatreweight/D');t.Branch('flag',flag,'flag/O')
for i in range(4):value[0]=i+.25;weight[0]=i+1;flag[0]=i%2;t.Fill()
t.Write();f.Close()
for drop in [False,True]:
 source=R.TFile.Open('/tmp/pthat_writer_input.root');out=R.TFile.Open('/tmp/pthat_writer_out_'+str(drop)+'.root','RECREATE');R.WriteReweightedTree(source.Get('ntmix'),out,R.std.vector('double')([0.,2.,3.,0.]),R.std.vector('float')([6.,12.,20.,8.]),R.std.vector('int')([5,10,15,10]),True,drop);out.Close();source.Close()
 with uproot.open('/tmp/pthat_writer_out_'+str(drop)+'.root') as ff:
  tt=ff['ntmix'];a=tt.arrays(library='np');mask=np.array([False,True,True,True]) if drop else np.ones(4,bool)
  assert tt.title=='kept title' and tt.classname=='TTree';assert np.array_equal(a['value'],np.arange(4,dtype=np.float32)[mask]+.25);assert np.array_equal(a['pThatreweight'],np.array([0.,2.,3.,0.])[mask]);assert np.array_equal(a['flag'],np.array([False,True,False,True])[mask]);assert tt.typenames()['pThats']=='int32_t';assert tt.typenames()['pthat']=='float'
print('Synthetic ROOT clone/weight/identity/drop tests PASS')
