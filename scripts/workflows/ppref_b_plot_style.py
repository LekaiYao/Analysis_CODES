"""Shared renderer extracted unchanged from accepted ppRef B baseline style."""
import math
import ROOT as R
labels = {"Bu":"B^{+}", "Bs":"B_{s}^{0}", "Bd":"B^{0}"}

def plot(channel, name, cfg, count, stats, destination, rootfiles, mass=False, selection_label="Raw input; no additional cuts"):
    lo,hi=cfg['range']; n=cfg['bins']
    h=R.TH1D(channel+'_'+name,'',n,lo,hi); h.SetDirectory(0); h.Sumw2()
    for i,v in enumerate(count,1): h.SetBinContent(i,int(v)); h.SetBinError(i,math.sqrt(int(v)))
    h.SetBinContent(0,stats['underflow']); h.SetBinError(0,math.sqrt(stats['underflow']))
    h.SetBinContent(n+1,stats['overflow']); h.SetBinError(n+1,math.sqrt(stats['overflow']))
    h.SetEntries(stats['finite'])
    c=R.TCanvas('canvas_'+channel+'_'+name,'',760,650)
    c.SetLeftMargin(.14); c.SetRightMargin(.04); c.SetBottomMargin(.13); c.SetTopMargin(.07)
    h.SetLineColor(R.kBlue+1); h.SetMarkerColor(R.kBlue+1); h.SetMarkerStyle(20); h.SetMarkerSize(.65); h.SetLineWidth(2)
    h.SetMinimum(0); h.SetMaximum(max(1,float(max(count)))*1.55)
    h.GetXaxis().SetTitle(cfg['title']); h.GetYaxis().SetTitle('Candidates / (10 MeV/c^{2})' if mass else 'Candidates / bin')
    h.GetYaxis().SetTitleOffset(1.5); h.Draw('E1')
    lab=R.TLatex(); lab.SetNDC(); lab.SetTextFont(42); lab.SetTextSize(.042)
    lab.DrawLatex(.18,.87,cfg.get('system_label','ppRef')+' '+labels[channel]+' DATA')
    lab.SetTextSize(.030); lab.DrawLatex(.18,.815,selection_label)
    lab.DrawLatex(.18,.77,f"N = {stats['total']:,}; in range = {int(sum(count)):,}")
    lab.SetTextSize(.025)
    fraction=(stats['underflow']+stats['overflow'])/max(1,stats['finite'])
    lab.DrawLatex(.18,.725,f"Outside range: {100*fraction:.2f}%; non-finite: {stats['nonfinite']:,}")
    if stats['min']==stats['max']: lab.DrawLatex(.18,.68,'Constant branch in this input')
    c.RedrawAxis()
    sideband_lines=[]
    if cfg.get('fiducial_label'):
        lab.SetTextSize(.026);lab.DrawLatex(.18,.675,cfg['fiducial_label'])
    if cfg.get('sidebands'):
        lab.SetTextSize(.022);lab.DrawLatex(.60,.87,'ML sideband boundaries');lab.DrawLatex(.60,.825,'(not applied as a cut)')
        for interval in cfg['sidebands']:
            for edge in interval:
                line=R.TLine(edge,0,edge,.58*h.GetMaximum())
                line.SetLineColor(cfg.get('sideband_color',R.kRed+1));line.SetLineStyle(2);line.SetLineWidth(2);line.Draw();sideband_lines.append(line)
    c.Print(str(destination)); c.Close()
    key='Bmass_ppRef' if mass else 'Variables_ppRef'; rootfiles[key].cd(); h.Write()
