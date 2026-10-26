#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Add raw DATA sidebands to accepted sPlot/MC histograms without refitting.

Produce in a new staging directory; publish only after independent validation.
Original histograms supply the accepted binning, normalization and errors.
"""
from pathlib import Path
from validation_archive_paths import validation_input
import argparse, hashlib, json, shutil, subprocess, re
import numpy as np
import ROOT as R

WINDOWS = {"psi": [[3.62, 3.66], [3.71, 3.75]], "x": [[3.80, 3.85], [3.89, 3.94]]}

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(8388608), b""):
            h.update(block)
    return h.hexdigest()

def dump(p, obj):
    p.write_text(json.dumps(obj, indent=2) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=["psi2s_ppref_prompt_validation", "v27_xeff50_splot", "v27_xeff50_prompt"])
    ap.add_argument("--source", type=Path, help="Explicit validated study directory, including staging")
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    repo = Path(__file__).resolve().parents[2]
    src = a.source.resolve() if a.source else repo / "plotER/validation_archive/results" / a.study
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    assert R.gROOT.GetVersion() == "6.32.02", "Use the project ROOT 6.32.02 environment"
    R.gROOT.SetBatch(True); R.gStyle.SetOptStat(0); R.gStyle.SetOptTitle(0)
    prompt = a.study == "psi2s_ppref_prompt_validation"
    datafile = src / ("artifacts/psi2s_sweights.root" if prompt else "artifacts/data_sweights.root")
    tree = "ntmix_PSI2S_sWeight" if prompt else "ntmix_sweights"
    summary = json.loads((src / "summary.json").read_text())
    sources = {str(datafile): sha(datafile), str(src / "summary.json"): sha(src / "summary.json")}
    report = {"study": a.study, "data_selection": summary["selection"], "score_threshold": summary.get("score_threshold"), "sideband_weight": "unit raw DATA, no sWeight, no subtraction", "intervals": "closed; union of the two disjoint windows", "ratio_display": [0, 5], "refit": False, "species": {}}
    for species in (["psi"] if prompt else ["x", "psi"]):
        directory = out if prompt else out / species
        fig = directory / "pdf"; fig.mkdir(parents=True)
        metricfile = validation_input(src / ("cdf_comparison.json" if prompt else species + "/variable_metrics.json"))
        metrics = json.loads(metricfile.read_text()); metrics = sorted(metrics, key=lambda row: row["variable"]) if prompt else metrics
        names = [row["variable"] for row in metrics]
        arrays = R.RDataFrame(tree, str(datafile)).AsNumpy(["Bmass", *names])
        mass = arrays["Bmass"].astype(float)
        windows = WINDOWS[species]
        masks = [(mass >= lo) & (mass <= hi) for lo, hi in windows]
        mask = masks[0] | masks[1]
        assert mask.sum() > 0 and not np.any(masks[0] & masks[1])
        histfile = src / ("artifacts/validation_histograms.root" if prompt else species + "/histograms.root")
        sources[str(histfile)] = sha(histfile); sources[str(metricfile)] = sha(metricfile)
        histout = directory / ("artifacts/validation_histograms.root" if prompt else "histograms.root")
        histout.parent.mkdir(exist_ok=True)
        shutil.copy2(histfile, histout)
        root = R.TFile.Open(str(histout), "UPDATE")
        rows = []
        for j, rec in enumerate(metrics):
            n = rec["variable"]; sp = root.Get(n + "_data"); mc = root.Get(n + "_mc")
            assert sp and mc
            lo, hi = rec["display_min"], rec["display_max"]; bins = rec["bins"]
            assert bins == 15 and sp.GetNbinsX() == bins
            edges = np.linspace(lo, hi, bins + 1)
            values = arrays[n][mask].astype(float); finite = np.isfinite(values); values = values[finite]
            display = np.abs(values) if rec["display_abs"] else values
            inside = (display >= lo) & (display < hi)
            raw = np.histogram(display[inside], edges)[0]; norm = int(raw.sum()) or len(values)
            assert norm > 0
            sb = R.TH1D(n + "_sideband", "", bins, lo, hi); sb.SetDirectory(0); sb.Sumw2()
            for i in range(bins): sb.SetBinContent(i + 1, raw[i] / norm); sb.SetBinError(i + 1, np.sqrt(raw[i]) / norm)
            for h, color, marker in [(sp, R.kRed+1, 24), (mc, R.kOrange+7, 1), (sb, R.kBlue+1, 20)]:
                h.SetLineColor(color); h.SetMarkerColor(color); h.SetMarkerStyle(marker); h.SetLineWidth(2)
            band = mc.Clone(n + "_band"); band.SetDirectory(0); band.SetFillColorAlpha(R.kOrange+7, .3); band.SetMarkerSize(0)
            c = R.TCanvas("c"+str(j), "", 760, 650)
            top = R.TPad("top"+str(j), "", 0, .30, 1, 1); bot = R.TPad("bot"+str(j), "", 0, 0, 1, .30)
            for pad in [top, bot]: pad.SetLeftMargin(.14); pad.SetRightMargin(.04); pad.Draw()
            top.SetBottomMargin(.01); bot.SetTopMargin(.01); bot.SetBottomMargin(.33); top.cd()
            high = max(h.GetBinContent(i)+h.GetBinError(i) for h in [sp,mc,sb] for i in range(1,16))
            low = min(0., min(h.GetBinContent(i)-h.GetBinError(i) for h in [sp,mc,sb] for i in range(1,16)))
            sp.SetMaximum(max(.1,high)*1.7); sp.SetMinimum(low*1.15); sp.GetYaxis().SetTitle("Normalized entries"); sp.GetXaxis().SetLabelSize(0); sp.GetXaxis().SetTitleSize(0)
            sp.Draw("E"); band.Draw("E2 SAME"); mc.Draw("HIST SAME"); sb.Draw("E SAME"); sp.Draw("E SAME")
            leg=R.TLegend(.60,.74,.95,.91); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(.033)
            leg.AddEntry(mc,"Prompt MC","l"); leg.AddEntry(sp,"sPlot DATA","lep"); leg.AddEntry(sb,"Sideband DATA (raw)","lep"); leg.Draw()
            lab=R.TLatex(); lab.SetNDC(); lab.SetTextSize(.047)
            lab.DrawLatex(.18,.87,"ppRef #psi(2S)" if prompt else ("PbPb23 X(3872)" if species=="x" else "PbPb23 #psi(2S)"))
            lab.SetTextSize(.031); lab.DrawLatex(.18,.80,"SV-PV significance < 2" if prompt else "v27 v1, X efficiency 50%")
            distance = rec["cdf_after"] if prompt else rec["cdf_distance"]
            lab.DrawLatex(.18,.735,"CDF D(sPlot,MC) = %.3f" % distance)
            if "extra_selection" in summary: lab.DrawLatex(.18,.61,"L_{xy} < 0.01 cm")
            lab.SetTextSize(.029); lab.DrawLatex(.18,.67,"SB: [%.2f,%.2f] #cup [%.2f,%.2f] GeV" % tuple(sum(windows, [])))
            bot.cd(); axis=R.TH1D(n+"_ratio_axis_sb","",bins,lo,hi); axis.SetDirectory(0); axis.SetMinimum(0); axis.SetMaximum(5)
            axis.GetYaxis().SetTitle("DATA / MC"); axis.GetYaxis().SetTitleSize(.09); axis.GetYaxis().SetLabelSize(.08); axis.GetYaxis().SetTitleOffset(.7); axis.GetYaxis().SetNdivisions(304)
            # Reuse the original variable title from the accepted plotted configuration.
            title=n
            for definition in (repo/"plotER/validation_archive/aux.h").read_text().splitlines():
                match=re.search(r'\{"([^"\n]+)",\s*"([^"\n]+)",\s*kNBins,',definition)
                if match and match[1]==n: title=match[2].split(";")[1]
            for q in R.getSignalVars("ntmix_PSI2S" if species=="psi" else "ntmix_X3872"):
                if str(q.expr)==n: title=str(q.title).split(";")[1]
            axis.GetXaxis().SetTitle(title); axis.GetXaxis().SetTitleSize(.11); axis.GetXaxis().SetLabelSize(.10); axis.Draw("AXIS")
            graphs=[]; clips=0
            for key,h,color,marker in [("splot",sp,R.kRed+1,24),("sideband",sb,R.kBlue+1,20)]:
                g=R.TGraphErrors(); g.SetName(n+"_"+key+"_ratio"); g.SetMarkerStyle(marker); g.SetMarkerColor(color); g.SetLineColor(color)
                for i in range(1,bins+1):
                    den=mc.GetBinContent(i); num=h.GetBinContent(i)
                    if den<=0: continue
                    val=num/den; err=np.hypot(h.GetBinError(i)/den,num*mc.GetBinError(i)/den**2)
                    k=g.GetN(); g.SetPoint(k,mc.GetBinCenter(i),val); g.SetPointError(k,0,float(err)); clips+=int(val<0 or val>5)
                g.Draw("P SAME"); graphs.append(g)
            line=R.TLine(lo,1,hi,1); line.SetLineStyle(2); line.Draw()
            if clips: lab.SetTextSize(.075); lab.DrawLatex(.18,.82,f"{clips} ratio points outside [0,5]")
            c.SaveAs(str(fig/(n+".pdf"))); root.cd(); sb.Write("",R.TObject.kOverwrite)
            for g in graphs:g.Write("",R.TObject.kOverwrite)
            rows.append({"variable":n,"finite_entries":len(values),"nonfinite":int((~finite).sum()),"inrange_entries":int(raw.sum()),"normalization":norm,"raw_bins":raw.tolist(),"ratio_points_outside_0_5":clips})
            c.Close()
        root.Close()
        book = out / ("validation_all.pdf" if prompt else ("x_validation_all.pdf" if species=="x" else "psi2s_validation_all.pdf"))
        subprocess.run(["pdfunite",*[str(fig/(n+".pdf")) for n in names],str(book)],check=True)
        report["species"][species]={"windows":windows,"left_entries":int(masks[0].sum()),"right_entries":int(masks[1].sum()),"union_entries":int(mask.sum()),"variables":len(names),"metrics":rows}
    for path,digest in sources.items():assert sha(Path(path))==digest
    report["source_sha256"]=sources; report["script_sha256"]=sha(Path(__file__)); report["root_version"]=R.gROOT.GetVersion()
    report["error_convention"]="Existing Sumw2; raw sideband sqrt(N); independent DATA/MC ratio propagation; no global normalization or model covariance"
    dump(out/"sideband_metadata.json",report)
    shutil.copy2(__file__,out/"update_validation_sideband.py")
    print("COMPLETE",a.study,{s:r["union_entries"] for s,r in report["species"].items()},flush=True)

if __name__=="__main__":
    R.gInterpreter.Declare('#include "'+str(Path(__file__).resolve().parents[2]/"plotER/validation_archive/aux.h")+'"')
    main()
