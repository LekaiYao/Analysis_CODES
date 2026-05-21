#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>

#include "aux/parameters.h"
#include "aux/masses.h"

void plot_Data(TString TREE ="ntmix", TString systemNAME = "PbPb23"){
    gStyle->SetOptStat(0);

    // Create a TChain and add all files from the directory
    TChain chain(Form("%s", TREE.Data()));

    if(systemNAME.Contains("PbPb23")){        //PbPb23 data
        //chain.Add("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb23_scored_DATA.root");
        chain.Add("/eos/user/k/kprince/x3872/DATA_PbPb_AANN.root");
    } else if(systemNAME.Contains("PbPb24")) {//PbPb24 data
        chain.Add("/eos/user/k/kprince/X3872_PbPb/DATA_24b_PbPb_AANN.root");
    }else if (systemNAME.Contains("PbPb")){
        //chain.Add("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_PbPb_scored_DATA.root");
        chain.Add("/eos/user/k/kprince/X3872_PbPb/DATA_PbPb_AANN.root");
        chain.Add("/eos/user/k/kprince/X3872_PbPb/DATA_24_PbPb_AANN.root");
    }else { //ppRef data
        if (TREE == "ntmix"){       chain.Add("/eos/user/k/kprince/X3872_pp_new/DATA_pp_AANN.root");}
        //if (TREE == "ntmix"){     chain.Add("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_xgboost/scored_samples/flat_ntmix_ppRef_scored_DATA.root");}
        else if (TREE == "ntKp") {  chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_DATA.root");}
        else if (TREE == "ntphi"){  chain.Add("/eos/user/h/hmarques/Analysis_CODES/plotER/Data_Bs.root");}
        else if (TREE == "ntKstar"){chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_DATA.root");}
    }

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
    canvas->SetLeftMargin(0.15); // or try 0.18 for more space

    // Define histogram parameters
    double hist_Xlow = 5;     // Minimum Bmass
    double hist_Xhigh = 5.8;  // Maximum Bmass
    if (TREE == "ntmix"){hist_Xlow = 3.6; hist_Xhigh = 4.0;}
    int nbinsmasshisto = 40;    
    double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

    TString SELECTIONcuts = "1";


    // ppRef:  abs(By) < 1.6 && Bpt > 10 && Prediction > 0.65 && BQvalue < 0.2 && 1
    // PbPb23: abs(By) < 1.6 && Bpt > 10 && Prediction > 0.90 && BQvalue < 0.15 && CentBin > 10 && CentBin < 80 && 1
    // PbPb24: abs(By) < 1.6 && Bpt > 10 && Prediction > 0.85 && BQvalue < 0.15 && CentBin > 10 && CentBin < 80 && 1
    if (TREE == "ntmix") {
        if (systemNAME.Contains("PbPb23")) {
            SELECTIONcuts = "abs(By) < 1.6 & Bpt > 10 & Prediction > 0.91 & BQvalue < 0.15 ";
        } else if (systemNAME.Contains("PbPb24")){
            SELECTIONcuts = "abs(By) < 1.6 & Bpt > 10 & Prediction > 0.91 & BQvalue < 0.2 ";

        } else if (systemNAME.Contains("PbPb"))  {
            //SELECTIONcuts = "abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.15 && CentBin > 20 && Btrk1dR < .25 && Btrk2dR < .25 && BtrkPtimb > 0.15";  // sample already has cuts applied
            SELECTIONcuts = "  abs(By) < 1.6 & Bpt > 10 & BQvalue < 0.25  & Prediction > 0.88";  
        } else {
            //SELECTIONcuts = "abs(By) < 1.6 && Bpt > 10 && Prediction > 0.65 && BQvalue < 0.2";
            SELECTIONcuts = "Prediction > 0.59 && Bpt > 10 && abs(By) < 1.6 && BQvalue < 0.15"; // 
        }
    }
    else if (TREE != "ntmix"){ SELECTIONcuts = "Bnorm_svpvDistance_2D > 4" ;}

    //SELECTIONcuts = "1";
    cout << "Applying selection cuts: " << SELECTIONcuts.Data() << std::endl;
    //std::cout << "DATA entries (after cuts): " << chain.GetEntries(SELECTIONcuts.Data()) << std::endl;

    TString Xlabel;
    if (TREE == "ntmix")      {Xlabel  = "m_{J/#Psi #pi^{+} #pi^{-}} (GeV)";} 
    else if (TREE == "ntphi") {Xlabel  = "m_{J/#Psi K^{+} K^{-}} (GeV)";    }
    else if (TREE == "ntKp")  {Xlabel  = "m_{J/#Psi K^{+}} (GeV)";          }
    else if (TREE == "ntKstar"){Xlabel = "m_{J/#Psi K^{+} #pi^{-}} (GeV)"; }

    // Create an histogram for Bmass
    TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; %s ; Entries / %.1f MeV", Xlabel.Data(), bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);

    chain.Draw("Bmass >> hist_Bmass", Form("%s ",SELECTIONcuts.Data()));
    // Customize the histogram
    hist_Bmass->SetLineColor(kBlack);
    hist_Bmass->SetLineWidth(1);
    hist_Bmass->SetFillColor(kBlack);    
    hist_Bmass->SetFillStyle(3017); 
    hist_Bmass->SetMinimum(0.0);
    hist_Bmass->Draw("");

    // Print system and number of entries in top right corner
    int nentries = hist_Bmass->GetEntries();
    TLatex* sys_entries = new TLatex(0.2, 0.88, Form("%s, N = %d", systemNAME.Data(), nentries));
	sys_entries->SetNDC();
	sys_entries->SetTextAlign(13);
	sys_entries->SetTextFont(42);
	sys_entries->SetTextSize(0.035);
	sys_entries->SetLineWidth(2);
	sys_entries->Draw();

    if (true && TREE == "ntmix") {
        const double yMin = 0.0;
        const double yMax = hist_Bmass->GetMaximum() * 1.02;

        TLine *lineX = new TLine(X3872_MASS, yMin, X3872_MASS, yMax);
        lineX->SetLineStyle(2);
        lineX->SetLineColor(kOrange-3);
        lineX->SetLineWidth(2);
        lineX->Draw("same");

        TLine *linePSI = new TLine(PSI2S_MASS, yMin, PSI2S_MASS, yMax);
        linePSI->SetLineStyle(2);
        linePSI->SetLineColor(kOrange-2);
        linePSI->SetLineWidth(2);
        linePSI->Draw("same");
    }

    gPad->Update();

    // Save the canvas as an image
    canvas->SaveAs(Form("DATA_%s_%s_Bmass.pdf", systemNAME.Data(), TREE.Data()));

    // Clean up
    delete hist_Bmass;
    delete canvas;
}
