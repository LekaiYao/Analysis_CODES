#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>

// Very simple 2D plotting macro following plot_dataMC logic
// Usage examples (from plotER folder):
//   root -l -b -q 'plot2D_dataMC.C'                                   // defaults
//   root -l -b -q 'plot2D_dataMC.C("Bmass","Bnorm_svpvDistance_2D","ntmix","ppRef",100,0,6,100,0,15)'
//   root -l -b -q 'plot2D_dataMC.C("Bmass","Bnorm_svpvDistance_2D","ntphi","ppRef",100,0,6,100,0,15)'


void plot2D_dataMC(TString TREE = "ntKstar", TString systemNAME = "ppRef")
{
    // Variables list and ranges (per variable) — order matters                                                 #{1.004,1.035}
    static const char* variables[]  = {  "BtrkminPt",  "Bmass",   "Btktkpt",    "Bujmass",   "Bpt",  "abs(By)",  "Balpha",     "BQvalue",      "Btktkmass", "Bcos_dtheta", "BtrkPtimb", "Bchi2cl", "Btrk1dR", "Btrk1Pt","Btrk2Pt", "Bnorm_svpvDistance_2D"};
    static const double ranges[][2] = {      {0.5,4},  {5 , 6},    {0,15},     {2.9,3.25}, {0, 50},   {0, 2.4}, {0, 3.15},     {0.5,1.5}, {0.75,1.05},          {0.999,1},                {0,0.9}, {0.05,1}, {0,1}, {0.5,5},{0.5,5}, {0,15}};
    int nVars = sizeof(variables)/sizeof(variables[0]);

    auto axisLabel = [&](const TString &var)->TString{
        if (var == "Bmass") {
            if (TREE == "ntKp")        return "m_{J/#Psi K^{+}} [GeV/c^{2}]";
            else if (TREE == "ntKstar")return "m_{J/#Psi K^{+} #pi^{-}} [GeV/c^{2}]";
            else if (TREE == "ntphi")  return "m_{J/#Psi K^{+} K^{-}} [GeV/c^{2}]";
            else if (TREE == "ntmix")  return "m_{J/#Psi #pi^{+} #pi^{-}} [GeV/c^{2}]";
        }
        if (var == "Bpt") return "p_{T} [GeV/c]";
        if (var == "abs(By)") return "|y|";
        return var;
    };

    TChain chain(TREE.Data());
    TTree *tree_MC = nullptr;
    TString path_to_MC;

    // Minimal file routing (extend as needed)
    if (systemNAME.Contains("PbPb")) {
        if (TREE == "ntmix") {
            chain.Add("/lstore/cms/hlegoinha/DATA_X3872_PbPb.root");
            path_to_MC = "/lstore/cms/lekai/X/MC/MC_X3872_pbpb.root";
        }
    } else {
        if (TREE == "ntmix") {
            chain.Add("/eos/user/k/kprince/x3872/DATA_pp_AANN.root");
            path_to_MC = "/eos/user/k/kprince/x3872/MC_X3872_pp_AANN.root";
        } else if (TREE == "ntphi") {
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_DATA.root");
            path_to_MC = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/MC_Bs.root";
        }
        else if (TREE == "ntKp") {
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_DATA.root");
            path_to_MC = "/eos/user/h/hmarques/Analysis_CODES/MC_Bu.root";
        }
        else if (TREE == "ntKstar") {
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_DATA.root");
            path_to_MC = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_MC.root";
        }
    }

    TFile *fmc = TFile::Open(path_to_MC.Data());
    fmc->GetObject(TREE.Data(), tree_MC); //TREE.Data()

    std::cout << "[plot2D] DATA entries: " << chain.GetEntries() << std::endl;
    std::cout << "[plot2D] MC   entries: " << tree_MC->GetEntries() << std::endl;

    TString ANYsel = "Bnorm_svpvDistance_2D > 4"; // customize if needed

    // Prepare output directory
    TString FOLDER_path = "./2Dplots";
    gSystem->mkdir(FOLDER_path, true);

    // Canvas setup
    TCanvas *c = new TCanvas("c2d", "", 800, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.15);
    gStyle->SetOptStat(0);

    // Also save a simple 1D plot of Bmass (DATA)
    {
        int nbinsMass = 80;
        double mmin = 0., mmax = 10.;
        if (TREE == "ntmix") { mmin = 3.6; mmax = 4.0; }
        else if (TREE == "ntphi" || TREE == "ntKp" || TREE == "ntKstar") { mmin = 5.0; mmax = 6.0; }

        TCanvas *c1d = new TCanvas("c1Dmass", "", 600, 500);
        TH1F *hmass_data = new TH1F("hmass_data", Form("; %s; Entries", axisLabel("Bmass").Data()), nbinsMass, mmin, mmax);
        chain.Draw("Bmass>>hmass_data", ANYsel, "goff");
        hmass_data->SetLineColor(kBlue);
        hmass_data->SetFillColor(kBlue);
        hmass_data->SetFillStyle(3358);
        hmass_data->Draw("HIST");
        c1d->SaveAs(Form("./2Dplots/%s_%s__Bmass__DATA_1D.pdf", TREE.Data(), systemNAME.Data()));
        delete hmass_data;
        delete c1d;
    }

    // Iterate all variable pairs (i != j) in given order
    for (int i = 0; i < nVars; ++i) {
        for (int j = 0; j < nVars; ++j) {
            if (i == j) continue;
            TString varX = variables[i];
            TString varY = variables[j];
            double xmin = ranges[i][0], xmax = ranges[i][1];
            double ymin = ranges[j][0], ymax = ranges[j][1];
            int nbinsX = 100, nbinsY = 100;

            if (TREE == "ntmix") {
                
                // MC signal: X(3872)
                TH2F *h2_sig = new TH2F(Form("h2_sig_%d_%d", i, j), "", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
                tree_MC->Draw(Form("%s:%s>>%s", varY.Data(), varX.Data(), h2_sig->GetName()), Form("%s && isX3872==1", ANYsel.Data()), "goff");
                h2_sig->GetXaxis()->SetTitle(axisLabel(varX));
                h2_sig->GetYaxis()->SetTitle(axisLabel(varY));
                h2_sig->Draw("COLZ");
                c->SaveAs(Form("./2Dplots/%s_%s__%s_VS_%s__MC_X3872.pdf", TREE.Data(), systemNAME.Data(), varY.Data(), varX.Data()));
                delete h2_sig;
                /*
                // MC specific background: Psi(2S)
                TH2F *h2_spec = new TH2F(Form("h2_spec_%d_%d", i, j), "", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
                tree_MC->Draw(Form("%s:%s>>%s", varY.Data(), varX.Data(), h2_spec->GetName()), Form("%s && isX3872==0", ANYsel.Data()), "goff");
                h2_spec->GetXaxis()->SetTitle(axisLabel(varX));
                h2_spec->GetYaxis()->SetTitle(axisLabel(varY));
                h2_spec->Draw("COLZ");
                c->SaveAs(Form("./2Dplots/%s_%s__%s_VS_%s__MC_Psi2S.pdf", TREE.Data(), systemNAME.Data(), varY.Data(), varX.Data()));
                delete h2_spec;
                
                // DATA (no sideband to keep simple)
                TH2F *h2_data = new TH2F(Form("h2_data_%d_%d", i, j), "", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
                chain.Draw(Form("%s:%s>>%s", varY.Data(), varX.Data(), h2_data->GetName()), ANYsel.Data(), "goff");
                h2_data->GetXaxis()->SetTitle(axisLabel(varX));
                h2_data->GetYaxis()->SetTitle(axisLabel(varY));
                h2_data->Draw("COLZ");
                c->SaveAs(Form("./2Dplots/%s_%s__%s_VS_%s__DATA.pdf", TREE.Data(), systemNAME.Data(), varY.Data(), varX.Data()));
                delete h2_data;
                */
            } else {

                // Generic MC 2D
                TH2F *h2_mc = new TH2F(Form("h2_mc_%d_%d", i, j), "", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
                tree_MC->Draw(Form("%s:%s>>%s", varY.Data(), varX.Data(), h2_mc->GetName()), ANYsel.Data(), "goff");
                h2_mc->GetXaxis()->SetTitle(axisLabel(varX));
                h2_mc->GetYaxis()->SetTitle(axisLabel(varY));
                h2_mc->Draw("COLZ");
                c->SaveAs(Form("./2Dplots/%s_%s__%s_VS_%s__MC.pdf", TREE.Data(), systemNAME.Data(), varY.Data(), varX.Data()));
                delete h2_mc;

                // DATA
                TH2F *h2_data = new TH2F(Form("h2_data_%d_%d", i, j), "", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
                chain.Draw(Form("%s:%s>>%s", varY.Data(), varX.Data(), h2_data->GetName()), ANYsel.Data(), "goff");
                h2_data->GetXaxis()->SetTitle(axisLabel(varX));
                h2_data->GetYaxis()->SetTitle(axisLabel(varY));
                h2_data->Draw("COLZ");
                c->SaveAs(Form("./2Dplots/%s/%s_%s__%s_VS_%s__DATA.pdf",TREE.Data(), TREE.Data(), systemNAME.Data(), varY.Data(), varX.Data()));
                delete h2_data;
            }
        }
    }


    delete c;
    fmc->Close();
}

int main() {
    // Default: iterate adjacent variable pairs using configured ranges
    plot2D_dataMC();
    return 0;
}
