#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

void plot2D_dataMC(TString TREE = "ntmix_PSI2S", TString systemNAME = "ppRef", TString cut = "1")
{
    TString dataPath;
    TString mcPath;
    TString dataTree;
    TString mcTree;
    TString outTag = TREE;
    TString mcLabel = "MC";
    const bool isNtmixSignal = (TREE == "ntmix_X3872" || TREE == "ntmix_PSI2S");

    if (TREE == "ntmix_X3872") {
        dataTree = "ntmix";
        mcTree = "ntmix_X3872";
        mcLabel = "X(3872) MC";
        if (systemNAME == "ppRef") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root";
        }
        if (systemNAME == "PbPb23") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root";
        }
        if (systemNAME == "PbPb24") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872.root";
        }
    }

    if (TREE == "ntmix_PSI2S") {
        dataTree = "ntmix";
        mcTree = "ntmix_PSI2S";
        mcLabel = "#psi(2S) MC";
        if (systemNAME == "ppRef") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root";
        }
        if (systemNAME == "PbPb23") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root";
        }
        if (systemNAME == "PbPb24") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root";
        }
    }

    if (TREE == "ntphi") {
        dataTree = "ntphi";
        mcTree = "ntphi";
        if (systemNAME == "ppRef") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/Data_Bs.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/MC_Bs.root";
        }
        if (systemNAME == "PbPb23") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntphi_PbPb23_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntphi_PbPb23_MC.root";
        }
        if (systemNAME == "PbPb24") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntphi_PbPb24_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntphi_PbPb24_MC.root";
        }
    }

    if (TREE == "ntKp") {
        dataTree = "ntKp";
        mcTree = "ntKp";
        if (systemNAME == "ppRef") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_MC.root";
        }
        if (systemNAME == "PbPb23") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntKp_PbPb23_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntKp_PbPb23_MC.root";
        }
        if (systemNAME == "PbPb24") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntKp_PbPb24_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntKp_PbPb24_MC.root";
        }
    }

    if (TREE == "ntKstar") {
        dataTree = "ntKstar";
        mcTree = "ntKstar";
        if (systemNAME == "ppRef") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKstar_ppRef_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKstar_ppRef_MC.root";
        }
        if (systemNAME == "PbPb23") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntKstar_PbPb23_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_ntKstar_PbPb23_MC.root";
        }
        if (systemNAME == "PbPb24") {
            dataPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntKstar_PbPb24_DATA.root";
            mcPath = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb24/flat_ntKstar_PbPb24_MC.root";
        }
    }

    const TString outDir = Form("2Dplots/%s_%s", outTag.Data(), systemNAME.Data());
    gSystem->mkdir(outDir, true);

    std::cout << "[plot2D_dataMC] TREE   = " << TREE << std::endl;
    std::cout << "[plot2D_dataMC] SYSTEM = " << systemNAME << std::endl;
    std::cout << "[plot2D_dataMC] CUT    = " << cut << std::endl;
    std::cout << "[plot2D_dataMC] DATA   = " << dataPath << " (" << dataTree << ")" << std::endl;
    std::cout << "[plot2D_dataMC] MC     = " << mcPath << " (" << mcTree << ")" << std::endl;

    TFile* dataFile = TFile::Open(dataPath, "READ");
    TFile* mcFile = TFile::Open(mcPath, "READ");
    TTree* treeData = nullptr;
    TTree* treeMC = nullptr;
    dataFile->GetObject(dataTree, treeData);
    mcFile->GetObject(mcTree, treeMC);

    const int nMaps = 2;
    TString xVar[nMaps] = {"Bpt", "Bmass"};
    TString yVar[nMaps] = {"abs(By)", "Btktkmass"};
    TString mapTag[nMaps] = {"absY_vs_Bpt", "Btktkmass_vs_Bmass"};
    TString xLabel[nMaps] = {"p_{T} [GeV/c]", "B mass [GeV/c^{2}]"};
    TString yLabel[nMaps] = {"|y|", "m_{trk trk} [GeV/c^{2}]"};
    int xBins[nMaps] = {60, 80};
    int yBins[nMaps] = {48, 80};
    double xMin[nMaps] = {0.0, isNtmixSignal ? 3.6 : 5.0};
    double xMax[nMaps] = {60.0, isNtmixSignal ? 4.0 : 6.0};
    double yMin[nMaps] = {0.0, 0.0};
    double yMax[nMaps] = {2.4, 1.2};
    if (TREE == "ntmix_PSI2S") {
        xMin[1] = 3.6;
        xMax[1] = 3.8;
    } else if (TREE == "ntmix_X3872") {
        xMin[1] = 3.8;
        xMax[1] = 4.0;
    }
    if (TREE == "ntphi") {
        yMin[1] = 0.95;
        yMax[1] = 1.10;
    }
    if (TREE == "ntKstar") {
        yMin[1] = 0.70;
        yMax[1] = 1.20;
    }
    if (isNtmixSignal) {
        yMin[1] = 0.2;
        yMax[1] = 1.2;
    }

    TString massWindowCut = "";
    if (TREE == "ntmix_PSI2S") {
        massWindowCut = "(Bmass >= 3.6 && Bmass <= 3.8)";
    } else if (TREE == "ntmix_X3872") {
        massWindowCut = "(Bmass >= 3.8 && Bmass <= 4.0)";
    }
    TString drawCut = cut;
    if (!massWindowCut.IsNull()) {
        drawCut = Form("(%s) && %s", cut.Data(), massWindowCut.Data());
    }

    TCanvas canvas("canvas_2d", "", 850, 720);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.16);
    canvas.SetBottomMargin(0.12);
    gStyle->SetOptStat(0);

    for (int i = 0; i < nMaps; ++i) {
        TH2F hData(Form("hData_%d", i), Form(";%s;%s;Entries", xLabel[i].Data(), yLabel[i].Data()), xBins[i], xMin[i], xMax[i], yBins[i], yMin[i], yMax[i]);
        treeData->Draw(Form("%s:%s>>%s", yVar[i].Data(), xVar[i].Data(), hData.GetName()), drawCut, "goff");
        hData.SetDirectory(nullptr);
        hData.Draw("COLZ");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.15, 0.91, Form("%s DATA", systemNAME.Data()));
        canvas.SaveAs(Form("%s/%s_%s_DATA_%s.pdf", outDir.Data(), outTag.Data(), systemNAME.Data(), mapTag[i].Data()));

        TH2F hMC(Form("hMC_%d", i), Form(";%s;%s;Entries", xLabel[i].Data(), yLabel[i].Data()), xBins[i], xMin[i], xMax[i], yBins[i], yMin[i], yMax[i]);
        treeMC->Draw(Form("%s:%s>>%s", yVar[i].Data(), xVar[i].Data(), hMC.GetName()), drawCut, "goff");
        hMC.SetDirectory(nullptr);
        hMC.Draw("COLZ");
        latex.DrawLatex(0.15, 0.91, Form("%s %s", systemNAME.Data(), mcLabel.Data()));
        canvas.SaveAs(Form("%s/%s_%s_MC_%s.pdf", outDir.Data(), outTag.Data(), systemNAME.Data(), mapTag[i].Data()));
    }

    dataFile->Close();
    mcFile->Close();
}

int main()
{
    plot2D_dataMC();
    return 0;
}
