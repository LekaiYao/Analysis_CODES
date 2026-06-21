#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "../plotER/aux/masses.h"

void Flat_TREEs( TString FILEIN="", TString NUN="", TString TREENAME="ntmix", TString SYSTEM="ppRef", TString KIND="MC", TString PARTICLE="", TString PVSNP="")
{
    bool isMC = (KIND == "MC");
    bool Fid_region  = true;       // apply fiducial region selection
    TString specCASES = PARTICLE + PVSNP;   // explicit case tag passed by the runner after filelist preselection

    std::cout << "Flattening tree: " << TREENAME << " , System: " << SYSTEM << " , kind: " << KIND << std::endl;
    bool isPP     = (SYSTEM == "ppRef") ;
    bool isPbPb   = (SYSTEM == "PbPb23" || SYSTEM == "PbPb24" || SYSTEM == "PbPb25");

    // --------------------------------------------------
    // Input files using TChain
    // --------------------------------------------------
    TChain *tin = new TChain("Bfinder/" + TREENAME);

    auto addFilesFromList = [&](const TString& listPath) {
        if (listPath == "") return;
        std::ifstream fin(listPath.Data());
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            tin->Add(line.c_str());
        }
    };

    addFilesFromList(FILEIN);
    std::cout << "Total files added: "      << tin->GetNtrees() << std::endl;
    std::cout << "Total entries in chain: " << tin->GetEntries() << std::endl;

    // --------------------------------------------------
    // Output file
    // --------------------------------------------------
    TString outputFile = "flat_" + TREENAME + "_" + SYSTEM + "_" + KIND + specCASES + NUN + ".root";  //
    TFile *fout = new TFile(outputFile, "RECREATE");
    TTree *tout = new TTree(TREENAME + PARTICLE, "Flattened tree");

    // --------------------------------------------------
    // Prepare gen-level input chain and output tree
    // --------------------------------------------------
    TChain *tgen = new TChain("Bfinder/ntGen");
    // Reuse exactly the same file list used by the reco chain
    TObjArray *fileList = tin->GetListOfFiles();
    if (fileList) {
        TIter next(fileList);
        TObject *obj = nullptr;
        while ((obj = next())) {
            const char *fname = obj->GetTitle();
            tgen->Add(fname);
        }
    }

    TTree *tgenOut = nullptr;

    // --------------------------------------------------
    // Input reconstructed branches 
    // --------------------------------------------------
    const Int_t MAXCAND = 40000;
    Int_t Bsize, nSelectedChargedTracks, CentBin;
    Float_t PVx, PVy, PVz, PVnchi2;
    static Float_t Bmass[MAXCAND], Bpt[MAXCAND], By[MAXCAND], Bgen[MAXCAND];
    static Float_t BvtxX[MAXCAND], BvtxY[MAXCAND];
    static Bool_t Bmu1isTriggered[MAXCAND];
    static Bool_t Bmu2isTriggered[MAXCAND];
    static Bool_t Bmu1SoftMuID[MAXCAND];
    static Bool_t Bmu2SoftMuID[MAXCAND];
    static Bool_t Bmu1HybridSoftMuID[MAXCAND];
    static Bool_t Bmu2HybridSoftMuID[MAXCAND];
    static Bool_t Bmu1isAcc[MAXCAND];
    static Bool_t Bmu2isAcc[MAXCAND];
    static Bool_t  Btrk1highPurity[MAXCAND];
    static Bool_t  Btrk2highPurity[MAXCAND];
    static Float_t Btrk1Pt[MAXCAND];
    static Float_t Btrk1Eta[MAXCAND];
    static Float_t Btrk1Phi[MAXCAND];
    static Float_t Btrk1PtErr[MAXCAND];
    static Float_t Btrk1Chi2ndf[MAXCAND];
    static Float_t Btrk2Pt[MAXCAND];
    static Float_t Btrk2Eta[MAXCAND];
    static Float_t Btrk2Phi[MAXCAND];
    static Float_t Btrk2PtErr[MAXCAND];
    static Float_t Btrk2Chi2ndf[MAXCAND];
    static Float_t Btktkmass[MAXCAND];
    static Float_t Bujmass[MAXCAND];
    static Float_t BujvProb[MAXCAND];
    static Float_t Bchi2Prob[MAXCAND];
    static Float_t BtrkPtimb[MAXCAND];
    static Float_t Btrk1dR[MAXCAND];
    static Float_t Btrk2dR[MAXCAND];
    static Float_t Btktkpt[MAXCAND];
    static Float_t BtktkvProb[MAXCAND];
    static Float_t BQvalue[MAXCAND];
    static Float_t Bnorm_svpvDistance_2D[MAXCAND];
    static Float_t BsvpvDistance_2D[MAXCAND];
    static Float_t BsvpvDisErr_2D[MAXCAND];
    static Float_t Balpha[MAXCAND];
    static Float_t Bdtheta[MAXCAND];
    static Float_t Bcos_dtheta[MAXCAND];
    static Float_t Bnorm_trk1Dxy[MAXCAND];
    static Float_t Bnorm_trk2Dxy[MAXCAND];
    static Float_t Btrk1nPixelLayer[MAXCAND];
    static Float_t Btrk1nStripLayer[MAXCAND];
    static Float_t Btrk2nPixelLayer[MAXCAND];  
    static Float_t Btrk2nStripLayer[MAXCAND];
    static Float_t BLxy[MAXCAND];
    static Float_t Bmu1y[MAXCAND];
    static Float_t Bmu2y[MAXCAND];
    static Float_t Bmu1pt[MAXCAND];
    static Float_t Bmu2pt[MAXCAND];
    static Float_t Bnorm_trk1Dz[MAXCAND];
    static Float_t Bnorm_trk2Dz[MAXCAND];


    tin->SetBranchAddress("Bsize", &Bsize);
    tin->SetBranchAddress("PVx", &PVx);
    tin->SetBranchAddress("PVy", &PVy);
    tin->SetBranchAddress("PVz", &PVz);
    tin->SetBranchAddress("PVnchi2", &PVnchi2);
    tin->SetBranchAddress("nSelectedChargedTracks", &nSelectedChargedTracks);
    tin->SetBranchAddress("CentBin", &CentBin);
    tin->SetBranchAddress("Bmass", Bmass);
    tin->SetBranchAddress("BvtxX", BvtxX);
    tin->SetBranchAddress("BvtxY", BvtxY);
    tin->SetBranchAddress("Bpt", Bpt);
    tin->SetBranchAddress("By", By);
    tin->SetBranchAddress("Bgen", Bgen);
    tin->SetBranchAddress("Bmu1isTriggered", Bmu1isTriggered);
    tin->SetBranchAddress("Bmu2isTriggered", Bmu2isTriggered);
    tin->SetBranchAddress("Bmu1SoftMuID", Bmu1SoftMuID);
    tin->SetBranchAddress("Bmu2SoftMuID", Bmu2SoftMuID);
    tin->SetBranchAddress("Bmu1HybridSoftMuID", Bmu1HybridSoftMuID);
    tin->SetBranchAddress("Bmu2HybridSoftMuID", Bmu2HybridSoftMuID);
    tin->SetBranchAddress("Bmu1y", Bmu1y);
    tin->SetBranchAddress("Bmu2y", Bmu2y);
    tin->SetBranchAddress("Bmu1isAcc", Bmu1isAcc);
    tin->SetBranchAddress("Bmu2isAcc", Bmu2isAcc);  
    tin->SetBranchAddress("Btrk1highPurity", Btrk1highPurity);
    tin->SetBranchAddress("Btrk2highPurity", Btrk2highPurity);
    tin->SetBranchAddress("Btrk1Pt", Btrk1Pt);
    tin->SetBranchAddress("Btrk1Eta", Btrk1Eta);
    tin->SetBranchAddress("Btrk1Phi", Btrk1Phi);
    tin->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
    tin->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf);
    tin->SetBranchAddress("Btrk2Pt", Btrk2Pt);
    tin->SetBranchAddress("Btrk2Eta", Btrk2Eta);
    tin->SetBranchAddress("Btrk2Phi", Btrk2Phi);
    tin->SetBranchAddress("Btrk2PtErr", Btrk2PtErr);
    tin->SetBranchAddress("Btrk2Chi2ndf", Btrk2Chi2ndf);
    tin->SetBranchAddress("Btktkmass", Btktkmass);
    tin->SetBranchAddress("Bujmass", Bujmass);
    tin->SetBranchAddress("BujvProb", BujvProb);
    tin->SetBranchAddress("Bchi2Prob", Bchi2Prob);
    tin->SetBranchAddress("BtrkPtimb", BtrkPtimb);
    tin->SetBranchAddress("Btrk1dR", Btrk1dR);
    tin->SetBranchAddress("Btrk2dR", Btrk2dR);
    tin->SetBranchAddress("Btktkpt", Btktkpt);
    tin->SetBranchAddress("BtktkvProb", BtktkvProb);
    tin->SetBranchAddress("BQvalue", BQvalue);
    tin->SetBranchAddress("Bnorm_svpvDistance_2D", Bnorm_svpvDistance_2D);
    tin->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D);
    tin->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D);
    tin->SetBranchAddress("Balpha", Balpha);
    tin->SetBranchAddress("Bdtheta", Bdtheta);
    tin->SetBranchAddress("Bcos_dtheta", Bcos_dtheta);
    tin->SetBranchAddress("Bnorm_trk1Dxy", Bnorm_trk1Dxy);
    tin->SetBranchAddress("Bnorm_trk2Dxy", Bnorm_trk2Dxy);
    tin->SetBranchAddress("Btrk1nPixelLayer", Btrk1nPixelLayer);
    tin->SetBranchAddress("Btrk1nStripLayer", Btrk1nStripLayer);
    tin->SetBranchAddress("Btrk2nPixelLayer", Btrk2nPixelLayer);  
    tin->SetBranchAddress("Btrk2nStripLayer", Btrk2nStripLayer);
    tin->SetBranchAddress("BLxy", BLxy);
    tin->SetBranchAddress("Bmu1pt", Bmu1pt);
    tin->SetBranchAddress("Bmu2pt", Bmu2pt);
    tin->SetBranchAddress("Bnorm_trk1Dz", Bnorm_trk1Dz);
    tin->SetBranchAddress("Bnorm_trk2Dz", Bnorm_trk2Dz);
    // --------------------------------------------------
    // Output branches (for Selection // Analysis)
    // --------------------------------------------------
    float PVx_out, PVy_out, PVz_out, PVnchi2_out;
    float Bmass_out, Bpt_out, By_out, Bchi2Prob_out, Btrk1dR_out, Btrk2dR_out, BtrkPtimb_out, Btktkpt_out, Bujmass_out, Bnorm_svpvDistance_2D_out;
    float BQvalue_out, Bnorm_trk1Dxy_out, Bnorm_trk2Dxy_out, Balpha_out, Bdtheta_out, Bcos_dtheta_out, Btktkmass_out, Btrk1Pt_out, Btrk2Pt_out;
    float Bgen_out, BsvpvDistance_2D_out, BtktkvProb_out, BLxy_out;
    float BvtxX_out, BvtxY_out, BsvpvDisErr_2D_out;
    float Btrk1Eta_out, Btrk2Eta_out, Btrk1Phi_out, Btrk2Phi_out, Btrk1PtErr_out, Btrk2PtErr_out;
    float BujvProb_out, Bmu1y_out, Bmu2y_out, Bmu1pt_out, Bmu2pt_out, Bnorm_trk1Dz_out, Bnorm_trk2Dz_out;
    int CentBin_out, nSelectedChargedTracks_out;

    tout->Branch("PVx", &PVx_out, "PVx/F");
    tout->Branch("PVy", &PVy_out, "PVy/F");
    tout->Branch("PVz", &PVz_out, "PVz/F");
    tout->Branch("PVnchi2", &PVnchi2_out, "PVnchi2/F");
    tout->Branch("nSelectedChargedTracks", &nSelectedChargedTracks_out, "nSelectedChargedTracks/I");
    tout->Branch("CentBin", &CentBin_out, "CentBin/I");
    tout->Branch("Bmass", &Bmass_out, "Bmass/F");
    tout->Branch("Bpt", &Bpt_out, "Bpt/F");
    tout->Branch("By", &By_out, "By/F");
    tout->Branch("Bchi2Prob", &Bchi2Prob_out, "Bchi2Prob/F");
    tout->Branch("Btrk1dR", &Btrk1dR_out, "Btrk1dR/F");
    tout->Branch("Btrk2dR", &Btrk2dR_out, "Btrk2dR/F");
    tout->Branch("BtrkPtimb", &BtrkPtimb_out, "BtrkPtimb/F");
    tout->Branch("Btktkpt", &Btktkpt_out, "Btktkpt/F");
    tout->Branch("Bujmass", &Bujmass_out, "Bujmass/F");
    tout->Branch("BujvProb", &BujvProb_out, "BujvProb/F");
    tout->Branch("Bnorm_svpvDistance_2D", &Bnorm_svpvDistance_2D_out, "Bnorm_svpvDistance_2D/F");
    tout->Branch("BsvpvDistance_2D", &BsvpvDistance_2D_out, "BsvpvDistance_2D/F");
    tout->Branch("BsvpvDisErr_2D", &BsvpvDisErr_2D_out, "BsvpvDisErr_2D/F");
    tout->Branch("BQvalue", &BQvalue_out, "BQvalue/F");
    tout->Branch("Bnorm_trk1Dxy", &Bnorm_trk1Dxy_out, "Bnorm_trk1Dxy/F");
    tout->Branch("Bnorm_trk2Dxy", &Bnorm_trk2Dxy_out, "Bnorm_trk2Dxy/F");
    tout->Branch("Balpha", &Balpha_out, "Balpha/F");
    tout->Branch("Bdtheta", &Bdtheta_out, "Bdtheta/F");
    tout->Branch("Bcos_dtheta", &Bcos_dtheta_out, "Bcos_dtheta/F");
    tout->Branch("Btktkmass", &Btktkmass_out, "Btktkmass/F");
    tout->Branch("Btrk1Pt", &Btrk1Pt_out, "Btrk1Pt/F");
    tout->Branch("Btrk2Pt", &Btrk2Pt_out, "Btrk2Pt/F");
    tout->Branch("Btrk1Eta", &Btrk1Eta_out, "Btrk1Eta/F");
    tout->Branch("Btrk2Eta", &Btrk2Eta_out, "Btrk2Eta/F");
    tout->Branch("Btrk1Phi", &Btrk1Phi_out, "Btrk1Phi/F");
    tout->Branch("Btrk2Phi", &Btrk2Phi_out, "Btrk2Phi/F");
    tout->Branch("Btrk1PtErr", &Btrk1PtErr_out, "Btrk1PtErr/F");
    tout->Branch("Btrk2PtErr", &Btrk2PtErr_out, "Btrk2PtErr/F");
    tout->Branch("BtktkvProb", &BtktkvProb_out, "BtktkvProb/F");
    tout->Branch("BLxy", &BLxy_out, "BLxy/F");
    tout->Branch("BvtxX", &BvtxX_out, "BvtxX/F");
    tout->Branch("BvtxY", &BvtxY_out, "BvtxY/F");
    tout->Branch("Bmu1y", &Bmu1y_out, "Bmu1y/F");
    tout->Branch("Bmu2y", &Bmu2y_out, "Bmu2y/F");
    tout->Branch("Bmu1pt", &Bmu1pt_out, "Bmu1pt/F");
    tout->Branch("Bmu2pt", &Bmu2pt_out, "Bmu2pt/F");
    tout->Branch("Bnorm_trk1Dz", &Bnorm_trk1Dz_out, "Bnorm_trk1Dz/F");
    tout->Branch("Bnorm_trk2Dz", &Bnorm_trk2Dz_out, "Bnorm_trk2Dz/F");
    if (isMC){
        tout->Branch("Bgen", &Bgen_out, "Bgen/F");
    }

    // --------------------------------------------------
    // Event loop
    // --------------------------------------------------
    Long64_t nentries = tin->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;
    for(Long64_t ev=0; ev<nentries; ++ev)
    {
        // keep track of progress
        if (ev % 100000 == 0) {
            std::cout << "Processing event " << ev << " / " << nentries 
                      << " (" << (100.0*ev/nentries) << "%)" << std::endl;
        }
        tin->GetEntry(ev);

        // Event-level variables
        PVx_out = PVx;
        PVy_out = PVy;
        PVz_out = PVz;
        PVnchi2_out = PVnchi2;
        CentBin_out = CentBin;
        nSelectedChargedTracks_out = nSelectedChargedTracks;

        // --------------------------------------------------
        // Candidate loop
        // --------------------------------------------------
        for(int i=0; i<Bsize; ++i)
        {
            // Pre-Selection (quality) Cuts
            // DATA is already filtered (except for the cuts marked with <-------- ** ); 
            // MC is not.
            if (true){
                // Muons
                // Muon Trigger matching   <-------- **
                if(isPP && !(Bmu1isTriggered[i] && Bmu2isTriggered[i] )) continue;
                // Soft vs HybridSoft Muon <-------- **
                if((isPP   && !(Bmu1SoftMuID[i] && Bmu2SoftMuID[i] )) ||
                   (isPbPb && !(Bmu1HybridSoftMuID[i] && Bmu2HybridSoftMuID[i]))) continue;
                // Muon Acceptance
                if( !(Bmu1isAcc[i] && Bmu2isAcc[i])) continue;
                // diMuon system 
                if( !(abs(Bujmass[i] - JPSI_MASS) < 0.15)) continue;
                if( !(BujvProb[i] > 0.01)) continue;

                // Tracks
                // single Track channel
                if(abs(Btrk1Eta[i]) > 2.4) continue;
                if((isPP && !(Btrk1Pt[i] > 0.5)) || (isPbPb && !(Btrk1Pt[i] > 0.9)) ) continue;
                if(((Btrk1PtErr[i] / Btrk1Pt[i]) >= 0.1)) continue;
                if(((Btrk1nPixelLayer[i] + Btrk1nStripLayer[i]) <= 10)) continue;
                if(((Btrk1Chi2ndf[i]/(Btrk1nPixelLayer[i] + Btrk1nStripLayer[i])) >= 0.18)) continue;
                if(!Btrk1highPurity[i]) continue;
                // diTrack channel 
                if (TREENAME !=  "ntKp"){ 
                    if(abs(Btrk2Eta[i]) > 2.4) continue;
                    if((isPP && !(Btrk2Pt[i] > 0.5)) || (isPbPb && !(Btrk2Pt[i] > 0.9)) ) continue;
                    if(((Btrk2PtErr[i] / Btrk2Pt[i]) >= 0.1)) continue;
                    if(((Btrk2nPixelLayer[i] + Btrk2nStripLayer[i]) <= 10)) continue;
                    if(((Btrk2Chi2ndf[i]/(Btrk2nPixelLayer[i] + Btrk2nStripLayer[i])) >= 0.18)) continue;
                    if(!Btrk2highPurity[i]) continue;
                    //diTrack system
                    if (TREENAME == "ntphi"   && !(abs(Btktkmass[i] - PHI_MASS)  <0.015)) continue;
                    if (TREENAME == "ntKstar" && !(abs(Btktkmass[i] - KSTAR_MASS)<0.150)) continue;
                }

                // Candidate
                if (TREENAME == "ntmix" && (Bmass[i] > 4.0 || Bmass[i] < 3.6)) continue;
                if(Bchi2Prob[i] < 0.005) continue;
                if(TREENAME == "ntmix" && Bpt[i]<4) continue; 
                else if (Bpt[i]<1) continue;
            }

            // Fiducial Region
            if (Fid_region){
                if (!((Bpt[i] >= 5) && abs(By[i]) <= 2.4)) continue;
            }

            // keep MC signal only
            Bgen_out = Bgen[i];
            if( isMC && !(Bgen_out==23333 || Bgen_out==24333 || Bgen_out==23433 || Bgen_out==24433 || Bgen_out==41000 )) continue;

            Bmass_out   = Bmass[i];
            Bpt_out     = Bpt[i];
            By_out      = By[i];
            Bchi2Prob_out = Bchi2Prob[i];
            Btrk1dR_out = Btrk1dR[i];
            Btrk2dR_out = Btrk2dR[i];
            BtrkPtimb_out = BtrkPtimb[i];
            Btktkpt_out   = Btktkpt[i];
            Bujmass_out   = Bujmass[i];
            BujvProb_out  = BujvProb[i];
            Bnorm_svpvDistance_2D_out = Bnorm_svpvDistance_2D[i];
            BsvpvDistance_2D_out = BsvpvDistance_2D[i];
            BsvpvDisErr_2D_out = BsvpvDisErr_2D[i];
            BQvalue_out = BQvalue[i];
            Bnorm_trk1Dxy_out = Bnorm_trk1Dxy[i];
            Bnorm_trk2Dxy_out = Bnorm_trk2Dxy[i];
            Balpha_out  = Balpha[i];
            BtktkvProb_out = BtktkvProb[i];
            Bdtheta_out = Bdtheta[i];
            Bcos_dtheta_out = Bcos_dtheta[i];
            Btktkmass_out = Btktkmass[i];
            Btrk1Pt_out   = Btrk1Pt[i];
            Btrk2Pt_out   = Btrk2Pt[i];
            Btrk1Eta_out = Btrk1Eta[i];
            Btrk2Eta_out = Btrk2Eta[i];
            Btrk1Phi_out = Btrk1Phi[i];
            Btrk2Phi_out = Btrk2Phi[i];
            Btrk1PtErr_out = Btrk1PtErr[i];
            Btrk2PtErr_out = Btrk2PtErr[i];
            BvtxX_out = BvtxX[i];
            BvtxY_out = BvtxY[i];
            Bmu1y_out = Bmu1y[i];
            Bmu2y_out = Bmu2y[i];
            BLxy_out  = BLxy[i];
            Bmu1pt_out = Bmu1pt[i];
            Bmu2pt_out = Bmu2pt[i];
            Bnorm_trk1Dz_out = Bnorm_trk1Dz[i];
            Bnorm_trk2Dz_out = Bnorm_trk2Dz[i];

            if(!std::isfinite(Bmass_out) || !std::isfinite(Bpt_out) || !std::isfinite(By_out) || !std::isfinite(Bnorm_trk1Dxy_out) ||
            !std::isfinite(CentBin_out) || !std::isfinite(Bchi2Prob_out) || !std::isfinite(Btrk1dR_out) || !std::isfinite(Bnorm_svpvDistance_2D_out)) continue;
            tout->Fill();
        }
    }

    // --------------------------------------------------
    // Copy ntGen tree (flat, one row per candidate in mass window)
    // --------------------------------------------------
    if (isMC && true) {
        fout->cd();

        // Gen-level input branches
        Int_t Gsize;
        static Int_t GpdgId[MAXCAND];
        static Float_t Gmu1eta[MAXCAND];
        static Float_t Gmu1pt[MAXCAND];
        static Float_t Gmu2eta[MAXCAND];
        static Float_t Gmu2pt[MAXCAND];
        static Float_t Gtk1pt[MAXCAND];
        static Float_t Gtk1eta[MAXCAND];
        static Float_t Gtk2pt[MAXCAND];
        static Float_t Gtk2eta[MAXCAND];
        static Float_t Gpt[MAXCAND];
        static Float_t Gy[MAXCAND];
        static Int_t GisSignal[MAXCAND];

        tgen->SetBranchAddress("Gsize", &Gsize);
        tgen->SetBranchAddress("GisSignal", GisSignal);
        tgen->SetBranchAddress("Gmu1eta", Gmu1eta);
        tgen->SetBranchAddress("Gmu1pt", Gmu1pt);
        tgen->SetBranchAddress("Gmu2eta", Gmu2eta);
        tgen->SetBranchAddress("Gmu2pt", Gmu2pt);
        tgen->SetBranchAddress("Gtk1pt", Gtk1pt);
        tgen->SetBranchAddress("Gtk1eta", Gtk1eta);
        tgen->SetBranchAddress("Gtk2pt", Gtk2pt);
        tgen->SetBranchAddress("Gtk2eta", Gtk2eta);
        tgen->SetBranchAddress("GpdgId", GpdgId);
        tgen->SetBranchAddress("Gpt", Gpt);
        tgen->SetBranchAddress("Gy", Gy);

        // Output tree (flat)
        tgenOut = new TTree("ntGen", "Gen-level (flat, filtered)");
        Float_t Gmu1eta_out, Gmu1pt_out, Gmu2eta_out, Gmu2pt_out;
        Float_t Gtk1pt_out, Gtk1eta_out, Gtk2pt_out, Gtk2eta_out, Gpt_out, Gy_out;
        Int_t GpdgId_out;
        tgenOut->Branch("Gmu1eta", &Gmu1eta_out, "Gmu1eta/F");
        tgenOut->Branch("Gmu1pt", &Gmu1pt_out, "Gmu1pt/F");
        tgenOut->Branch("Gmu2eta", &Gmu2eta_out, "Gmu2eta/F");
        tgenOut->Branch("Gmu2pt", &Gmu2pt_out, "Gmu2pt/F");
        tgenOut->Branch("Gtk1pt", &Gtk1pt_out, "Gtk1pt/F");
        tgenOut->Branch("Gtk1eta", &Gtk1eta_out, "Gtk1eta/F");
        tgenOut->Branch("Gtk2pt", &Gtk2pt_out, "Gtk2pt/F");
        tgenOut->Branch("Gtk2eta", &Gtk2eta_out, "Gtk2eta/F");
        tgenOut->Branch("GpdgId", &GpdgId_out, "GpdgId/I");
        tgenOut->Branch("Gpt", &Gpt_out, "Gpt/F");
        tgenOut->Branch("Gy", &Gy_out, "Gy/F");

        const Long64_t ngen = tgen->GetEntries();
        for (Long64_t ev=0; ev<ngen; ++ev) {
            tgen->GetEntry(ev);
            for (int i=0; i<Gsize; ++i) {

                if (TREENAME == "ntmix" && GisSignal[i] != 7) continue;               // keep only signal ntmix candidates for the selected MC sample
                else if (TREENAME == "ntphi"   && abs(GpdgId[i]) != 531) continue;  
                else if (TREENAME == "ntKstar" && abs(GpdgId[i]) != 511) continue;  
                else if (TREENAME == "ntKp"    && abs(GpdgId[i]) != 521) continue;

                // Fiducial Region
                if (Fid_region){
                    if ( !( (Gpt[i] >= 5) && abs(Gy[i]) <= 2.4) ) continue; 
                }

                GpdgId_out  = GpdgId[i];
                Gmu1eta_out = Gmu1eta[i];
                Gmu1pt_out  = Gmu1pt[i];
                Gmu2eta_out = Gmu2eta[i];
                Gmu2pt_out  = Gmu2pt[i];
                Gtk1pt_out  = Gtk1pt[i];
                Gtk1eta_out = Gtk1eta[i];
                Gtk2pt_out  = Gtk2pt[i];
                Gtk2eta_out = Gtk2eta[i];
                Gpt_out = Gpt[i];
                Gy_out  = Gy[i];
                tgenOut->Fill();

            }
        }
    }

    std::cout << "Output tree has " << tout->GetEntries() << " entries" << std::endl;
    // Write trees (only two top-level trees)
    tout->Write();
    if (isMC && tgenOut) tgenOut->Write();
    fout->Close();
    delete tin;
    std::cout << "Done & Saved -> " << outputFile << "\n";
}
