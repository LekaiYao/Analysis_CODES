#include <TFile.h>
#include <TList.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <memory>

void PrepareXMergedFitCache(const char* pb23Path,const char* pb23Tree,const char* pb23Selection,
 const char* pb24Path,const char* pb24Tree,const char* pb24Selection,
 const char* outputPath,const char* outputTree,const char* metadataPath) {
    TFile f23(pb23Path,"READ"),f24(pb24Path,"READ");
    auto* t23=dynamic_cast<TTree*>(f23.Get(pb23Tree));
    auto* t24=dynamic_cast<TTree*>(f24.Get(pb24Tree));
    if(!t23||!t24){std::cerr<<"missing input cache tree\n";gSystem->Exit(1);return;}
    gROOT->cd();
    std::unique_ptr<TTree> s23(t23->CopyTree(pb23Selection));
    std::unique_ptr<TTree> s24(t24->CopyTree(pb24Selection));
    if(!s23||!s24){std::cerr<<"selection failed\n";gSystem->Exit(2);return;}
    const Long64_t n23=s23->GetEntries(),n24=s24->GetEntries();
    TList trees; trees.Add(s23.get()); trees.Add(s24.get());
    std::unique_ptr<TTree> merged(TTree::MergeTrees(&trees));
    if(!merged){std::cerr<<"tree merge failed\n";gSystem->Exit(3);return;}
    merged->SetName(outputTree); merged->SetTitle("PbPb23+PbPb24 xeff25 merged fit cache");
    TFile output(outputPath,"RECREATE"); merged->Write(); output.Close();
    std::ofstream metadata(metadataPath);
    metadata<<"{\n  \"pb23_entries\": "<<n23<<",\n  \"pb24_entries\": "<<n24
            <<",\n  \"merged_entries\": "<<merged->GetEntries()<<"\n}\n";
    std::cout<<"[merged cache] pb23="<<n23<<" pb24="<<n24
             <<" merged="<<merged->GetEntries()<<std::endl;
}
