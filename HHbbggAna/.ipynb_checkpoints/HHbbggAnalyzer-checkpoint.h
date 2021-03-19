#ifndef HHbbggAnalyzer_h
#define HHbbggAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <map>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "TRandom.h"
#include "MainEvent.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif
// Header file for the classes stored in the TTree if any.

class HHbbggAnalyzer : public MainEvent {
 public :
   HHbbggAnalyzer(const TString &inputFileList="foo.txt", const char *outFileName="histo.root", TString dataset="data",const char *isData="F", TString year_num="2017");
   virtual ~HHbbggAnalyzer();
   void Analyze(bool isData, int option, string outputFileName, string label);
   
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *);
  //declare any specific function required
  
  void clearTreeVectors();
  void BookTreeBranches();
  bool DataIs;
  TString year;
  std::string yearst;
  std::map<std::string,float> muon_pt_cut;
  std::map<std::string,float> btag_cut;

  TH1D *h_sumOfgw = new TH1D("h_sumOfgenWeight","h_sumOfgenWeight",1,0,1);
  TH1D *h_sumOfgpw = new TH1D("h_sumOfgenpuWeight","h_sumOfgenpuWeight",1,0,1);
    
  TFile *oFile;
  //TFile *ohistFile;
  TTree* tree;
  uint t_run;
  uint t_luminosityBlock;
  ulong t_event;
  float leading_photon_pt;
  float leading_photon_eta;
  float leading_photon_phi;
  float subleading_photon_pt;
  float subleading_photon_eta;
  float subleading_photon_phi;
  float diphoton_pt;
  float diphoton_eta;
  float diphoton_mass;
    
};

#endif

#ifdef HHbbggAnalyzer_cxx
HHbbggAnalyzer::HHbbggAnalyzer(const TString &inputFileList, const char *outFileName, TString dataset, const char *isData, TString year_num) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if(*isData!='T') DataIs = false;
  else DataIs = true;
  year = year_num;
  yearst = std::string(year.Data());

  h_sumOfgw->SetBinContent(1,0.0);
  h_sumOfgpw->SetBinContent(1,0.0);

  //muon pT selection
  muon_pt_cut["2016"] = 26.0;
  muon_pt_cut["2017"] = 29.0;
  muon_pt_cut["2018"] = 26.0;   
  //b-tag score selection
  btag_cut["2016"] = 0.6321; 
  btag_cut["2017"] = 0.4941;
  btag_cut["2018"] = 0.4184;

  TChain *tree = new TChain("Events");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  MainEvent::Init(tree);
  oFile = new TFile(outFileName, "recreate");
  TString histname(outFileName);
  //ohistFile = new TFile("hist_"+histname, "recreate");
  BookTreeBranches();
}

bool HHbbggAnalyzer::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;

  while ( getline (infile,buffer) )
  {
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
  }
  infile.close();

  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

 HHbbggAnalyzer::~HHbbggAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   h_sumOfgw->Write();
   h_sumOfgpw->Write();
   oFile->Write();
   oFile->Close();
}

Long64_t HHbbggAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HHbbggAnalyzer::clearTreeVectors(){
  t_run=0;
  t_luminosityBlock=0;
  t_event=0;
  leading_photon_pt=0;
  leading_photon_eta=0;
  leading_photon_phi=0;
  subleading_photon_pt=0;
  subleading_photon_eta=0;
  subleading_photon_phi=0;
  diphoton_pt=0;
  diphoton_eta=0;
  diphoton_mass=0;
}

void HHbbggAnalyzer::BookTreeBranches(){
  tree = new TTree("tree","tree");
  tree->SetAutoSave(10000);

  tree->Branch("run", &t_run,"run/i");
  tree->Branch("lumi", &t_luminosityBlock,"lumi/i");
  tree->Branch("event", &t_event,"event/l");
  
  //photon
  tree->Branch("leading_photon_pt", &leading_photon_pt,"leading_photon_pt/f");
  tree->Branch("leading_photon_eta", &leading_photon_eta,"leading_photon_eta/f");
  tree->Branch("leading_photon_phi", &leading_photon_phi,"leading_photon_phi/f");  
  tree->Branch("subleading_photon_pt", &subleading_photon_pt,"subleading_photon_pt/f");
  tree->Branch("subleading_photon_eta", &subleading_photon_eta,"subleading_photon_eta/f");
  tree->Branch("subleading_photon_phi", &subleading_photon_phi,"subleading_photon_phi/f"); 
  tree->Branch("diphoton_pt", &diphoton_pt,"diphoton_pt/f"); 
  tree->Branch("diphoton_mass", &diphoton_mass,"diphoton_mass/f"); 
  tree->Branch("diphoton_eta", &diphoton_eta,"diphoton_eta/f");
    
  //jets
}
#endif // #ifdef HHbbggAnalyzer_cxx
