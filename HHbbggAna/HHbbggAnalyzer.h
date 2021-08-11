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
   HHbbggAnalyzer(const TString &inputFileList="foo.txt", const char *outFileName="histo.root", TString dataset="data", const char *isData="F", TString year_num="2017");
   virtual ~HHbbggAnalyzer();
   void Analyze(bool isData, int option, string outputFileName, string label);
   
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(string , const char *, const char *);
  void     cal_sumOfgw(string , const char *);
  //declare any specific function required
  
  void clearTreeVectors();
  void BookTreeBranches();
  bool DataIs;
  TString year;
  std::string yearst;
  std::map<std::string,float> muon_pt_cut;
  std::map<std::string,float> btag_cut;
  std::map<std::string, float> luminosity;
  std::map<std::string, float> xs;
  std::map<std::string, std::map<std::string, float>> sumOfgenw;
  
  TH1D *h_sumOfgw = new TH1D("h_sumOfgenWeight","h_sumOfgenWeight",1,0,1);
  TH1D *h_sumOfgpw = new TH1D("h_sumOfgenpuWeight","h_sumOfgenpuWeight",1,0,1);
    
  TFile *oFile;
  //TFile *ohistFile;
  TTree* tree;
  uint t_run;
  uint t_luminosityBlock;
  ulong t_event;
    
  int trig_decision;
  int pv_pass;
  float leading_photon_pt;
  float leading_photon_eta;
  float leading_photon_phi;
  float subleading_photon_pt;
  float subleading_photon_eta;
  float subleading_photon_phi;
  float diphoton_pt;
  float diphoton_eta;
  float diphoton_mass;
  float photon_delR;
  // added for bjet reconstruction
  float leading_bjet_pt;
  float leading_bjet_eta;
  float leading_bjet_phi;
  float subleading_bjet_pt;
  float subleading_bjet_eta;
  float subleading_bjet_phi;
  float dibjet_pt;
  float dibjet_eta;
  float dibjet_mass;
  float bjet_delR;
  int nbjet;
    
  // vbf jet reconstruction
  float leading_vbfjet_pt;
  float leading_vbfjet_eta;
  float leading_vbfjet_phi;
  float subleading_vbfjet_pt;
  float subleading_vbfjet_eta;
  float subleading_vbfjet_phi;
  float divbfjet_pt;
  float divbfjet_eta;
  float divbfjet_mass;
  float vbfjet_delR;
  float vbfjet_del_eta;
  int nvbfjet;
  
    
  // bjet corrections
  float leading_bjet_pt_corr;
  float subleading_bjet_pt_corr;
  float dibjet_pt_corr;
  float dibjet_eta_corr;
  float dibjet_mass_corr;
    
  //gen information
  float genHH_pt;
  float genHH_eta;
  float genHH_phi;
  float genHH_mass;
  float genLeadingH_pt;
  float genLeadingH_eta;
  float genLeadingH_phi;
  float genLeadingH_mass;
  float gensubLeadingH_pt;
  float gensubLeadingH_eta;
  float gensubLeadingH_phi;
  float gensubLeadingH_mass;
  float genLeadingPho_pt;
  float genLeadingPho_eta;
  float genLeadingPho_phi;
  float genLeadingPho_mass;
  float gensubLeadingPho_pt;
  float gensubLeadingPho_eta;
  float gensubLeadingPho_phi;
  float gensubLeadingPho_mass;
    
  float gen_matched_LeadingPho_pt;
  float gen_matched_LeadingPho_eta;
  float gen_matched_LeadingPho_phi;
  float gen_matched_subLeadingPho_pt;
  float gen_matched_subLeadingPho_eta;
  float gen_matched_subLeadingPho_phi;
    
  float genLeadingBjet_pt;
  float genLeadingBjet_eta;
  float genLeadingBjet_phi;
  float genLeadingBjet_mass;
  float gensubLeadingBjet_pt;
  float gensubLeadingBjet_eta;
  float gensubLeadingBjet_phi;
  float gensubLeadingBjet_mass;
  float genLeadingGenjet_pt;
  float genLeadingGenjet_eta;
  float genLeadingGenjet_phi;
  float genLeadingGenjet_mass;
  float gensubLeadingGenjet_pt;
  float gensubLeadingGenjet_eta;
  float gensubLeadingGenjet_phi;
  float gensubLeadingGenjet_mass;
    
  // added bjet matching
  float gen_matched_LeadingBjet_pt;
  float gen_matched_LeadingBjet_eta;
  float gen_matched_LeadingBjet_phi;
  float gen_matched_subLeadingBjet_pt;
  float gen_matched_subLeadingBjet_eta;
  float gen_matched_subLeadingBjet_phi;
  float gen_dibjet_mass;
  float gen_dibjet_pt;
  float gen_dibjet_eta;
    
  float genPho_deltaR;
  float genBjet_deltaR;
  float genBjet_H_pt;
  float genPho_H_pt;
  float genweight;   
    
  float recon;
  float bjet_recon;
  float photon_recon;
  float boostedCat;
  float VBFHH_recon;
    
  // added ML vars
  float leadingDeepBscore;
  float subleadingDeepBscore;
  float sumDeepBscore;
  float leading_pho_pt_over_dimass;
  float leading_bjet_pt_over_dimass;
  float leading_bjet_pt_over_dimass_corr;

  //add boosted object vars
  float fatJetPt;
  float fatJetEta;
  float fatJetPhi;
  float fatJetMassSD_UnCorrected;
  float fatJetbtagDDBvL;
  
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
 
  //luminosity for each year
  //https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
  luminosity["2016"] = 35.9;
  luminosity["2017"] = 41.5;
  luminosity["2018"] = 59.8;
  
  //xs in unit of fb
  //e.g. 
  //float xsHH = 31.05; //fb
  //float BRHbb = 5.824E-01;
  //float BRHgg = 2.270E-03;
  //ref https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MetaData/data/cross_sections.json   
  //ref https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGHH?redirectedfrom=LHCPhysics.LHCHXSWGHH#HHjj_VBF
  xs["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 31.05*5.824E-01*2.270E-03*2;
  xs["GluGluToHHTo2B2G_node_cHHH1_TuneCUETP8M1_PSWeights_13TeV-powheg-pythia8"] = 31.05*5.824E-01*2.270E-03*2;
  xs["VBFHHTo2B2G_CV_1_C2V_1_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8"] = 1.726*5.824E-01*2.270E-03*2;
  xs["VBFHHTo2B2G_CV_1_C2V_1_C3_1_13TeV-madgraph"] = 1.726*5.824E-01*2.270E-03*2;
  xs["VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"] = 2.257*0.00227*1.06*1000.;
  xs["ttHToGG_M125_TuneCP5_PSweights_13TeV-powheg-pythia8"] = 0.5071*0.00227*1.06*1000.;
  xs["ttHToGG_M125_13TeV_powheg_pythia8"] = 0.5071*0.00227*1.06*1000.;
  xs["ttHToGG_M125_13TeV_powheg_pythia8_v2"] = 0.5071*0.00227*1.06*1000.;
  xs["VBFHToGG_M125_13TeV_amcatnlo_pythia8"] = 3.7820*0.00227*1.06*1000.;
  xs["VBFHToGG_M-125_13TeV_powheg_pythia8"] = 3.7820*0.00227*1.06*1000.;
  xs["VBFHToGG_M125_13TeV_amcatnlo_pythia8_v2"] = 3.7820*0.00227*1.06*1000.;
  xs["GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 48.5800*0.00227*1.06*1000.;
  xs["GluGluHToGG_M-125_13TeV_powheg_pythia8"] = 48.5800*0.00227*1.06*1000.;
  xs["GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia"] = 48.5800*0.00227*1.06*1000.;
  xs["TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 831.76*1000.;
  xs["TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8"] = 831.76*1000.;
  xs["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 232.9*1000.;
  xs["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 220.0*1000.;
  xs["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 878.1*1000.;
  xs["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 850.8*1000.;
  xs["DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-Sherpa"] = 0.494*1000.;
  xs["DiPhotonJetsBox2BJets_MGG-80toInf_TuneSherpa_13TeV-Sherpa"] = 0.494*1000.;
  xs["DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-Sherpa"] = 0.8674276*1000.;
  xs["DiPhotonJetsBox1BJet_MGG-80toInf_TuneSherpa_13TeV-Sherpa"] = 0.8674276*1000.;
  xs["DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa"] = 84.4*1000.;
  xs["TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8"]=0.01687*1000.;
  xs["TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8"]=0.01731*1000.;
  xs["TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"]=4.078*1000.;
  xs["TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"]=3.819*1000.;
  xs["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 118100.0*1000.;
  xs["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 113400.0*1000.;
      
  std::map<std::string, float> sumOfgenw_2016, sumOfgenw_2017, sumOfgenw_2018;
    
  sumOfgenw_2018["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 5402.244961268351-122.629;
  sumOfgenw_2018["VBFHHTo2B2G_CV_1_C2V_1_C3_1_TuneCP5_PSWeights_13TeV-madgraph-pythia8"] = 392173.7132100609;
  sumOfgenw_2018["VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"] = 3800454.1885023806;
  sumOfgenw_2018["ttHToGG_M125_TuneCP5_PSweights_13TeV-powheg-pythia8"] = 526575.1838219386;
  sumOfgenw_2018["VBFHToGG_M125_13TeV_amcatnlo_pythia8"] = 7681928.604487639;
  sumOfgenw_2018["GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 213408144.61690226;
  sumOfgenw_2018["TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 638545951469.1536;
  sumOfgenw_2018["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 14366641.0;
  sumOfgenw_2018["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 10205533.0;
  sumOfgenw_2018["DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-Sherpa"] = 161251.90000000005;
  sumOfgenw_2018["DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-Sherpa"] = 169252.80000000002;
  sumOfgenw_2018["DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa"] = 6423331.299999999;
  sumOfgenw_2018["TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8"]= 50221.85598815847;
  sumOfgenw_2018["TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"]= 33438763.162006084;
  sumOfgenw_2018["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 10781997.056885984;    

  sumOfgenw_2017["DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-Sherpa"] = 179926.90000000002;
  sumOfgenw_2017["DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-Sherpa"] = 167403.10000000003;
  sumOfgenw_2017["DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa"] = 21638416.099999994;
  sumOfgenw_2017["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 17678702.0;
  sumOfgenw_2017["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 79243357.0;
  sumOfgenw_2017["GluGluHToGG_M-125_13TeV_powheg_pythia8"] = 20469841.991200007;
  sumOfgenw_2017["GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8"] = 5342.4187968381175 - 0.272939 + 8.68992;
  sumOfgenw_2017["VBFHHTo2B2G_CV_1_C2V_1_C3_1_13TeV-madgraph"] = 99078.000000;
  sumOfgenw_2017["TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8"] = 25145.499083880008;
  sumOfgenw_2017["TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"] = 51104339.85992301;
  sumOfgenw_2017["TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 318686315779.0374;
  sumOfgenw_2017["VBFHToGG_M-125_13TeV_powheg_pythia8"] = 3846447.1165199997;
  sumOfgenw_2017["VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"] = 4100171.411232001;
  sumOfgenw_2017["ttHToGG_M125_13TeV_powheg_pythia8"] = 504098.30106273957;
  sumOfgenw_2017["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8"] = 20622034.0;

  sumOfgenw_2016["DiPhotonJetsBox1BJet_MGG-80toInf_TuneSherpa_13TeV-Sherpa"] = 161957.1;
  sumOfgenw_2016["VBFHHTo2B2G_CV_1_C2V_1_C3_1_13TeV-madgraph"] = 3099141.4427206507;
  sumOfgenw_2016["DiPhotonJetsBox2BJets_MGG-80toInf_TuneSherpa_13TeV-Sherpa"] = 165812.60000000003;
  sumOfgenw_2016["DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa"] = 27856298.400000002;
  sumOfgenw_2016["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 23302380.0;
  sumOfgenw_2016["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 70829685.0;
  sumOfgenw_2016["GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia"] = 89212192.12336;
  sumOfgenw_2016["GluGluToHHTo2B2G_node_cHHH1_TuneCUETP8M1_PSWeights_13TeV-powheg-pythia8"] = 5372.566753895401 - 12.9627;
  sumOfgenw_2016["TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8"] = 25098.827877816002;
  sumOfgenw_2016["TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] = 92616271.19923201;
  sumOfgenw_2016["TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8"] = 96235105494.28423;
  sumOfgenw_2016["VBFHToGG_M125_13TeV_amcatnlo_pythia8_v2"] = 3908737.28136;
  sumOfgenw_2016["VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8"] = 1862202.435519;
  sumOfgenw_2016["ttHToGG_M125_13TeV_powheg_pythia8_v2"] = 431817.25939200015;
  sumOfgenw_2016["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 20816188.0;

  sumOfgenw["2016"] = sumOfgenw_2016;
  sumOfgenw["2017"] = sumOfgenw_2017;
  sumOfgenw["2018"] = sumOfgenw_2018;
  
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
  t_run = 0;
  t_luminosityBlock = 0;
  t_event = 0;
  trig_decision = 0;
  pv_pass = 0;
  leading_photon_pt = -999.;
  leading_photon_eta = -999.;
  leading_photon_phi = -999.;
  subleading_photon_pt = -999.;
  subleading_photon_eta = -999.;
  subleading_photon_phi = -999.;
  diphoton_pt = -999.;
  diphoton_eta = -999.;
  diphoton_mass = -999.;
  photon_delR = -999.;  
  // added for bjet reconstruction
  leading_bjet_pt = -999.;
  leading_bjet_eta = -999.;
  leading_bjet_phi = -999.;
  subleading_bjet_pt = -999.;
  subleading_bjet_eta = -999.;
  subleading_bjet_phi = -999.;
  dibjet_pt = -999.;
  dibjet_eta = -999.;
  dibjet_mass = -999.;
  bjet_delR = -999.;
  nbjet = -999.;
  // bjet corrections  
  leading_bjet_pt_corr = -999.;
  subleading_bjet_pt_corr = -999.;
  dibjet_pt_corr = -999.;
  dibjet_eta_corr = -999.;
  dibjet_mass_corr = -999.;

  // vbf jet
  leading_vbfjet_pt = -999.;
  leading_vbfjet_eta = -999.;
  leading_vbfjet_phi = -999.;
  subleading_vbfjet_pt = -999.;
  subleading_vbfjet_eta = -999.;
  subleading_vbfjet_phi = -999.;
  divbfjet_pt = -999.;
  divbfjet_eta = -999.;
  divbfjet_mass = -999.;
  vbfjet_delR = -999.;
  vbfjet_del_eta = -999.;
  nbjet = -999.;
    
  genHH_pt = -999.;
  genHH_eta = -999.;
  genHH_phi = -999.;
  genHH_mass = -999.;
  genLeadingH_pt = -999.;
  genLeadingH_eta = -999.;
  genLeadingH_phi = -999.;
  genLeadingH_mass = -999.;
  gensubLeadingH_pt = -999.;
  gensubLeadingH_eta = -999.;
  gensubLeadingH_phi = -999.;
  gensubLeadingH_mass = -999.;
  genLeadingPho_pt = -999.;
  genLeadingPho_eta = -999.;
  genLeadingPho_phi = -999.;
  genLeadingPho_mass = -999.;
  gensubLeadingPho_pt = -999.;
  gensubLeadingPho_eta = -999.;
  gensubLeadingPho_phi = -999.;
  gensubLeadingPho_mass = -999.;
  gen_matched_LeadingPho_pt = -999.;
  gen_matched_LeadingPho_eta = -999.;
  gen_matched_LeadingPho_phi = -999.;
  gen_matched_subLeadingPho_pt = -999.;
  gen_matched_subLeadingPho_eta = -999.;
  gen_matched_subLeadingPho_phi = -999.;
  genLeadingBjet_pt = -999.;
  genLeadingBjet_eta = -999.;
  genLeadingBjet_phi = -999.;
  genLeadingBjet_mass = -999.;
  gensubLeadingBjet_pt = -999.;
  gensubLeadingBjet_eta = -999.;
  gensubLeadingBjet_phi = -999.;
  gensubLeadingBjet_mass = -999.;
  genLeadingGenjet_pt = -999.;
  genLeadingGenjet_eta = -999.;
  genLeadingGenjet_phi = -999.;
  genLeadingGenjet_mass = -999.;
  gensubLeadingGenjet_pt = -999.;
  gensubLeadingGenjet_eta = -999.;
  gensubLeadingGenjet_phi = -999.;
  gensubLeadingGenjet_mass = -999.;
  // added bjet matching
  gen_matched_LeadingBjet_pt = -999.;
  gen_matched_LeadingBjet_eta = -999.;
  gen_matched_LeadingBjet_phi = -999.;
  gen_matched_subLeadingBjet_pt = -999.;
  gen_matched_subLeadingBjet_eta = -999.;
  gen_matched_subLeadingBjet_phi = -999.;
  gen_dibjet_mass = -999.;
  gen_dibjet_pt = -999.;
  gen_dibjet_eta = -999.;
  genPho_deltaR = -999.;
  genBjet_deltaR = -999.;
  genBjet_H_pt = -999.;
  genPho_H_pt = -999.;
  genweight = -999.;
    
  // recon?
  recon = 0.;
  bjet_recon = 0.;
  photon_recon = 0.;
  boostedCat = 0.;
  VBFHH_recon = 0.;
    
  // added ML vars
  leadingDeepBscore = -999.;
  subleadingDeepBscore = -999.;
  sumDeepBscore = -999.;
  leading_pho_pt_over_dimass = -999.;
  leading_bjet_pt_over_dimass = -999.;
  leading_bjet_pt_over_dimass_corr = -999.;

  //boosted category vars
  fatJetPt = -999.;
  fatJetEta = -999.;
  fatJetPhi = -999.;
  fatJetMassSD_UnCorrected = -999.;
  fatJetbtagDDBvL = -999.;

}

void HHbbggAnalyzer::BookTreeBranches(){
  tree = new TTree("tree","tree");
  //tree->SetAutoSave(10000);

  tree->Branch("run", &t_run,"run/i");
  tree->Branch("lumi", &t_luminosityBlock,"lumi/i");
  tree->Branch("event", &t_event,"event/l");
  tree->Branch("trig_decision", &trig_decision, "trig_decision/i");
  tree->Branch("pv_pass", &pv_pass, "pv_pass/i");    
      
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
  tree->Branch("photon_delR", &photon_delR,"photon_delR/f");
    
    
  //phone gen matched to reco information
  tree->Branch("gen_matched_LeadingPho_pt", &gen_matched_LeadingPho_pt, "gen_matched_LeadingPho_pt/f");
  tree->Branch("gen_matched_LeadingPho_eta", &gen_matched_LeadingPho_eta, "gen_matched_LeadingPho_eta/f");
  tree->Branch("gen_matched_LeadingPho_phi", &gen_matched_LeadingPho_phi, "gen_matched_LeadingPho_phi/f");
  tree->Branch("gen_matched_subLeadingPho_pt", &gen_matched_subLeadingPho_pt, "gen_matched_subLeadingPho_pt/f");
  tree->Branch("gen_matched_subLeadingPho_eta", &gen_matched_subLeadingPho_eta, "gen_matched_subLeadingPho_eta/f");
  tree->Branch("gen_matched_subLeadingPho_phi", &gen_matched_subLeadingPho_phi, "gen_matched_subLeadingPho_phi/f");
    
  //jets

  // added bjet
  tree->Branch("leading_bjet_pt", &leading_bjet_pt,"leading_bjet_pt/f");
  tree->Branch("leading_bjet_eta", &leading_bjet_eta,"leading_bjet_eta/f");
  tree->Branch("leading_bjet_phi", &leading_bjet_phi,"leading_bjet_phi/f");  
  tree->Branch("subleading_bjet_pt", &subleading_bjet_pt,"subleading_bjet_pt/f");
  tree->Branch("subleading_bjet_eta", &subleading_bjet_eta,"subleading_bjet_eta/f");
  tree->Branch("subleading_bjet_phi", &subleading_bjet_phi,"subleading_bjet_phi/f"); 
  tree->Branch("dibjet_pt", &dibjet_pt,"dibjet_pt/f"); 
  tree->Branch("dibjet_mass", &dibjet_mass,"dibjet_mass/f"); 
  tree->Branch("dibjet_eta", &dibjet_eta,"dibjet_eta/f");
  //bjet corrections
  tree->Branch("leading_bjet_pt_corr", &leading_bjet_pt_corr,"leading_bjet_pt_corr/f");
  tree->Branch("subleading_bjet_pt_corr", &subleading_bjet_pt_corr,"subleading_bjet_pt_corr/f");
  tree->Branch("dibjet_pt_corr", &dibjet_pt_corr,"dibjet_pt_corr/f"); 
  tree->Branch("dibjet_mass_corr", &dibjet_mass_corr,"dibjet_mass_corr/f"); 
  tree->Branch("dibjet_eta_corr", &dibjet_eta_corr,"dibjet_eta_corr/f");
  tree->Branch("bjet_delR", &bjet_delR,"bjet_delR/f");
  tree->Branch("nbjet", &nbjet,"nbjet/i");
    
  //bjet gen matched to reco information
  tree->Branch("gen_matched_LeadingBjet_pt", &gen_matched_LeadingBjet_pt, "gen_matched_LeadingBjet_pt/f");
  tree->Branch("gen_matched_LeadingBjet_eta", &gen_matched_LeadingBjet_eta, "gen_matched_LeadingBjet_eta/f");
  tree->Branch("gen_matched_LeadingBjet_phi", &gen_matched_LeadingBjet_phi, "gen_matched_LeadingBjet_phi/f");
  tree->Branch("gen_matched_subLeadingBjet_pt", &gen_matched_subLeadingBjet_pt, "gen_matched_subLeadingBjet_pt/f");
  tree->Branch("gen_matched_subLeadingBjet_eta", &gen_matched_subLeadingBjet_eta, "gen_matched_subLeadingBjet_eta/f");
  tree->Branch("gen_matched_subLeadingBjet_phi", &gen_matched_subLeadingBjet_phi, "gen_matched_subLeadingBjet_phi/f");
  tree->Branch("gen_dibjet_mass", &gen_dibjet_mass, "gen_dibjet_mass/f");
  tree->Branch("gen_dibjet_pt", &gen_dibjet_pt, "gen_dibjet_pt/f");
  tree->Branch("gen_dibjet_eta", &gen_dibjet_eta, "gen_dibjet_eta/f");

  // vbf jet
  tree->Branch("leading_vbfjet_pt", &leading_vbfjet_pt,"leading_vbfjet_pt/f");
  tree->Branch("leading_vbfjet_eta", &leading_vbfjet_eta,"leading_vbfjet_eta/f");
  tree->Branch("leading_vbfjet_phi", &leading_vbfjet_phi,"leading_vbfjet_phi/f");  
  tree->Branch("subleading_vbfjet_pt", &subleading_vbfjet_pt,"subleading_vbfjet_pt/f");
  tree->Branch("subleading_vbfjet_eta", &subleading_vbfjet_eta,"subleading_vbfjet_eta/f");
  tree->Branch("subleading_vbfjet_phi", &subleading_vbfjet_phi,"subleading_vbfjet_phi/f"); 
  tree->Branch("divbfjet_pt", &divbfjet_pt,"divbfjet_pt/f"); 
  tree->Branch("divbfjet_mass", &dibjet_mass,"divbfjet_mass/f"); 
  tree->Branch("divbfjet_eta", &dibjet_eta,"divbfjet_eta/f");
  tree->Branch("vbfjet_delR", &vbfjet_delR,"vbfjet_delR/f");
  tree->Branch("vbfjet_del_eta", &vbfjet_del_eta,"vbfjet_del_eta/f");
  tree->Branch("nvbfjet", &nvbfjet,"nvbfjet/i");
    
  //Gen information
  tree->Branch("genHH_pt", &genHH_pt,"genHH_pt/f"); 
  tree->Branch("genHH_mass", &genHH_mass,"genHH_mass/f"); 
  tree->Branch("genHH_eta", &genHH_eta,"genHH_eta/f");
  tree->Branch("genHH_phi", &genHH_phi,"genHH_phi/f");
  tree->Branch("genLeadingH_pt", &genLeadingH_pt,"genLeadingH_pt/f");
  tree->Branch("genLeadingH_mass", &genLeadingH_mass,"genLeadingH_mass/f");
  tree->Branch("genLeadingH_eta", &genLeadingH_eta,"genLeadingH_eta/f");
  tree->Branch("genLeadingH_phi", &genLeadingH_phi,"genLeadingH_phi/f");
  tree->Branch("gensubLeadingH_pt", &gensubLeadingH_pt,"gensubLeadingH_pt/f");
  tree->Branch("gensubLeadingH_mass", &gensubLeadingH_mass,"gensubLeadingH_mass/f");
  tree->Branch("gensubLeadingH_eta", &gensubLeadingH_eta,"gensubLeadingH_eta/f");
  tree->Branch("gensubLeadingH_phi", &gensubLeadingH_phi,"gensubLeadingH_phi/f");
  tree->Branch("genLeadingPho_pt", &genLeadingPho_pt,"genLeadingPho_pt/f");
  tree->Branch("genLeadingPho_mass", &genLeadingPho_mass,"genLeadingPho_mass/f");
  tree->Branch("genLeadingPho_eta", &genLeadingPho_eta,"genLeadingPho_eta/f");
  tree->Branch("genLeadingPho_phi", &genLeadingPho_phi,"genLeadingPho_phi/f");
  tree->Branch("gensubLeadingPho_pt", &gensubLeadingPho_pt,"gensubLeadingPho_pt/f");
  tree->Branch("gensubLeadingPho_mass", &gensubLeadingPho_mass,"gensubLeadingPho_mass/f");
  tree->Branch("gensubLeadingPho_eta", &gensubLeadingPho_eta,"gensubLeadingPho_eta/f");
  tree->Branch("gensubLeadingPho_phi", &gensubLeadingPho_phi,"gensubLeadingPho_phi/f");
  tree->Branch("genLeadingBjet_pt", &genLeadingBjet_pt,"genLeadingBjet_pt/f");
  tree->Branch("genLeadingBjet_mass", &genLeadingBjet_mass,"genLeadingBjet_mass/f");
  tree->Branch("genLeadingBjet_eta", &genLeadingBjet_eta,"genLeadingBjet_eta/f");
  tree->Branch("genLeadingBjet_phi", &genLeadingBjet_phi,"genLeadingBjet_phi/f");
  tree->Branch("gensubLeadingBjet_pt", &gensubLeadingBjet_pt,"gensubLeadingBjet_pt/f");
  tree->Branch("gensubLeadingBjet_mass", &gensubLeadingBjet_mass,"gensubLeadingBjet_mass/f");
  tree->Branch("gensubLeadingBjet_eta", &gensubLeadingBjet_eta,"gensubLeadingBjet_eta/f");
  tree->Branch("gensubLeadingBjet_phi", &gensubLeadingBjet_phi,"gensubLeadingBjet_phi/f");
  tree->Branch("genLeadingGenjet_pt", &genLeadingGenjet_pt,"genLeadingGenjet_pt/f");
  tree->Branch("genLeadingGenjet_mass", &genLeadingGenjet_mass,"genLeadingGenjet_mass/f");
  tree->Branch("genLeadingGenjet_eta", &genLeadingGenjet_eta,"genLeadingGenjet_eta/f");
  tree->Branch("genLeadingGenjet_phi", &genLeadingGenjet_phi,"genLeadingGenjet_phi/f");
  tree->Branch("gensubLeadingGenjet_pt", &gensubLeadingGenjet_pt,"gensubLeadingGenjet_pt/f");
  tree->Branch("gensubLeadingGenjet_mass", &gensubLeadingGenjet_mass,"gensubLeadingGenjet_mass/f");
  tree->Branch("gensubLeadingGenjet_eta", &gensubLeadingGenjet_eta,"gensubLeadingGenjet_eta/f");
  tree->Branch("gensubLeadingGenjet_phi", &gensubLeadingGenjet_phi,"gensubLeadingGenjet_phi/f");
  tree->Branch("genPho_deltaR", &genPho_deltaR, "genPho_deltaR/f");
  tree->Branch("genBjet_deltaR", &genBjet_deltaR, "genBjet_deltaR/f");
  tree->Branch("genBjet_H_pt", &genBjet_H_pt, "genbBjet_H_pt/f");
  tree->Branch("genPho_H_pt", &genPho_H_pt, "genPho_H_pt/f");
  
  //event weight (genweight) and other reco level information
  tree->Branch("genweight", &genweight, "genweight/f");
  tree->Branch("recon", &recon, "recon/f");
  tree->Branch("bjet_recon", &bjet_recon, "bjet_recon/f");
  tree->Branch("photon_recon", &photon_recon, "photon_recon/f");
  tree->Branch("boostedCat", &boostedCat, "boostedCat/f");
  tree->Branch("VBFHH_recon", &VBFHH_recon, "VBFHH_recon/f");
  tree->Branch("leadingDeepBscore", &leadingDeepBscore, "leadingDeepBscore/f");
  tree->Branch("subleadingDeepBscore", &subleadingDeepBscore, "subleadingDeepBscore/f");
  tree->Branch("sumDeepBscore", &sumDeepBscore, "sumDeepBscore/f");
  tree->Branch("leading_pho_pt_over_dimass", &leading_pho_pt_over_dimass, "leading_pho_pt_over_dimass/f");
  tree->Branch("leading_bjet_pt_over_dimass", &leading_bjet_pt_over_dimass, "leading_bjet_pt_over_dimass/f");
  tree->Branch("leading_bjet_pt_over_dimass_corr", &leading_bjet_pt_over_dimass_corr, "leading_bjet_pt_over_dimass_corr/f");

  //boosted category  
  tree->Branch("fatJetPt", &fatJetPt, "fatJetPt/f");
  tree->Branch("fatJetEta", &fatJetEta, "fatJetEta/f");
  tree->Branch("fatJetPhi", &fatJetPhi, "fatJetPhi/f");
  tree->Branch("fatJetMassSD_UnCorrected", &fatJetMassSD_UnCorrected, "fatJetMassSD_UnCorrected/f");
  tree->Branch("fatJetbtagDDBvL", &fatJetbtagDDBvL, "fatJetbtagDDBvL"); 
}
#endif // #ifdef HHbbggAnalyzer_cxx
