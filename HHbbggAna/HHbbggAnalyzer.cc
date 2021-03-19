#define HHbbggAnalyzer_cxx
#include "HHbbggAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

#ifdef MAKECINT
#pragma link C++ class vector<float>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<int>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<bool>+;
#endif

int main(int argc, char* argv[])
{

  if(argc < 4) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type and year"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  TString year_num          = argv[5];
  HHbbggAnalyzer Hmm(inputFileList, outFileName, data, isData, year_num);
  cout << "dataset " << data << " year " <<year_num<< endl;
  Hmm.EventLoop(data, isData);

  return 0;
}

void HHbbggAnalyzer::EventLoop(const char *data,const char *isData)
{ 
  if (fChain == 0) return;
  //clearTreeVectors();
  //cout<<"cleared tree vectors\n";
  //BookTreeBranches();
  //cout<<"booked tree branches\n";
  float muon_mass = 0.1056583745;
  bool datafile = true;
  if(*isData=='F') datafile = false;
 
  //btag SF
  BTagCalibration calib("deepcsv","data/btagSF/DeepCSV_94XSF_V3_B_F.csv");
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,  // operating point
			       "central",             // central sys type
			       {"up", "down"});      // other sys types

  reader.load(calib,                // calibration instance
	      BTagEntry::FLAV_B,    // btag flavour
	      "comb")       ;        // measurement type
  
   Long64_t nentries = fChain->GetEntriesFast();
   cout <<"total entries: "<<nentries<<endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%5000==0) cout <<"entry: "<<jentry<<endl;
      clearTreeVectors();

      //sum of genWeight
      float value_h_sumOfgw = h_sumOfgw->GetBinContent(1);
      if(*isData=='F')   value_h_sumOfgw = value_h_sumOfgw + genWeight;
      else value_h_sumOfgw = value_h_sumOfgw + 1.0;
      h_sumOfgw->SetBinContent(1,value_h_sumOfgw);

      //sum of genWeight and pileupweight
      float value_h_sumOfgpw = h_sumOfgpw->GetBinContent(1);
      if(*isData=='F')   value_h_sumOfgpw = value_h_sumOfgpw + genWeight;
      else value_h_sumOfgpw = value_h_sumOfgpw + 1.0;
      h_sumOfgpw->SetBinContent(1,value_h_sumOfgpw);

      bool trig_decision = false;
      //• 2016 : HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v* 2017 : HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*
      //trigger for 2016 and 2017 to be changed in the code
      if( year=="2016" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
      if( year=="2017" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
      if( year=="2018" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
       
      //photon selection
      int index_ph1(-999), index_ph2(-999);
      vector<int> index_photon;
      index_photon.clear();
      bool Event_sel = false;
      for(int i=0;i<nPhoton;i++){
          if(Photon_mvaID_WP90[i]==1 && (Photon_pt[i] > 25) );
          bool eta_cut = (year=="2016" ? (fabs(Photon_eta[i]) < 2.4) : (fabs(Photon_eta[i]) < 2.5)) && ((fabs(Photon_eta[i])<1.44 || fabs(Photon_eta[i])>1.57)); 
          if(eta_cut) index_photon.push_back(i);
      }
      if(index_photon.size()<2) continue;
      else if(index_photon.size()==2){
          TLorentzVector photon_1, photon_2, diphoton;
          photon_1.SetPtEtaPhiM(Photon_pt[0],Photon_eta[0],Photon_phi[0],0);
          photon_2.SetPtEtaPhiM(Photon_pt[1],Photon_eta[1],Photon_phi[1],0);
          diphoton = photon_1 + photon_2;
          if(Photon_pt[0]/diphoton.M()>1/3 && Photon_pt[1]/diphoton.M()>1/4){
              index_ph1 = index_photon.at(0);
              index_ph2 = index_photon.at(1);
              diphoton_pt = diphoton.Pt();
              diphoton_mass = diphoton.M();
              diphoton_eta = diphoton.Eta();
          }
          else continue;
      }
      else{
          float tmp_diphoton_pt = 0.;
          for(int j=0; j<index_photon.size(); j++){
              for(int k=j+1; k<index_photon.size(); k++){
                  TLorentzVector photon_1, photon_2, diphoton;
                  photon_1.SetPtEtaPhiM(Photon_pt[j],Photon_eta[j],Photon_phi[j],0);
                  photon_2.SetPtEtaPhiM(Photon_pt[k],Photon_eta[k],Photon_phi[k],0);
                  diphoton = photon_1 + photon_2;
                  if(Photon_pt[j]/diphoton.M()>1/3 && Photon_pt[k]/diphoton.M()>1/4 && tmp_diphoton_pt < diphoton.Pt()){
                      tmp_diphoton_pt = diphoton.Pt();
                      index_ph1 = j;
                      index_ph2 = k;
                      diphoton_pt = tmp_diphoton_pt;
                      diphoton_mass = diphoton.M();
                      diphoton_eta = diphoton.Eta();
                  }
              }
          }
      }
      //jet selection
      
      t_run =run;
      t_luminosityBlock=luminosityBlock;
      t_event=event;
      leading_photon_pt=Photon_pt[index_ph1];
      leading_photon_eta=Photon_eta[index_ph1];
      leading_photon_phi=Photon_phi[index_ph1];
      subleading_photon_pt=Photon_pt[index_ph2];
      subleading_photon_eta=Photon_eta[index_ph2];
      subleading_photon_phi=Photon_phi[index_ph2];
    
      tree->Fill();      
   }
}
