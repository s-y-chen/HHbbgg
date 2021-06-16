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
  const char *isRunGen      = argv[5];
  TString year_num          = argv[6];
  HHbbggAnalyzer Hmm(inputFileList, outFileName, data, isData, year_num);
  cout << "dataset " << data << " year " <<year_num<< endl;
  Hmm.EventLoop(data, isData, isRunGen);

  return 0;
}

void HHbbggAnalyzer::EventLoop(const char *data, const char *isData, const char *isRunGen)
{ 
  if (fChain == 0) return;
  //clearTreeVectors();
  //cout<<"cleared tree vectors\n";
  //BookTreeBranches();
  //cout<<"booked tree branches\n";
  float muon_mass = 0.1056583745;

  float xsHH = 31.05; //fb
  float BRHbb = 5.824E-01;
  float BRHgg = 2.270E-03;
  float lumi = 300; //fb-1
      
  bool datafile = true;
  if(*isData=='F') datafile = false;

  bool debug = false;
    
  const int nHpTbin = 10;
  float HpT_bounds[nHpTbin+1] = {0, 100, 200, 300, 400, 500, 600, 700, 1000, 2000, 3000};
  float all_events[nHpTbin] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  float pass_events[nHpTbin] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 
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
      if(jentry%1000==0) cout <<"entry: "<<jentry<<endl;
      clearTreeVectors();

      bool skip = false;
      genweight = genWeight*xsHH*BRHbb*BRHgg*2.*lumi;
      if(fabs(genweight)>5){
          cout <<"event with abnormal large weight: event run "<<event<<" "<<run<<endl;
          continue;
      } 
      
      if(!(event==20522 && run==1) && debug) continue;
     
      //gen information: 
      //pdg ID scheme: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
      //0<= Jet_genJetIdx < nGenJet.
      //if(*isRunGen=='T'){
          int h1_index_tmp = -999, h2_index_tmp = -999;
          int photon1_index_tmp = -999, photon2_index_tmp = -999;
          int bjet1_index_tmp = -999, bjet2_index_tmp = -999, genjet1_index_tmp = -999, genjet2_index_tmp = -999;
       
          vector<int> index_higgs, index_diphoton, index_dibjet, index_digenjet;
          index_higgs.clear();
          index_diphoton.clear();
          index_dibjet.clear();
          index_digenjet.clear();
       
          for(int igenpart=0; igenpart<nGenPart; igenpart++){    
              if(GenPart_pdgId[igenpart]==22 && GenPart_pdgId[GenPart_genPartIdxMother[igenpart]]==25 && GenPart_status[GenPart_genPartIdxMother[igenpart]]==62){ 
                  genPho_H_pt = GenPart_pt[GenPart_genPartIdxMother[igenpart]];
                  for(int i=0; i<nHpTbin; i++){
                      if(HpT_bounds[i]<genPho_H_pt &&  HpT_bounds[i+1]>genPho_H_pt) all_events[i] = all_events[i]  + genweight;
                  }
                  
                  index_diphoton.push_back(igenpart);
                  if(debug){
                      cout<<"photon status "<<GenPart_status[igenpart]
                      <<" IdxMother "<<GenPart_genPartIdxMother[igenpart]
                      <<" pt "<<GenPart_pt[igenpart]
                      <<" mass "<<GenPart_mass[igenpart]<<endl;
                      cout <<" mother pdgid"<< GenPart_pdgId[GenPart_genPartIdxMother[igenpart]]<< " "<<GenPart_pt[GenPart_genPartIdxMother[igenpart]] <<" "<<GenPart_eta[GenPart_genPartIdxMother[igenpart]]<<endl;
                  }
              }
              if(fabs(GenPart_pdgId[igenpart])==5 && GenPart_pdgId[GenPart_genPartIdxMother[igenpart]]==25 && GenPart_status[GenPart_genPartIdxMother[igenpart]]==62){ 
                  genBjet_H_pt = GenPart_pt[GenPart_genPartIdxMother[igenpart]];
                  index_dibjet.push_back(igenpart);
                  float dr_min = 999.;
                  int index_getjet = -1;
                  for(int igenjet=0; igenjet<nGenJet; igenjet++){
                      float dr_tmp = DeltaR(GenJet_eta[igenjet], GenJet_phi[igenjet], GenPart_eta[igenpart], GenPart_phi[igenpart]);
                      if(dr_min>dr_tmp){
                          dr_min = dr_tmp;
                          index_getjet = igenjet;
                      }
                  }
                  index_digenjet.push_back(index_getjet);
                      
                  if(debug){
                      cout<<"b-quark status "<<GenPart_status[igenpart]
                      <<" IdxMother "<<GenPart_genPartIdxMother[igenpart]
                      <<" pt "<<GenPart_pt[igenpart]
                      <<" eta "<<GenPart_eta[igenpart]
                      <<" phi "<<GenPart_phi[igenpart]
                      <<" mass "<<GenPart_mass[igenpart]<<endl;
                      cout <<" mother pdgid"<< GenPart_pdgId[GenPart_genPartIdxMother[igenpart]]<< " "<<GenPart_pt[GenPart_genPartIdxMother[igenpart]] <<" "<<GenPart_eta[GenPart_genPartIdxMother[igenpart]]<<endl;
                  }
               }
              if(GenPart_pdgId[igenpart]==25 && GenPart_status[igenpart]==62) index_higgs.push_back(igenpart);
          }
          //fill in the gen information for the two Higgs bosons
          if(index_higgs.size()==2){
              if(GenPart_pt[index_higgs.at(0)] > GenPart_pt[index_higgs.at(1)]){
                  h1_index_tmp = index_higgs.at(0);
                  h2_index_tmp = index_higgs.at(1);
              }
              else{
                  h1_index_tmp = index_higgs.at(1);
                  h2_index_tmp = index_higgs.at(0);
              }
          }
          else{
              cout <<"N(GEN-Higgs)!=2: "<<index_higgs.size()<<endl;
              skip = true;
          }
          if(h1_index_tmp > -999 && h2_index_tmp > -999){
              genLeadingH_pt = GenPart_pt[h1_index_tmp];
              genLeadingH_eta = GenPart_eta[h1_index_tmp];
              genLeadingH_phi = GenPart_phi[h1_index_tmp];
              genLeadingH_mass = GenPart_mass[h1_index_tmp];
              
              gensubLeadingH_pt = GenPart_pt[h2_index_tmp];
              gensubLeadingH_eta = GenPart_eta[h2_index_tmp];
              gensubLeadingH_phi = GenPart_phi[h2_index_tmp];
              gensubLeadingH_mass = GenPart_mass[h2_index_tmp];
              
              TLorentzVector H1, H2, HH;
              H1.SetPtEtaPhiM(GenPart_pt[h1_index_tmp], GenPart_eta[h1_index_tmp], GenPart_phi[h1_index_tmp], GenPart_mass[h1_index_tmp]);
              H2.SetPtEtaPhiM(GenPart_pt[h2_index_tmp], GenPart_eta[h2_index_tmp], GenPart_phi[h2_index_tmp], GenPart_mass[h2_index_tmp]);
              HH = H1 + H2;
              genHH_pt = HH.Pt();
              genHH_phi = HH.Phi();
              genHH_eta = HH.Eta();
              genHH_mass = HH.M();                        
          }
          //fill in the gen information for the photon pairs
          if(index_diphoton.size()==2){
              if(GenPart_pt[index_diphoton.at(0)] > GenPart_pt[index_diphoton.at(1)]){
                  photon1_index_tmp = index_diphoton.at(0);
                  photon2_index_tmp = index_diphoton.at(1);
              }
              else{
                  photon1_index_tmp = index_diphoton.at(1);
                  photon2_index_tmp = index_diphoton.at(0);
              }
          }
          else{
              cout <<"N(GEN-Photon)!=2: "<<index_diphoton.size()<<endl;
              skip = true;
          }
          if(photon1_index_tmp > -999 && photon2_index_tmp > -999){
              genLeadingPho_pt = GenPart_pt[photon1_index_tmp];
              genLeadingPho_eta = GenPart_eta[photon1_index_tmp];
              genLeadingPho_phi = GenPart_phi[photon1_index_tmp];
              genLeadingPho_mass = GenPart_mass[photon1_index_tmp];        
              gensubLeadingPho_pt = GenPart_pt[photon2_index_tmp];
              gensubLeadingPho_eta = GenPart_eta[photon2_index_tmp];
              gensubLeadingPho_phi = GenPart_phi[photon2_index_tmp];
              gensubLeadingPho_mass = GenPart_mass[photon2_index_tmp]; 
              genPho_deltaR = DeltaR(GenPart_eta[photon1_index_tmp],GenPart_phi[photon1_index_tmp], GenPart_eta[photon2_index_tmp], GenPart_phi[photon2_index_tmp]);
          }
          //fill in the gen information for the bjet pairs
          if(index_dibjet.size()==2){
              if(GenPart_pt[index_dibjet.at(0)] > GenPart_pt[index_dibjet.at(1)]){
                  bjet1_index_tmp = index_dibjet.at(0);
                  bjet2_index_tmp = index_dibjet.at(1);
                  genjet1_index_tmp = index_digenjet.at(0);
                  genjet2_index_tmp = index_digenjet.at(1);
              }
              else{
                  bjet1_index_tmp = index_dibjet.at(1);
                  bjet2_index_tmp = index_dibjet.at(0);
                  genjet1_index_tmp = index_digenjet.at(1);
                  genjet2_index_tmp = index_digenjet.at(0);
              }
          }
          else{
              cout <<"N(GEN-BJet)!=2: "<<index_dibjet.size()<<endl;
              skip = true;
          }
          if(bjet1_index_tmp > -999 && bjet2_index_tmp > -999){
              genLeadingBjet_pt = GenPart_pt[bjet1_index_tmp];
              genLeadingBjet_eta = GenPart_eta[bjet1_index_tmp];
              genLeadingBjet_phi = GenPart_phi[bjet1_index_tmp];
              genLeadingBjet_mass = GenPart_mass[bjet1_index_tmp];       
              gensubLeadingBjet_pt = GenPart_pt[bjet2_index_tmp];
              gensubLeadingBjet_eta = GenPart_eta[bjet2_index_tmp];
              gensubLeadingBjet_phi = GenPart_phi[bjet2_index_tmp];
              gensubLeadingBjet_mass = GenPart_mass[bjet2_index_tmp];
              genBjet_deltaR = DeltaR(genLeadingBjet_eta, genLeadingBjet_phi, gensubLeadingBjet_eta, gensubLeadingBjet_phi);
              if(genjet1_index_tmp>=0){
                  genLeadingGenjet_pt = GenJet_pt[genjet1_index_tmp];
                  genLeadingGenjet_eta = GenJet_eta[genjet1_index_tmp];
                  genLeadingGenjet_phi = GenJet_phi[genjet1_index_tmp];
                  genLeadingGenjet_mass = GenJet_mass[genjet1_index_tmp];
              }
              if(genjet2_index_tmp>=0){ 
                  gensubLeadingGenjet_pt = GenJet_pt[genjet2_index_tmp];
                  gensubLeadingGenjet_eta = GenJet_eta[genjet2_index_tmp];
                  gensubLeadingGenjet_phi = GenJet_phi[genjet2_index_tmp];
                  gensubLeadingGenjet_mass = GenJet_mass[genjet2_index_tmp];
              }
          }  
       
      if(skip) continue;
       
      //sum of genWeight
      float value_h_sumOfgw = h_sumOfgw->GetBinContent(1);
      if(*isData=='F' && fabs(genWeight)<1.)   value_h_sumOfgw = value_h_sumOfgw + genWeight;
      else value_h_sumOfgw = value_h_sumOfgw + 1.0;
      h_sumOfgw->SetBinContent(1,value_h_sumOfgw);

      //sum of genWeight and pileupweight
      float value_h_sumOfgpw = h_sumOfgpw->GetBinContent(1);
      if(*isData=='F' && fabs(genWeight)<10)   value_h_sumOfgpw = value_h_sumOfgpw + genWeight;
      else value_h_sumOfgpw = value_h_sumOfgpw + 1.0;
      h_sumOfgpw->SetBinContent(1,value_h_sumOfgpw);
       
      //}
      //else{ 
          
      bool trig_decision = false;
      //• 2016 : HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v* 2017 : HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*
      //trigger for 2016 and 2017 to be changed in the code
      if( year=="2016" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
      if( year=="2017" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
      if( year=="2018" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision =true;
       
       
      //photon selection
      int index_ph1(-999), index_ph2(-999);
      //match reco photon to genPhoton
      int gen_index1_matched_reco_photon = -999, gen_index2_matched_reco_photon = -999;
       
      vector<int> index_photon;
      index_photon.clear();
      bool Event_sel = false;
      for(int i=0;i<nPhoton;i++){
          if(Photon_mvaID_WP90[i]==1 && Photon_pt[i] > 25){ 
              bool eta_cut = (year=="2016" ? (fabs(Photon_eta[i]) < 2.4) : (fabs(Photon_eta[i]) < 2.5)) && ((fabs(Photon_eta[i])<1.44 || fabs(Photon_eta[i])>1.57)); 
              if(eta_cut) index_photon.push_back(i);
          } 
      }
      //if(index_photon.size()<2) continue;
      if(index_photon.size()>1 && trig_decision){
          if(index_photon.size()==2){
              TLorentzVector photon_1, photon_2, diphoton;
              photon_1.SetPtEtaPhiM(Photon_pt[0],Photon_eta[0],Photon_phi[0],0);
              photon_2.SetPtEtaPhiM(Photon_pt[1],Photon_eta[1],Photon_phi[1],0);
              diphoton = photon_1 + photon_2;
              if(Photon_pt[0]/diphoton.M()>1/3 && Photon_pt[1]/diphoton.M()>1/4 && 100 < diphoton.M() && diphoton.M() < 180){
                  index_ph1 = index_photon.at(0);
                  index_ph2 = index_photon.at(1);
                  diphoton_pt = diphoton.Pt();
                  diphoton_mass = diphoton.M();
                  diphoton_eta = diphoton.Eta();
              }
          }
          //else continue;
      
          else if(index_photon.size()>2){
              float tmp_diphoton_pt = 0.;
              for(int j=0; j<index_photon.size(); j++){
                  for(int k=j+1; k<index_photon.size(); k++){
                      TLorentzVector photon_1, photon_2, diphoton;
                      photon_1.SetPtEtaPhiM(Photon_pt[j],Photon_eta[j],Photon_phi[j],0);
                      photon_2.SetPtEtaPhiM(Photon_pt[k],Photon_eta[k],Photon_phi[k],0);
                      diphoton = photon_1 + photon_2;
                      if(Photon_pt[j]/diphoton.M()>1/3 && Photon_pt[k]/diphoton.M()>1/4 && tmp_diphoton_pt < diphoton.Pt() && 100 < diphoton.M() && diphoton.M() < 180){
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
      }
       
      //the first index is gen, second is reco
      float dR11 = DeltaR(GenPart_eta[photon1_index_tmp], GenPart_phi[photon1_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
      float dR12 = DeltaR(GenPart_eta[photon1_index_tmp], GenPart_phi[photon1_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
      float dR21 = DeltaR(GenPart_eta[photon2_index_tmp], GenPart_phi[photon2_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
      float dR22 = DeltaR(GenPart_eta[photon2_index_tmp], GenPart_phi[photon2_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
       
      if(dR11<0.4 && dR11<dR21){
          gen_index1_matched_reco_photon = photon1_index_tmp;
          if(dR22<0.4) gen_index2_matched_reco_photon = photon2_index_tmp;
          
      } 
      else if(dR21<0.4 && dR11>dR21){
          gen_index1_matched_reco_photon = photon2_index_tmp;
          if(dR12<0.4) gen_index2_matched_reco_photon = photon1_index_tmp;
      }    
      //efficiency vs Higgs pT
      for(int i=0; i<nHpTbin; i++){
          if(HpT_bounds[i]<genPho_H_pt &&  HpT_bounds[i+1]>genPho_H_pt) pass_events[i] = pass_events[i]  + genweight;
      }
      
      //jet selection
      int index_bj1(-999), index_bj2(-999);
      //match reco bjet to genbjet
      int gen_index1_matched_reco_bjet = -999, gen_index2_matched_reco_bjet = -999;
       
      vector<int> index_bjet;
      index_bjet.clear();
      Event_sel = false;
      for(int i=0;i<nJet;i++){
          if(Jet_btagDeepB[i]>= 0.5 && (Jet_pt[i] > 25)){
              bool eta_cut = (year=="2016" ? (fabs(Jet_eta[i]) < 2.4) : (fabs(Jet_eta[i]) < 2.5)); 
              if(eta_cut) index_bjet.push_back(i);
          }
      }
      //if(index_bjet.size()<2) continue;
      if(index_bjet.size()>1 && trig_decision){
          if(index_bjet.size()==2){
              TLorentzVector bjet_1, bjet_2, dibjet;
              bjet_1.SetPtEtaPhiM(Jet_pt[0],Jet_eta[0],Jet_phi[0],0);
              bjet_2.SetPtEtaPhiM(Jet_pt[1],Jet_eta[1],Jet_phi[1],0);
              dibjet = bjet_1 + bjet_2;
              if(70 < dibjet.M() && dibjet.M() < 190){
                  index_bj1 = index_bjet.at(0);
                  index_bj2 = index_bjet.at(1);
                  dibjet_pt = dibjet.Pt();
                  dibjet_mass = dibjet.M();
                  dibjet_eta = dibjet.Eta();
              }
          }
          else if(index_bjet.size()>2){
              // want highest b tagging scores
              float b_score_sum = 0.;
              for(int j=0; j<index_bjet.size(); j++){
                  for(int k=j+1; k<index_bjet.size(); k++){
                      float b_score_sum_jk = 0.;
                      b_score_sum_jk = Jet_btagDeepB[j] + Jet_btagDeepB[k];
                      TLorentzVector bjet_1, bjet_2, dibjet;
                      bjet_1.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],0);
                      bjet_2.SetPtEtaPhiM(Jet_pt[k],Jet_eta[k],Jet_phi[k],0);
                      dibjet = bjet_1 + bjet_2;
                      if(b_score_sum_jk > b_score_sum){
                          b_score_sum = b_score_sum_jk;
                          index_bj1 = j;
                          index_bj2 = k;
                          dibjet_pt = dibjet.Pt();
                          dibjet_mass = dibjet.M();
                          dibjet_eta = dibjet.Eta();
                      }
                  }
              }
          }
      }   
          
      //the first index is gen bjet, second is photon
      float dR11bg = DeltaR(GenPart_eta[bjet1_index_tmp], GenPart_phi[bjet1_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
      float dR12bg = DeltaR(GenPart_eta[bjet1_index_tmp], GenPart_phi[bjet1_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
      float dR21bg = DeltaR(GenPart_eta[bjet2_index_tmp], GenPart_phi[bjet2_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
      float dR22bg = DeltaR(GenPart_eta[bjet2_index_tmp], GenPart_phi[bjet2_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
          
      // first index if gen bjet, second is recon bjet
      float dR11bb = DeltaR(GenPart_eta[bjet1_index_tmp], GenPart_phi[bjet1_index_tmp], Jet_eta[index_bj1], Jet_phi[index_bj1]);
      float dR12bb = DeltaR(GenPart_eta[bjet1_index_tmp], GenPart_phi[bjet1_index_tmp], Jet_eta[index_bj2], Jet_phi[index_bj2]);
      float dR21bb = DeltaR(GenPart_eta[bjet2_index_tmp], GenPart_phi[bjet2_index_tmp], Jet_eta[index_bj1], Jet_phi[index_bj1]);
      float dR22bb = DeltaR(GenPart_eta[bjet2_index_tmp], GenPart_phi[bjet2_index_tmp], Jet_eta[index_bj2], Jet_phi[index_bj2]);
      if(dR11bg > 0.4 && dR12bg > 0.4 && dR21bg > 0.4 && dR22bg > 0.4){
          if(dR11bb<dR21bb){
              gen_index1_matched_reco_bjet = bjet1_index_tmp;
              gen_index2_matched_reco_bjet = bjet2_index_tmp;
          
          } 
          else {
              gen_index1_matched_reco_bjet = bjet2_index_tmp;
              gen_index2_matched_reco_bjet = bjet1_index_tmp;
          }  
      }
      
      //trigger if statement
       
      t_run =run;
      t_luminosityBlock=luminosityBlock;
      t_event=event;
      if(index_ph1>=0){
          leading_photon_pt=Photon_pt[index_ph1];
          leading_photon_eta=Photon_eta[index_ph1];
          leading_photon_phi=Photon_phi[index_ph1];
      }
      if(index_ph2>=0){
          subleading_photon_pt=Photon_pt[index_ph2];
          subleading_photon_eta=Photon_eta[index_ph2];
          subleading_photon_phi=Photon_phi[index_ph2]; 
      }
      if(gen_index1_matched_reco_photon>=0){
            gen_matched_LeadingPho_pt = GenPart_pt[gen_index1_matched_reco_photon];
            gen_matched_LeadingPho_eta = GenPart_eta[gen_index1_matched_reco_photon];
            gen_matched_LeadingPho_phi = GenPart_phi[gen_index1_matched_reco_photon];
      }
      if(gen_index2_matched_reco_photon>=0){
            gen_matched_subLeadingPho_pt = GenPart_pt[gen_index2_matched_reco_photon];
            gen_matched_subLeadingPho_eta = GenPart_eta[gen_index2_matched_reco_photon];
            gen_matched_subLeadingPho_phi = GenPart_phi[gen_index2_matched_reco_photon];           
      }
      // bjet
      if(index_bj1>=0){
          leading_bjet_pt=Jet_pt[index_bj1];
          leading_bjet_eta=Jet_eta[index_bj1];
          leading_bjet_phi=Jet_phi[index_bj1];
      }
      if(index_bj2>=0){
          subleading_bjet_pt=Jet_pt[index_bj2];
          subleading_bjet_eta=Jet_eta[index_bj2];
          subleading_bjet_phi=Jet_phi[index_bj2]; 
      }
      if(gen_index1_matched_reco_photon>=0){
            gen_matched_LeadingBjet_pt = GenPart_pt[gen_index1_matched_reco_bjet];
            gen_matched_LeadingBjet_eta = GenPart_eta[gen_index1_matched_reco_bjet];
            gen_matched_LeadingBjet_phi = GenPart_phi[gen_index1_matched_reco_bjet];
      }
      if(gen_index2_matched_reco_photon>=0){
            gen_matched_subLeadingBjet_pt = GenPart_pt[gen_index2_matched_reco_bjet];
            gen_matched_subLeadingBjet_eta = GenPart_eta[gen_index2_matched_reco_bjet];
            gen_matched_subLeadingBjet_phi = GenPart_phi[gen_index2_matched_reco_bjet];           
      } 
   //}
   tree->Fill(); 
   } // supposed to be here or above this line? 
   for(int i=0; i<nHpTbin; i++){
       cout <<"H pT bin "<<HpT_bounds[i]<<" - "<<HpT_bounds[i+1]<<" eff: "<<pass_events[i]/all_events[i] <<" pass: "<<pass_events[i]<<" all: "<<all_events[i]<<endl;   
    }
}
