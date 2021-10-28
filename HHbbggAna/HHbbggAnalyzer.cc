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
    string samplename         = argv[3];
    const char *isData        = argv[4];
    const char *isRunGen      = argv[5];
    TString year          = argv[6];
    HHbbggAnalyzer HHbbgg(inputFileList, outFileName, samplename, isData, year);
    cout << "samplename " << samplename << " year " <<year<< endl;
    HHbbgg.EventLoop(samplename, isData, isRunGen);
    //HHbbgg.cal_sumOfgw(samplename, isData);
    return 0;
}

void HHbbggAnalyzer::cal_sumOfgw(string samplename, const char *isData)
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    cout <<"total entries: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    double sum = 0.;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%10000==0) cout <<"entry: "<<jentry<<endl;
        sum = sum + genWeight;
    }
    h_sumOfgw->SetBinContent(1,sum);            
}
    
void HHbbggAnalyzer::EventLoop(string samplename, const char *isData, const char *isRunGen)
{ 
    if (fChain == 0) return;
    clearTreeVectors();
    cout<<"cleared tree vectors\n";
    BookTreeBranches();
    cout<<"booked tree branches\n";
    float muon_mass = 0.1056583745;

    float lumi = luminosity[yearst]; //fb-1
    cout <<"year "<<yearst<<" lumi "<<lumi<<endl;
      
    bool datafile = true;
    if(*isData=='F') datafile = false;

    bool debug = false;
    
    const int nHpTbin = 10;
    float HpT_bounds[nHpTbin+1] = {0, 100, 200, 300, 400, 500, 600, 700, 1000, 2000, 3000};
    float all_events[nHpTbin] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float pass_events[nHpTbin] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    /* 
    //btag SF
    BTagCalibration calib("deepcsv","data/btagSF/DeepCSV_94XSF_V3_B_F.csv");
    BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,  // operating point
			       "central",             // central sys type
			       {"up", "down"});      // other sys types
    reader.load(calib,                // calibration instance
	      BTagEntry::FLAV_B,    // btag flavour
	      "comb")       ;        // measurement type
    */
    
    //xs
    double sumOfgenweight = 1.0;
    if(!datafile){
        sumOfgenweight = sumOfgenw[yearst][samplename];
        cout <<"samplename: "<<samplename<< " xs: "<<xs[samplename]<<" sumOfgenweight "<<sumOfgenweight<<endl;
    }
    bool skip = true;    
    bool checkGenWeight = false;
    if(samplename.find("GluGluToHHTo2B2G_node_cHHH1") != std::string::npos){
        skip = false;
        checkGenWeight = true; 
    }

    //event loop
    Long64_t nentries = fChain->GetEntriesFast();
    cout <<"total entries: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    //debug
    //for (Long64_t jentry=0; jentry<10; jentry++){
    for (Long64_t jentry=0; jentry<nentries;jentry++){
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        if(jentry%10000==0) cout <<"entry: "<<jentry<<endl;
        clearTreeVectors();
       
        genweight = 1.;
       
        //gen information: 
        //pdg ID scheme: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
        //0<= Jet_genJetIdx < nGenJet.
        int h1_index_tmp = -999, h2_index_tmp = -999;
        int photon1_index_tmp = -999, photon2_index_tmp = -999;
        int bjet1_index_tmp = -999, bjet2_index_tmp = -999, genjet1_index_tmp = -999, genjet2_index_tmp = -999;
       
        vector<int> index_higgs, index_diphoton, index_dibjet, index_digenjet;
        index_higgs.clear();
        index_diphoton.clear();
        index_dibjet.clear();
        index_digenjet.clear();
          
        //start of reco selection
        
        //HLT cut 
        //2016 : HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v* 2017 and 2018: HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*
        //see AN2017_286 Section 3.1 Trigger Selection
        if( yearst=="2016" && HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 ) trig_decision = 1;
        if( yearst=="2017" && HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 ) trig_decision = 1; 
        if( yearst=="2018" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision = 1;
        //if(trig_decision<1) cout <<"trig_decision "<<trig_decision<<endl;

        
        //photon selection
        int index_ph1(-999), index_ph2(-999);
        //match reco photon to genPhoton
        int gen_index1_matched_reco_photon = -999, gen_index2_matched_reco_photon = -999;
       
        //declare jet selection variables as well
        int index_bj1(-999), index_bj2(-999); 
        int index_vbfj1(-999), index_vbfj2(-999);
        //match reco bjet to genbjet
        int gen_index1_matched_reco_bjet = -999, gen_index2_matched_reco_bjet = -999;
        
        vector<int> index_bjet;
        index_bjet.clear();
        
        vector<int> index_vbfjet_0;
        index_vbfjet_0.clear();
        
        vector<int> index_vbfjet_1;
        index_vbfjet_1.clear();
        
        vector<int> index_photon;
        index_photon.clear();
        for(int i=0;i<nPhoton;i++){
            if(Photon_mvaID_WP90[i]==1 && Photon_pt[i] > 25){ 
                bool eta_cut = (fabs(Photon_eta[i]) < 2.5) && ((fabs(Photon_eta[i])<1.44 || fabs(Photon_eta[i])>1.57)); 
                if(eta_cut) index_photon.push_back(i);
            } 
        }
        nphoton = index_photon.size();
        if( (nphoton>1) && (trig_decision==1) ){ 
                TLorentzVector photon_1, photon_2, diphoton;
                photon_1.SetPtEtaPhiM(Photon_pt[0],Photon_eta[0],Photon_phi[0],0);
                photon_2.SetPtEtaPhiM(Photon_pt[1],Photon_eta[1],Photon_phi[1],0);
                diphoton = photon_1 + photon_2;
                    index_ph1 = index_photon.at(0);
                    index_ph2 = index_photon.at(1);
                    diphoton_pt = diphoton.Pt();
                    diphoton_mass = diphoton.M();
                    diphoton_eta = diphoton.Eta();
                    diphoton_pt_over_diphoton_mass = diphoton_pt / diphoton_mass;
                    photon_delR = DeltaR(photon_1.Eta(), photon_1.Phi(), photon_2.Eta(), photon_2.Phi());
             
            }
            
           
          
        
            
            for(int i=0;i<nJet;i++){
            //medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation medium is 1% misidentify rate for light-flavor (udsg) jet with b-jet efficiency 75%
                //if(Jet_btagDeepB[i]>= 0.4184 && (Jet_pt[i] > 25)){  
                if(Jet_btagDeepB[i]>= 0 && Jet_pt[i] > 25){  
                    index_vbfjet_0.push_back(i);
                    }
                }
            nbjet = index_bjet.size();
            
            if(nbjet>1){
                TLorentzVector bjet_1, bjet_2, dibjet;
                TLorentzVector bjet_1_corr, bjet_2_corr, dibjet_corr;
                
                    bjet_1.SetPtEtaPhiM(Jet_pt[index_bjet[0]],Jet_eta[index_bjet[0]],Jet_phi[index_bjet[0]],Jet_mass[index_bjet[0]]);
                    bjet_2.SetPtEtaPhiM(Jet_pt[index_bjet[1]],Jet_eta[index_bjet[1]],Jet_phi[index_bjet[1]],Jet_mass[index_bjet[1]]);
                    dibjet = bjet_1 + bjet_2;

                    // bjet energy regression corrections
                    bjet_1_corr.SetPtEtaPhiM(Jet_pt[index_bjet[0]] * Jet_bRegCorr[index_bjet[0]], Jet_eta[index_bjet[0]],Jet_phi[index_bjet[0]],Jet_mass[index_bjet[0]]);
                    bjet_2_corr.SetPtEtaPhiM(Jet_pt[index_bjet[1]] * Jet_bRegCorr[index_bjet[1]], Jet_eta[index_bjet[1]],Jet_phi[index_bjet[1]],Jet_mass[index_bjet[1]]);
                    dibjet_corr = bjet_1_corr + bjet_2_corr;
                    
                    index_bj1 = index_bjet[0];
                    index_bj2 = index_bjet[1];
                  
                
                          

              
                    dibjet_pt = dibjet.Pt();
                    dibjet_mass = dibjet.M();
                    dibjet_eta = dibjet.Eta();
                    dibjet_pt_over_dibjet_mass = dibjet_pt / dibjet_mass;
                    
                    dibjet_pt_corr = dibjet_corr.Pt();
                    dibjet_mass_corr = dibjet_corr.M();
                    dibjet_eta_corr = dibjet_corr.Eta();
                    dibjet_pt_over_dibjet_mass_corr = dibjet_pt_corr / dibjet_mass_corr;
                    bjet_delR = DeltaR(bjet_1.Eta(), bjet_1.Phi(), bjet_2.Eta(), bjet_2.Phi());
                    sumDeepBscore = Jet_btagDeepB[index_bj1] + Jet_btagDeepB[index_bj2];
                    
           
            }//end of at least two b-jet loop
           
            
        
        if(index_ph1>=0){
          leading_photon_pt=Photon_pt[index_ph1];
          leading_photon_eta=Photon_eta[index_ph1];
          leading_photon_phi=Photon_phi[index_ph1];
          leading_pho_pt_over_dimass = leading_photon_pt / diphoton_mass;
        }
        if(index_ph2>=0){
          subleading_photon_pt=Photon_pt[index_ph2];
          subleading_photon_eta=Photon_eta[index_ph2];
          subleading_photon_phi=Photon_phi[index_ph2]; 
          subleading_pho_pt_over_dimass = subleading_photon_pt / diphoton_mass;
        }
        // bjet
        if(index_bj1>=0){
            leading_bjet_pt=Jet_pt[index_bj1];
            leading_bjet_pt_corr=Jet_pt[index_bj1] * Jet_bRegCorr[index_bj1];
            leading_bjet_eta=Jet_eta[index_bj1];
            leading_bjet_phi=Jet_phi[index_bj1];
            leadingDeepBscore = Jet_btagDeepB[index_bj1];
            leading_bjet_pt_over_dimass = leading_bjet_pt / dibjet_mass;
            leading_bjet_pt_over_dimass_corr = leading_bjet_pt_corr / dibjet_mass_corr;
            dphi_met_leading_bjet = DeltaPhi(leading_bjet_phi, MET_phi);
        }
        if(index_bj2>=0){
            subleading_bjet_pt=Jet_pt[index_bj2];
            subleading_bjet_pt_corr=Jet_pt[index_bj2] * Jet_bRegCorr[index_bj2];
            subleading_bjet_eta=Jet_eta[index_bj2];
            subleading_bjet_phi=Jet_phi[index_bj2]; 
            subleadingDeepBscore = Jet_btagDeepB[index_bj2];
            subleading_bjet_pt_over_dimass = subleading_bjet_pt / dibjet_mass;
            subleading_bjet_pt_over_dimass_corr = subleading_bjet_pt_corr / dibjet_mass_corr;
            dphi_met_subleading_bjet = DeltaPhi(subleading_bjet_phi, MET_phi);
        }

            
        if(leading_photon_pt > 0){
            photon_recon = 1;
            
        }
        if(leading_bjet_pt > 0){
            bjet_recon = 1;
        }
        if (photon_recon=1){
            tree->Fill();
        }
        else if (bjet_recon=1){
            tree->Fill();
        }
        
    }//end of event loop
    for(int i=0; i<nHpTbin; i++){
        cout <<"H pT bin "<<HpT_bounds[i]<<" - "<<HpT_bounds[i+1]<<" eff: "<<pass_events[i]/all_events[i] <<" pass: "<<pass_events[i]<<" all: "<<all_events[i]<<endl;   
    }
}
