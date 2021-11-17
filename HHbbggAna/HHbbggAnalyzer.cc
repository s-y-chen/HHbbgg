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
       
        if(datafile) genweight = 1.;
        else{ 
            if(checkGenWeight && yearst=="2018" && (event==20522) && (run==1)){
                cout <<"large genWeight" <<event<< " " <<genWeight<<endl;
                continue;
            } 
            else if(checkGenWeight && yearst=="2017" && (fabs(genWeight)-0.03)>0.1){
                //event 35873 genWeight 0.272939
                //event 179625 genWeight -8.68992
                cout <<"large genWeight" <<event<<" "<<genWeight<<endl;
                continue;
            }
            else if(checkGenWeight && yearst=="2016" && (fabs(genWeight)-0.03)>0.1){
                //event 149762 genWeight 12.9627
                cout <<"large genWeight" <<event<<" "<<genWeight<<endl;
                continue;
            }
            genweight = genWeight*xs[samplename]*lumi/sumOfgenweight;
        }
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
          
        if(*isRunGen=='T' && (!skip)){
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
            //else{
                //cout <<"N(GEN-Photon)!=2: "<<index_diphoton.size()<<endl;
                //skip = true;
            //}
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
                genLeadingBjet_pt = GenPart_pt[bjet1_index_tmp];
                genLeadingBjet_eta = GenPart_eta[bjet1_index_tmp];
                genLeadingBjet_phi = GenPart_phi[bjet1_index_tmp];
                genLeadingBjet_mass = GenPart_mass[bjet1_index_tmp];       
                gensubLeadingBjet_pt = GenPart_pt[bjet2_index_tmp];
                gensubLeadingBjet_eta = GenPart_eta[bjet2_index_tmp];
                gensubLeadingBjet_phi = GenPart_phi[bjet2_index_tmp];
                gensubLeadingBjet_mass = GenPart_mass[bjet2_index_tmp];
                genBjet_deltaR = DeltaR(genLeadingBjet_eta, genLeadingBjet_phi, gensubLeadingBjet_eta, gensubLeadingBjet_phi);
                genLeadingGenjet_pt = GenJet_pt[genjet1_index_tmp];
                genLeadingGenjet_eta = GenJet_eta[genjet1_index_tmp];
                genLeadingGenjet_phi = GenJet_phi[genjet1_index_tmp];
                genLeadingGenjet_mass = GenJet_mass[genjet1_index_tmp];
                
                gensubLeadingGenjet_pt = GenJet_pt[genjet2_index_tmp];
                gensubLeadingGenjet_eta = GenJet_eta[genjet2_index_tmp];
                gensubLeadingGenjet_phi = GenJet_phi[genjet2_index_tmp];
                gensubLeadingGenjet_mass = GenJet_mass[genjet2_index_tmp];
                
            }
        }//end of gen information proccessing
          
        //start of reco selection
        
        //HLT cut 
        //2016 : HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v* 2017 and 2018: HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*
        //see AN2017_286 Section 3.1 Trigger Selection
        if( yearst=="2016" && HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 ) trig_decision = 1;
        if( yearst=="2017" && HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 ) trig_decision = 1; 
        if( yearst=="2018" && (HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90==1 || HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95==1) ) trig_decision = 1;
        //if(trig_decision<1) cout <<"trig_decision "<<trig_decision<<endl;
        
        //Primary vertex cut
        if(PV_npvsGood > 0) pv_pass = 1;
        
        //photon selection
        int index_ph1(-999), index_ph2(-999);
        //match reco photon to genPhoton
        int gen_index1_matched_reco_photon = -999, gen_index2_matched_reco_photon = -999;
       
        //declare jet selection variables as well
        int index_bj1(-999), index_bj2(-999); 
        int index_vbfj1(-999), index_vbfj2(-999);
        //match reco bjet to genbjet
        int gen_index1_matched_reco_bjet = -999, gen_index2_matched_reco_bjet = -999;
        
        
        //AK8 fat jet
        int index_fatJet(-999);
        
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
            if(nphoton==2){
                TLorentzVector photon_1, photon_2, diphoton;
                photon_1.SetPtEtaPhiM(Photon_pt[0],Photon_eta[0],Photon_phi[0],0);
                photon_2.SetPtEtaPhiM(Photon_pt[1],Photon_eta[1],Photon_phi[1],0);
                diphoton = photon_1 + photon_2;
                if(Photon_pt[0]/diphoton.M()>1./3. && Photon_pt[1]/diphoton.M()>1./4. && 100 < diphoton.M() && diphoton.M() < 180){
                    index_ph1 = index_photon.at(0);
                    index_ph2 = index_photon.at(1);
                    diphoton_pt = diphoton.Pt();
                    diphoton_mass = diphoton.M();
                    diphoton_eta = diphoton.Eta();
                    diphoton_pt_over_diphoton_mass = diphoton_pt / diphoton_mass;
                    photon_delR = DeltaR(photon_1.Eta(), photon_1.Phi(), photon_2.Eta(), photon_2.Phi());
                }
            }
            
            else if(nphoton>2){
                float tmp_diphoton_pt = 0.;
                for(int j=0; j<nphoton; j++){
                    for(int k=j+1; k<nphoton; k++){
                        TLorentzVector photon_1, photon_2, diphoton;
                        photon_1.SetPtEtaPhiM(Photon_pt[j],Photon_eta[j],Photon_phi[j],0);
                        photon_2.SetPtEtaPhiM(Photon_pt[k],Photon_eta[k],Photon_phi[k],0);
                        diphoton = photon_1 + photon_2;
                        if(Photon_pt[j]/diphoton.M()>1./3. && Photon_pt[k]/diphoton.M()>1./4. && tmp_diphoton_pt < diphoton.Pt() && 100 < diphoton.M() && diphoton.M() < 180){
                            tmp_diphoton_pt = diphoton.Pt();
                            index_ph1 = j;
                            index_ph2 = k;
                            diphoton_pt = tmp_diphoton_pt;
                            diphoton_mass = diphoton.M();
                            diphoton_eta = diphoton.Eta();
                            diphoton_pt_over_diphoton_mass = diphoton_pt / diphoton_mass;
                            photon_delR = DeltaR(photon_1.Eta(), photon_1.Phi(), photon_2.Eta(), photon_2.Phi());
                        }
                    }
                }
            }
            //the first index is gen, second is reco
            float dR11 = DeltaR(GenPart_eta[photon1_index_tmp], GenPart_phi[photon1_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
            float dR12 = DeltaR(GenPart_eta[photon1_index_tmp], GenPart_phi[photon1_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
            float dR21 = DeltaR(GenPart_eta[photon2_index_tmp], GenPart_phi[photon2_index_tmp], Photon_eta[index_ph1], Photon_phi[index_ph1]);
            float dR22 = DeltaR(GenPart_eta[photon2_index_tmp], GenPart_phi[photon2_index_tmp], Photon_eta[index_ph2], Photon_phi[index_ph2]);
            //compare w/ gen is just for matching
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
                  
            for(int i=0;i<nJet;i++){
            //medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X and https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation medium is 1% misidentify rate for light-flavor (udsg) jet with b-jet efficiency 75%
                //if(Jet_btagDeepB[i]>= 0.4184 && (Jet_pt[i] > 25)){  
                if(Jet_btagDeepB[i]>= 0 && Jet_pt[i] > 25){  
                    float dRbg1 = DeltaR(Jet_eta[i], Jet_phi[i], Photon_eta[index_ph1], Photon_phi[index_ph1]);
                    float dRbg2 = DeltaR(Jet_eta[i], Jet_phi[i], Photon_eta[index_ph2], Photon_phi[index_ph2]);
                    if (dRbg1 > 0.4 && dRbg2 > 0.4){
                        bool eta_cut = (year=="2016" ? (fabs(Jet_eta[i]) < 2.4) : (fabs(Jet_eta[i]) < 2.5)); 
                        bool jet_puid_cut = jet_puid_sel("loose", Jet_puId[i], Jet_pt[i]);
                        if(eta_cut && jet_puid_cut) index_bjet.push_back(i);
                    }
                }
                if(Jet_pt[i] > 30){  
                    float dRbg1 = DeltaR(Jet_eta[i], Jet_phi[i], Photon_eta[index_ph1], Photon_phi[index_ph1]);
                    float dRbg2 = DeltaR(Jet_eta[i], Jet_phi[i], Photon_eta[index_ph2], Photon_phi[index_ph2]);
                    if (dRbg1 > 0.4 && dRbg2 > 0.4){
                        bool eta_cut = fabs(Jet_eta[i]) < 4.7; 
                        bool jet_puid_cut = jet_puid_sel("tight", Jet_puId[i], Jet_pt[i]);
                        if(eta_cut && jet_puid_cut) index_vbfjet_0.push_back(i);
                    }
                }
            }//end of nJet loop
            nbjet = index_bjet.size();
            
            if(nbjet>1){
                TLorentzVector bjet_1, bjet_2, dibjet;
                TLorentzVector bjet_1_corr, bjet_2_corr, dibjet_corr;
                
                if(nbjet==2){
                    bjet_1.SetPtEtaPhiM(Jet_pt[index_bjet[0]],Jet_eta[index_bjet[0]],Jet_phi[index_bjet[0]],Jet_mass[index_bjet[0]]);
                    bjet_2.SetPtEtaPhiM(Jet_pt[index_bjet[1]],Jet_eta[index_bjet[1]],Jet_phi[index_bjet[1]],Jet_mass[index_bjet[1]]);
                    dibjet = bjet_1 + bjet_2;

                    // bjet energy regression corrections
                    bjet_1_corr.SetPtEtaPhiM(Jet_pt[index_bjet[0]] * Jet_bRegCorr[index_bjet[0]], Jet_eta[index_bjet[0]],Jet_phi[index_bjet[0]],Jet_mass[index_bjet[0]]);
                    bjet_2_corr.SetPtEtaPhiM(Jet_pt[index_bjet[1]] * Jet_bRegCorr[index_bjet[1]], Jet_eta[index_bjet[1]],Jet_phi[index_bjet[1]],Jet_mass[index_bjet[1]]);
                    dibjet_corr = bjet_1_corr + bjet_2_corr;
                    
                    index_bj1 = index_bjet[0];
                    index_bj2 = index_bjet[1];
                    
                }//end of if(nbjet==2){}
                
                else if(nbjet>2){
                    //b-jet selection -- want highest b tagging scores
                    float b_score_sum = -999.;
                    for(int j=0; j<nbjet; j++){
                        for(int k=j+1; k<nbjet; k++){
                            float b_score_sum_jk = 0.;
                            b_score_sum_jk = Jet_btagDeepB[index_bjet[j]] + Jet_btagDeepB[index_bjet[k]];
                            if(b_score_sum_jk > b_score_sum){
                                
                                //before corr
                                bjet_1.SetPtEtaPhiM(Jet_pt[index_bjet[j]],Jet_eta[index_bjet[j]],Jet_phi[index_bjet[j]],Jet_mass[index_bjet[j]]);
                                bjet_2.SetPtEtaPhiM(Jet_pt[index_bjet[k]],Jet_eta[index_bjet[k]],Jet_phi[index_bjet[k]],Jet_mass[index_bjet[k]]);
                                dibjet = bjet_1 + bjet_2;
                                
                                //corr
                                bjet_1_corr.SetPtEtaPhiM(Jet_pt[index_bjet[j]] * Jet_bRegCorr[index_bjet[j]],Jet_eta[index_bjet[j]],Jet_phi[index_bjet[j]],Jet_mass[index_bjet[j]]);
                                bjet_2_corr.SetPtEtaPhiM(Jet_pt[index_bjet[k]] * Jet_bRegCorr[index_bjet[k]] ,Jet_eta[k],Jet_phi[index_bjet[k]],Jet_mass[index_bjet[k]]);
                                dibjet_corr = bjet_1_corr + bjet_2_corr;
                               
                                index_bj1 = index_bjet[j];
                                index_bj2 = index_bjet[k];
                                b_score_sum = b_score_sum_jk;
                            }
                        }
                    }  
                }//end of else if(nbjet>2){}            

                // fill both
                if(70 < dibjet_corr.M() && dibjet_corr.M() < 190){
                    
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
                    
                    // first index if gen bjet, second is recon bjet; the goal is to match gen bjet and reco-bjet
                    float dR11bb = DeltaR(GenJet_eta[genjet1_index_tmp], GenJet_phi[genjet1_index_tmp], Jet_eta[index_bj1], Jet_phi[index_bj1]);
                    float dR12bb = DeltaR(GenJet_eta[genjet1_index_tmp], GenJet_phi[genjet1_index_tmp], Jet_eta[index_bj2], Jet_phi[index_bj2]);
                    float dR21bb = DeltaR(GenJet_eta[genjet2_index_tmp], GenJet_phi[genjet2_index_tmp], Jet_eta[index_bj1], Jet_phi[index_bj1]);
                    float dR22bb = DeltaR(GenJet_eta[genjet2_index_tmp], GenJet_phi[genjet2_index_tmp], Jet_eta[index_bj2], Jet_phi[index_bj2]);
       
                    if(dR11bb<dR21bb && dR11bb<0.4){
                        gen_index1_matched_reco_bjet = genjet1_index_tmp;
                        if(dR22bb<0.4) gen_index2_matched_reco_bjet = genjet2_index_tmp;
                    } 
                    else if(dR11bb>dR21bb && dR21bb<0.4){
                        gen_index1_matched_reco_bjet = genjet2_index_tmp;
                        if(dR12bb<0.4) gen_index2_matched_reco_bjet = genjet1_index_tmp;
                    }   
                }//end of di-bjet mass window if  
                else{
                    index_bj1 = -800;
                    index_bj2 = -800;
                    for(int ifatj=0; ifatj<nFatJet; ifatj++){
                        if( FatJet_btagDDBvL[ifatj]>0.5 && FatJet_pt[ifatj] > 200 && FatJet_msoftdrop[ifatj]<190. && FatJet_msoftdrop[ifatj]>70.) boostedCat = 1.0;
                    }
                }
            }//end of at least two b-jet loop
            float tmp_delR_pho_bjet = 0.; // find minimum delR between photon and jet
            for (int j=0; j < nbjet; j++){
                for (int k=0; k < nphoton; k++){
                    tmp_delR_pho_bjet = DeltaR(Jet_eta[index_bjet[j]],Jet_phi[index_bjet[j]], Photon_eta[index_photon[k]], Photon_phi[index_photon[k]]);
                    if (tmp_delR_pho_bjet < all_pho_bjet_min_dR){
                        all_pho_bjet_min_dR = tmp_delR_pho_bjet;
                    }
                }
            }
            
            if(!index_vbfjet_0.empty()){
                // remove b-jets from vbf-jet vector
                if(find(index_vbfjet_0.begin(), index_vbfjet_0.end(), index_bj1) != index_vbfjet_0.end()) {
                    index_vbfjet_0.erase(find(index_vbfjet_0.begin(), index_vbfjet_0.end(), index_bj1));
                }
                if(find(index_vbfjet_0.begin(), index_vbfjet_0.end(), index_bj2) != index_vbfjet_0.end()) {
                    index_vbfjet_0.erase(find(index_vbfjet_0.begin(), index_vbfjet_0.end(), index_bj2));
                }
            }
            nvbfjet = index_vbfjet_0.size();
            if(nvbfjet > 1){
                // check delR separation with b-jets
                for (int j=0; j<nvbfjet; j++){
                    float dRvb1 = DeltaR(Jet_eta[index_vbfjet_0[j]], Jet_phi[index_vbfjet_0[j]], Jet_eta[index_bj1], Jet_phi[index_bj1]);
                    if (dRvb1 > 0.4){
                        index_vbfjet_1.push_back(index_vbfjet_0[j]);
                    } 
                }
                nvbfjet = index_vbfjet_1.size();
                TLorentzVector vbfjet_1, vbfjet_2, divbfjet;
                if(nvbfjet == 2){
                    if(Jet_pt[index_vbfjet_1.at(0)] > 40 && Jet_pt[index_vbfjet_1.at(1)] > 30){
                        index_vbfj1 = index_vbfjet_1.at(0);
                        index_vbfj2 = index_vbfjet_1.at(1);
                    }
                } // end of nvbfjet == 2 loop
                else{
                    float divbjet_mass_tmp = -999.;
                     for(int j=0; j<nvbfjet; j++){
                        for(int k=j+1; k<nvbfjet; k++){
                            vbfjet_1.SetPtEtaPhiM(Jet_pt[index_vbfjet_1[j]], Jet_eta[index_vbfjet_1[j]], Jet_phi[index_vbfjet_1[j]], Jet_mass[index_vbfjet_1[j]]);
                            vbfjet_2.SetPtEtaPhiM(Jet_pt[index_vbfjet_1[k]], Jet_eta[index_vbfjet_1[k]], Jet_phi[index_vbfjet_1[k]], Jet_mass[index_vbfjet_1[k]]);
                            divbfjet = vbfjet_1 + vbfjet_2;
                            if(divbfjet.M() > divbjet_mass_tmp){
                                if(vbfjet_1.Pt() > 40 && vbfjet_2.Pt() > 30){
                                    divbjet_mass_tmp = divbfjet.M();
                                    index_vbfj1 = index_vbfjet_1[j];
                                    index_vbfj2 = index_vbfjet_1[k];
                                }
                            }
                        }
                    }
                } // end of nvbfjet > 2 loop 
            } // end of at least 2 vbf jet loop
        }//this is the end of the condition: if(index_photon.size()>1 && trig_decision)  
        
        //trigger if statement
        t_run =run;
        t_luminosityBlock=luminosityBlock;
        t_event=event;
        t_njet = nJet;
        t_MET_pt = MET_pt;
        t_MET_phi = MET_phi;
        t_MET_sumEt = MET_sumEt;
        
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
        if(gen_index1_matched_reco_bjet>=0){
            gen_matched_LeadingBjet_pt = GenJet_pt[gen_index1_matched_reco_bjet];
            gen_matched_LeadingBjet_eta = GenJet_eta[gen_index1_matched_reco_bjet];
            gen_matched_LeadingBjet_phi = GenJet_phi[gen_index1_matched_reco_bjet];
        }
        if(gen_index2_matched_reco_bjet>=0){
            gen_matched_subLeadingBjet_pt = GenJet_pt[gen_index2_matched_reco_bjet];
            gen_matched_subLeadingBjet_eta = GenJet_eta[gen_index2_matched_reco_bjet];
            gen_matched_subLeadingBjet_phi = GenJet_phi[gen_index2_matched_reco_bjet];           
        } 
        
        // matched dibjet variables
        if(gen_index1_matched_reco_bjet>=0 && gen_index2_matched_reco_bjet>=0){
            TLorentzVector gen_bjet_1, gen_bjet_2, gen_dibjet;
            gen_bjet_1.SetPtEtaPhiM(GenPart_pt[gen_index1_matched_reco_bjet],GenPart_eta[gen_index1_matched_reco_bjet],GenPart_phi[gen_index1_matched_reco_bjet],0);
            gen_bjet_2.SetPtEtaPhiM(GenPart_pt[gen_index2_matched_reco_bjet],GenPart_eta[gen_index2_matched_reco_bjet],GenPart_phi[gen_index2_matched_reco_bjet],0);
            gen_dibjet = gen_bjet_1 + gen_bjet_2;
            gen_dibjet_mass = gen_dibjet.M();
            gen_dibjet_pt = gen_dibjet.Pt();
            gen_dibjet_eta = gen_dibjet.Eta();
        }
        
        // vbf jet
         if(index_vbfj1>=0){
            leading_vbfjet_pt=Jet_pt[index_vbfj1];
            leading_vbfjet_eta=Jet_eta[index_vbfj1];
            leading_vbfjet_phi=Jet_phi[index_vbfj1];
            subleading_vbfjet_pt=Jet_pt[index_vbfj2];
            subleading_vbfjet_eta=Jet_eta[index_vbfj2];
            subleading_vbfjet_phi=Jet_phi[index_vbfj2]; 
            TLorentzVector vbfjet_1, vbfjet_2, divbfjet;
            vbfjet_1.SetPtEtaPhiM(leading_vbfjet_pt, leading_vbfjet_eta, leading_vbfjet_phi, Jet_mass[index_vbfj1]);
            vbfjet_2.SetPtEtaPhiM(subleading_vbfjet_pt, subleading_vbfjet_eta, subleading_vbfjet_phi, Jet_mass[index_vbfj2]);
            divbfjet = vbfjet_1 + vbfjet_2;
            divbfjet_mass = divbfjet.M();
            divbfjet_pt = divbfjet.Pt();
            divbfjet_eta = divbfjet.Eta();
            leading_vbfjet_pt_over_dimass = leading_vbfjet_pt / divbfjet_mass;
            subleading_vbfjet_pt_over_dimass = subleading_vbfjet_pt / divbfjet_mass;
            divbfjet_pt_over_dimass = divbfjet_pt / divbfjet_mass;
            vbfjet_delR = DeltaR(leading_vbfjet_eta, leading_vbfjet_phi, subleading_vbfjet_eta, subleading_vbfjet_phi);
            vbfjet_del_eta = fabs(vbfjet_1.Eta() - vbfjet_2.Eta());
        }
        
        if (boostedCat > 0){
            fatJetPt = FatJet_pt[index_fatJet];
            fatJetEta = FatJet_eta[index_fatJet];
            fatJetPhi = FatJet_phi[index_fatJet];
            fatJetMassSD_UnCorrected = FatJet_msoftdrop[index_fatJet];
            fatJetbtagDDBvL = FatJet_btagDDBvL[index_fatJet];
        }
        if (leading_photon_pt > 0 && leading_bjet_pt > 0){
            float dRbg11 = DeltaR(leading_bjet_eta, leading_bjet_phi, leading_photon_eta, leading_photon_phi);
            float dRbg12 = DeltaR(leading_bjet_eta, leading_bjet_phi, subleading_photon_eta, subleading_photon_phi);
            float dRbg21 = DeltaR(subleading_bjet_eta, subleading_bjet_phi, leading_photon_eta, leading_photon_phi);
            float dRbg22 = DeltaR(subleading_bjet_eta, subleading_bjet_phi, subleading_photon_eta, subleading_photon_phi);
            rec_pho_bjet_min_dR = min({dRbg11, dRbg12, dRbg21, dRbg22});
            TLorentzVector bjet_1, bjet_2, photon_1, photon_2, four_body;
            bjet_1.SetPtEtaPhiM(leading_bjet_pt,leading_bjet_eta,leading_bjet_phi,Jet_mass[index_bj1]);
            bjet_2.SetPtEtaPhiM(subleading_bjet_pt,subleading_bjet_eta,subleading_bjet_phi,Jet_mass[index_bj2]);
            photon_1.SetPtEtaPhiM(leading_photon_pt,leading_photon_eta,leading_photon_phi, 0);
            photon_2.SetPtEtaPhiM(subleading_photon_pt,subleading_photon_eta,subleading_photon_phi, 0);
            four_body = bjet_1 + bjet_2 + photon_1 + photon_2;
            jjyy_mass = four_body.M();
            mass_x = jjyy_mass - dibjet_mass_corr - diphoton_mass + 250;
            
            recon = 1;
            if(leading_vbfjet_pt > 0){
                VBFHH_recon = 1;
            }
            else{
                ggHH_recon = 1;
            }
            tree->Fill();
        } 
        /*if(leading_photon_pt > 0){
            photon_recon = 1;
        }
        if(leading_bjet_pt > 0){
            bjet_recon = 1;
        }*/
        
        //tree->Fill();
    }//end of event loop
    for(int i=0; i<nHpTbin; i++){
        cout <<"H pT bin "<<HpT_bounds[i]<<" - "<<HpT_bounds[i+1]<<" eff: "<<pass_events[i]/all_events[i] <<" pass: "<<pass_events[i]<<" all: "<<all_events[i]<<endl;   
    }
}
