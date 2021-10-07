#include "Riostream.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

void CSV2Root() {

//    TString dir = gSystem->UnixPathName("DNN_Trees/combine_sequential_DNN/post_ms/data_result.csv");
//    dir.ReplaceAll("UScsvToRoot.C","");
//    dir.ReplaceAll("/./","/");

   TFile *f = new TFile("DNN_Trees/combine_sequential_DNN/post_ms/data_result.root","RECREATE");
   TTree *tree = new TTree("tree","data from csv file");
   tree->ReadFile("DNN_Trees/combine_sequential_DNN/post_ms/data_result_alt.csv","index/I:leading_photon_eta/D:leading_photon_phi/D:subleading_photon_eta/D:subleading_photon_phi/D:leading_bjet_pt_corr/D:leading_bjet_eta/D:leading_bjet_phi/D:subleading_bjet_pt_corr/D:subleading_bjet_eta/D:subleading_bjet_phi/D:leadingDeepBscore/D:subleadingDeepBscore/D:sumDeepBscore/D:leading_pho_pt_over_dimass/D:subleading_pho_pt_over_dimass/D:diphoton_pt_over_diphoton_mass/D:dibjet_pt_over_dibjet_mass_corr/D:rec_pho_bjet_min_dR/D:all_pho_bjet_min_dR/D:dphi_met_leading_bjet/D:dphi_met_subleading_bjet/D:MET_pt/D:MET_phi/D:MET_sumEt/D:nbjet/I:nphoton/I:njet/I:dibjet_mass_corr/D:diphoton_mass/D:event/D:genweight/D:DNN_score/F:mass_sculpt_cut/L:mass_sculpt_cut_sm/L",',');
   f->Write();
}