#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
using namespace RooFit ;
using namespace std ;


void fit_mass(string rt_file_name, string rt_file_path, string obs_var){
    TString path = "plots/sigfit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/"+rt_file_path;
 
    TString min = "110";
    TString max = "140";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = obs_var;
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;
    
    // Declare observable x
    RooRealVar* mgg = new RooRealVar(obs,obs,125,mind,maxd) ;

    //gauss1 pdf
    RooRealVar* mean1 = new RooRealVar("mean1","mean of gaussian",125,mind,maxd) ;
    RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussian",10,0.1,100) ;
    RooAbsPdf* gauss1 = new RooGaussian("gauss1","gauss1",*mgg,*mean1,*sigma1) ; 
    
    //gauss2 pdf
    RooRealVar* mean2 = new RooRealVar("mean2","mean of gaussian",120,mind,maxd) ;
    RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussian",10,0.1,100) ;
    RooAbsPdf* gauss2 = new RooGaussian("gauss2","gauss2",*mgg,*mean2,*sigma2) ; 
    
    //gauss3 pdf
    RooRealVar* mean3 = new RooRealVar("mean3","mean of gaussian",130,mind,maxd) ;
    RooRealVar* sigma3 = new RooRealVar("sigma3","width of gaussian",10,0.1,100) ;
    RooAbsPdf* gauss3 = new RooGaussian("gauss3","gauss3",*mgg,*mean3,*sigma3) ; 
 
    RooRealVar* frac1 = new RooRealVar("frac1","fraction of gauss1",0.1,0.0,1) ;
    RooRealVar* frac2 = new RooRealVar("frac2","fraction of gauss2",0.1,0.0,1) ;
    RooAbsPdf* model = new RooAddPdf("model","g1+g2",RooArgList(*gauss1,*gauss2, *gauss3),RooArgList(*frac1, *frac2)) ;
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e-10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mgg);
    obsAndWeight.add(*evWeight);

    RooDataSet* data = new RooDataSet("ggfmc","ggfmc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*sigTree),Cut(cuttree)) ;
    data->Print() ;

    RooNLLVar* nll = (RooNLLVar*)model->createNLL(*data);
  
    RooMinimizer minim(*nll);
    minim.setStrategy(1);
    minim.setPrintLevel(1);
    minim.setEps(1);
    int status = minim.minimize("Minuit2", "Migrad");
  
    /*
    RooArgSet POIs(*mean1, *sigma1, *mean2, *sigma2, *frac1, *frac2);
    RooMinuit minim(*nll);
    minim.migrad() ;
    minim.minos(POIs) ;
    */
    RooFitResult *result = minim.save("fitResult","Fit Results");
    result->Write();
  
    RooPlot* dtframe = mgg->frame(Range(mind,maxd,kTRUE),Title("mass"));
    data->plotOn(dtframe);
    model->plotOn(dtframe);
    model->plotOn(dtframe, Components(*gauss1), LineStyle(kDashed), LineColor(kGreen));
    model->plotOn(dtframe, Components(*gauss2), LineStyle(kDashed), LineColor(kRed));
    model->plotOn(dtframe, Components(*gauss3), LineStyle(kDashed), LineColor(kOrange));
    model->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* c1 = new TCanvas();
    dtframe->Draw();
    c1->SaveAs(path+filename+obs+".png");
    c1->SaveAs(path+filename+obs+".pdf");
}

void dimass_fit(){

     string names [10] = {"glu_glu_to_HH_signal", "VH_to_GG", "ttH_to_GG", "VBFH_to_GG", "glu_glu_H_to_GG", "gjet_small_pt", "gjet_big_pt", "diphoton_jets_box", "diphoton_jets_box_1B", "diphoton_jets_box_2B"};
    
    string paths [10] = {"job_1_ntuple0625v1/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.root", "job_2_ntuple0625v1/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root", "job_3_ntuple0625v1/ttHToGG_M125_TuneCP5_PSweights_13TeV-powheg-pythia8.root", "job_4_ntuple0625v1/VBFHToGG_M125_13TeV_amcatnlo_pythia8.root", "job_5_ntuple0625v1/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", "job_6_ntuple0625v1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", "job_7_ntuple0625v1/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", "job_8_ntuple0625v1/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-Sherpa.root", "job_9_ntuple0625v1/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-Sherpa.root", "job_10_ntuple0625v1/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root"};
   
    for (int i=0; i< 10; i++){
        fit_mass(names[i], paths[i], "diphoton_mass");
        fit_mass(names[i], paths[i], "dibjet_condition_corr_mass");
    }
}