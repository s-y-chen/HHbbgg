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
#include "RooCrystalBall.h"
using namespace RooFit ;
using namespace std ;


void gauss_fit(string rt_file_name, string rt_file_path, string obs_var){
    TString path = "plots/sigfit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/"+rt_file_path;
 
    TString min = "115";
    TString max = "135";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = obs_var;
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;
    
    // Declare observable x
    RooRealVar* mgg = new RooRealVar(obs,obs,125,mind,maxd) ;

    //gauss1 pdf
    RooRealVar* mean1 = new RooRealVar("mean1","mean of gaussian",125,mind,maxd) ;
    RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussian",5,0.1,50) ;
    RooAbsPdf* gauss1 = new RooGaussian("gauss1","gauss1",*mgg,*mean1,*sigma1) ; 
    
    //gauss2 pdf
    RooRealVar* mean2 = new RooRealVar("mean2","mean of gaussian",120,mind,maxd) ;
    RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussian",5,0.1,50) ;
    RooAbsPdf* gauss2 = new RooGaussian("gauss2","gauss2",*mgg,*mean2,*sigma2) ; 
    
    //gauss3 pdf
    RooRealVar* mean3 = new RooRealVar("mean3","mean of gaussian",130,mind,maxd) ;
    RooRealVar* sigma3 = new RooRealVar("sigma3","width of gaussian",5,0.1,50) ;
    RooAbsPdf* gauss3 = new RooGaussian("gauss3","gauss3",*mgg,*mean3,*sigma3) ; 
 
    RooRealVar* frac1 = new RooRealVar("frac1","fraction of gauss1",0.1,0.0,1) ;
    RooRealVar* frac2 = new RooRealVar("frac2","fraction of gauss2",0.1,0.0,1) ;
    RooAbsPdf* model = new RooAddPdf("model","g1+g2",RooArgList(*gauss1,*gauss2, *gauss3),RooArgList(*frac1, *frac2)) ;
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e10,1e10) ;

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
    c1->SaveAs(path+filename+obs+"_gauss.png");
    c1->SaveAs(path+filename+obs+"_gauss.pdf");
}


void CB_fit(string rt_file_name, string rt_file_path, string obs_var)
{
    gSystem->Load("RooCrystalBall_cxx.so");
    
    TString path = "plots/sigfit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/" +rt_file_path;
 
    TString min = "115";
    TString max = "135";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = obs_var; 
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;
    //TString cuttree = obs + " < " + max + " && " + obs + " > " + min + " && genLeadingH_pt > 400 ";

    // Declare observable x
    RooRealVar* mjj = new RooRealVar(obs,obs,125,mind,maxd) ;
    //RooRealVar* hpt = new RooRealVar("genLeadingH_pt","genLeadingH_pt",125,0,100000) ;

    // CB
    RooRealVar m0("m0","m0", 125, mind, maxd);

    RooRealVar alphaL("alphaL","alphaL", 2, 0.1,10);
    RooRealVar nL("nL","nL", 2,0.1,50);
    RooRealVar sigmaL("sigmaL","sigmaL", 5,0.1,50);

    RooRealVar alphaR("alphaR","alphaR", 2,0.1,10);
    RooRealVar nR("nR","nR", 2,0.1,50);
    RooRealVar sigmaR("sigmaR","sigmaR", 5,0.1,50);
    
    RooAbsPdf* cb1 = new RooCrystalBall("cb1", "cb1", *mjj, m0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mjj);
    //obsAndWeight.add(*hpt);
    obsAndWeight.add(*evWeight);

    RooDataSet* data = new RooDataSet("ggfmc","ggfmc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*sigTree),Cut(cuttree)) ;
    data->Print() ;

    RooNLLVar* nll = (RooNLLVar*)cb1->createNLL(*data);
  
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
    TFile f(path+"fit_result_cb.root","recreate");
    RooFitResult *result = minim.save("fitResult","Fit Results");
    result->Write();
    f.Close();
 
    RooPlot* dtframe = mjj->frame(Range(mind,maxd,kTRUE),Title("mass"));
    data->plotOn(dtframe);
    cb1->plotOn(dtframe);
    cb1->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* can1 = new TCanvas();
    dtframe->Draw();
    can1->SaveAs(path+filename+obs+"_cb.png");
    can1->SaveAs(path+filename+obs+"_cb.pdf");
}


void Bernstein_fit(string rt_file_name, string rt_file_path, string obs_var)
{
    TString path = "plots/sigfit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/"  +rt_file_path;
 
    TString min = "100";
    TString max = "150";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = obs_var;
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;
    
    // Declare observable x
    RooRealVar* mjj = new RooRealVar(obs,obs,125,mind,maxd) ;
    
    // Bernstein1 pdf
    RooRealVar a1("a1","a1", 1,0.1,100);
    RooRealVar a2("a2","a2", 2,0.1,100);
    RooRealVar a3("a3","a3", 3,0.1,100);
    RooRealVar a4("a4","a4", 4,0.1,100);
    //RooRealVar a5("a5","a5", 5,0.1,1000);
    //RooRealVar a6("a6","a6", 6,0.1,1000);
    RooAbsPdf* bern1 = new RooBernstein("bern1", "bern1", *mjj, RooArgList(a1, a2, a3, a4));
    
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mjj);
    obsAndWeight.add(*evWeight);

    RooDataSet* data = new RooDataSet("ggfmc","ggfmc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*sigTree),Cut(cuttree)) ;
    data->Print() ;

    RooNLLVar* nll = (RooNLLVar*)bern1->createNLL(*data);
  
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
  
    RooPlot* dtframe = mjj->frame(Range(mind,maxd,kTRUE),Title("mass"));
    data->plotOn(dtframe);
    bern1->plotOn(dtframe);
    bern1->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* can1 = new TCanvas();
    dtframe->Draw();
    can1->SaveAs(path+filename+obs+"_bern.png");
    can1->SaveAs(path+filename+obs+"_bern.pdf");

}



void exponential_fit(string rt_file_name, string rt_file_path, string obs_var){
    
    TString path = "plots/sigfit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/" + rt_file_path;
 
    TString min = "100";
    TString max = "150";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = obs_var; 
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;

    // Declare observable x
    RooRealVar* mjj = new RooRealVar(obs,obs,125,mind,maxd) ;
    
    // Exponential
    RooRealVar k("k", "k", -3, -100, 0.1);
    RooAbsPdf* ex1 = new RooExponential("ex1", "ex1", *mjj, k);
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mjj);
    obsAndWeight.add(*evWeight);

    RooDataSet* data = new RooDataSet("ggfmc","ggfmc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*sigTree),Cut(cuttree)) ;
    data->Print() ;

    RooNLLVar* nll = (RooNLLVar*)ex1->createNLL(*data);
  
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
  
    RooPlot* dtframe = mjj->frame(Range(mind,maxd,kTRUE),Title("mass"));
    data->plotOn(dtframe);
    ex1->plotOn(dtframe);
    ex1->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* can1 = new TCanvas();
    dtframe->Draw();
    can1->SaveAs(path+filename+obs+"_expon.png");
    can1->SaveAs(path+filename+obs+"_expon.pdf");
}

void dimass_fit(){

     string names [10] = {"glu_glu_to_HH_signal", "VH_to_GG", "ttH_to_GG", "VBFH_to_GG", "glu_glu_H_to_GG", "gjet_small_pt", "gjet_big_pt", "diphoton_jets_box", "diphoton_jets_box_1B", "diphoton_jets_box_2B"};
    
    string paths [10] = {"job_1_ntuple0625v1/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.root", "job_2_ntuple0625v1/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root", "job_3_ntuple0625v1/ttHToGG_M125_TuneCP5_PSweights_13TeV-powheg-pythia8.root", "job_4_ntuple0625v1/VBFHToGG_M125_13TeV_amcatnlo_pythia8.root", "job_5_ntuple0625v1/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", "job_6_ntuple0625v1/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", "job_7_ntuple0625v1/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root", "job_8_ntuple0625v1/DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-Sherpa.root", "job_9_ntuple0625v1/DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-Sherpa.root", "job_10_ntuple0625v1/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa.root"};
   
    for (int i=0; i<1; i++){
        //gauss_fit(names[i], paths[i], "diphoton_mass");
        //CB_fit(names[i], paths[i], "diphoton_mass");
    }
    
    //gauss_fit(names[0], paths[0], "dibjet_condition_corr_mass"); //ggHH signal
    CB_fit(names[0], paths[0], "dibjet_mass_corr"); //ggHH signal
    /*
    Bernstein_fit(names[3], paths[3], "dibjet_condition_corr_mass"); //VBFH
    Bernstein_fit(names[4], paths[4], "dibjet_condition_corr_mass"); // ggH
    gauss_fit(names[2], paths[2], "dibjet_condition_corr_mass"); // ttH
    exponential_fit(names[1], paths[1], "dibjet_condition_corr_mass"); // VH
    */
}
