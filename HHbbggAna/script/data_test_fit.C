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


void exponential_fit(string rt_file_name, string rt_file_path, string obs_var){
    
    TString path = "plots/MC_dimass_fit/"; 
    TString filename =  rt_file_name; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/skimming" + rt_file_path + ".root";
 
    TString min = "50";
    TString max = "300";
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


void data_test_fit()
{
    gauss_fit(names_90_ggHH[i], paths_90_ggHH[i], "diphoton_mass");
}