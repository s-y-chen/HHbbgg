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


void glugluH_to_gg_mjj_fit()
{
    TString path = "plots/sigfit/"; 
    TString filename =  "GluGluHToGG"; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/job_5_ntuple0625v1/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root";
 
    TString min = "0";
    TString max = "200";
    double mind = min.Atof();
    double maxd = max.Atof();
    TString obs = "dibjet_condition_corr_mass";
    TString cuttree = obs + " < " + max + " && " + obs + " > " + min;
    
    // Declare observable x
    RooRealVar* mjj = new RooRealVar(obs,obs,125,mind,maxd) ;
    
    // Bernstein1 pdf
    RooRealVar a1("a1","a1", 1,0.1,1000);
    RooRealVar a2("a2","a2", 2,0.1,1000);
    RooRealVar a3("a3","a3", 3,0.1,1000);
    RooRealVar a4("a4","a4", 4,0.1,1000);
    RooRealVar a5("a5","a5", 5,0.1,1000);
    RooRealVar a6("a6","a6", 6,0.1,1000);
    RooAbsPdf* bern1 = new RooBernstein("bern1", "bern1", *mjj, RooArgList(a1, a2, a3, a4, a5, a6));
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e-10,1e10) ;

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
  
    RooPlot* dtframe = mjj->frame(Range(mind,maxd,kTRUE),Title("GluGluToHH dijet mass"));
    data->plotOn(dtframe);
    bern1->plotOn(dtframe);
    bern1->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* can1 = new TCanvas();
    dtframe->Draw();
    can1->SaveAs(path+filename+obs+".png");
    can1->SaveAs(path+filename+obs+".pdf");

}
