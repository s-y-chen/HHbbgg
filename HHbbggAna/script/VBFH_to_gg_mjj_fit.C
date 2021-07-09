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


void VBFH_to_gg_mjj_fit()
{
    TString path = "plots/sigfit/"; 
    TString filename =  "VBFHToGG"; // change to match sample name
 
    //change this to the path of signal sample you want to work with
    TString signalfile = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/HHbbggAna/condor/output/job_4_ntuple0625v1/VBFHToGG_M125_13TeV_amcatnlo_pythia8.root";
 
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
    //RooRealVar a4("a4","a4", 4,0.1,1000);
    //RooRealVar a5("a5","a5", 5,0.1,1000);
    //RooRealVar a6("a6","a6", 6,0.1,1000);
    RooAbsPdf* bern1 = new RooBernstein("bern1", "bern1", *mjj, RooArgList(a1, a2, a3));

    // Bernstein2 pdf
    RooRealVar b1("b1","b1", 7,0.1,100);
    RooRealVar b2("b2","b2", 8,0.1,100);
    RooRealVar b3("b3","b3", 9,0.1,100);
    //RooRealVar b4("b4","b4", 10,0.1,100);
    //RooRealVar b5("b5","b5", 11,0.1,100);
    //RooRealVar b6("b6","b6", 12,0.1,100);
    RooAbsPdf* bern2 = new RooBernstein("bern2", "bern2", *mjj, RooArgList(b1, b2, b3));
    
    // Bernstein3 pdf
    RooRealVar c1("c1","c1", 13,0.1,100);
    RooRealVar c2("c2","c2", 14,0.1,100);
    RooRealVar c3("c3","c3", 15,0.1,100);
    //RooRealVar c4("c4","c4", 16,0.1,100);
    //RooRealVar c5("c5","c5", 17,0.1,100);
    //RooRealVar c6("c6","c6", 18,0.1,100);
    RooAbsPdf* bern3 = new RooBernstein("bern3", "bern3", *mjj, RooArgList(c1, c2, c3));
    
    RooRealVar* frac1 = new RooRealVar("frac1","fraction of bern1",0.1,0.0,1) ;
    RooRealVar* frac2 = new RooRealVar("frac2","fraction of bern2",0.1,0.0,1) ;
    RooAbsPdf* model = new RooAddPdf("model","bn1+bn2",RooArgList(*bern1, *bern2, *bern3),RooArgList(*frac1, *frac2)) ;
       
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*)sigFile.Get("tree");
    RooRealVar* evWeight = new RooRealVar("genweight","genweight",1,-1e-10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mjj);
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
  
    RooPlot* dtframe = mjj->frame(Range(mind,maxd,kTRUE),Title("GluGluToHH dijet mass"));
    data->plotOn(dtframe);
    model->plotOn(dtframe);
    model->plotOn(dtframe, Components(*bern1), LineStyle(kDashed), LineColor(kGreen));
    model->plotOn(dtframe, Components(*bern2), LineStyle(kDashed), LineColor(kRed));
    model->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* can1 = new TCanvas();
    dtframe->Draw();
    can1->SaveAs(path+filename+obs+".png");
    can1->SaveAs(path+filename+obs+".pdf");

}
