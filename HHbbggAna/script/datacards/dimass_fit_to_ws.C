#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "RooCrystalBall.h"

using namespace RooFit ;
using namespace std ;

RooAbsPdf* gauss_func(RooRealVar* obs, double mind, double maxd){

    //gauss1 pdf
    RooRealVar* mean1 = new RooRealVar("mean1","mean of gaussian",125,mind,maxd) ;
    RooRealVar* sigma1 = new RooRealVar("sigma1","width of gaussian",2.5,0.1,50) ;
    RooAbsPdf* gauss1 = new RooGaussian("gauss1","gauss1",*obs,*mean1,*sigma1) ; 
    
    //gauss2 pdf
    RooRealVar* mean2 = new RooRealVar("mean2","mean of gaussian",120,mind,maxd) ;
    RooRealVar* sigma2 = new RooRealVar("sigma2","width of gaussian",2.5,0.1,50) ;
    RooAbsPdf* gauss2 = new RooGaussian("gauss2","gauss2",*obs,*mean2,*sigma2) ; 
    
    //gauss3 pdf
    RooRealVar* mean3 = new RooRealVar("mean3","mean of gaussian",130,mind,maxd) ;
    RooRealVar* sigma3 = new RooRealVar("sigma3","width of gaussian",2.5,0.1,50) ;
    RooAbsPdf* gauss3 = new RooGaussian("gauss3","gauss3",*obs,*mean3,*sigma3) ; 
 
    RooRealVar* frac1 = new RooRealVar("frac1","fraction of gauss1",0.1,0.0,1) ;
    RooRealVar* frac2 = new RooRealVar("frac2","fraction of gauss2",0.1,0.0,1) ;
    RooRealVar* frac3 = new RooRealVar("frac3","fraction of gauss3",0.1,0.0,1) ;
    RooAbsPdf* model = new RooAddPdf("model","g1+g2",RooArgList(*gauss1,*gauss2, *gauss3),RooArgList(*frac1, *frac2, *frac3)) ;
    
    return model;
}

RooAbsPdf* CB_func(RooRealVar* obs, double mind, double maxd){
    gSystem->Load("RooCrystalBall_cxx.so");
    
    RooRealVar* m0 = new RooRealVar("m0","m0", 125, mind, maxd);
    RooRealVar* alphaL = new RooRealVar("alphaL","alphaL", 2, 0.1,10);
    RooRealVar* nL = new RooRealVar("nL","nL", 2,0.1,50);
    RooRealVar* sigmaL = new RooRealVar("sigmaL","sigmaL", 5,0.1,50);
    RooRealVar* alphaR = new RooRealVar("alphaR","alphaR", 2,0.1,10);
    RooRealVar* nR = new RooRealVar("nR","nR", 2,0.1,50);
    RooRealVar* sigmaR = new RooRealVar("sigmaR","sigmaR", 5,0.1,50);   
    RooAbsPdf* model = new RooCrystalBall("model", "CB", *obs, *m0, *sigmaL, *sigmaR, *alphaL, *nL, *alphaR, *nR);       
    
    return model;   
}
    
RooAbsPdf* Bernstein_func(RooRealVar* obs, double mind, double maxd){
    
    // Bernstein1 pdf
    RooRealVar* a1 = new RooRealVar("a1","a1", 0.5,0,1);
    RooRealVar* a2 = new RooRealVar("a2","a2", 0.5,0,1);
    RooRealVar* a3 = new RooRealVar("a3","a3", 0.5,0,1);
    //RooRealVar* a4 = new RooRealVar("a4","a4", 1,0,1);
    //RooRealVar a5("a5","a5", 5,0.1,1000);
    //RooRealVar a6("a6","a6", 6,0.1,1000);
    RooAbsPdf* model = new RooBernstein("model", "Bern", *obs, RooArgList(*a1, *a2, *a3));
      
    return model;   
}

RooAbsPdf* exponential_func(RooRealVar* obs, double mind, double maxd){

    RooRealVar* k = new RooRealVar("k", "k", -3, -100, 0.1);
    RooAbsPdf* model = new RooExponential("model", "Expo", *obs, *k);
     
    return model;
}

RooAbsPdf* get_func(TString func_name, RooRealVar* obs, double mind, double maxd){
    RooAbsPdf* model = NULL;
    if(func_name=="Gaussian"){
        cout <<"Guassian function is chosen"<<endl;
        model = gauss_func(obs,  mind, maxd);
    } 
    else if(func_name=="CB"){
        cout <<"CB function is chosen"<<endl; 
        model = CB_func(obs, mind, maxd);
    } 
    else if(func_name=="Bern"){
        cout <<"Bern function is chosen"<<endl; 
        model = Bernstein_func(obs, mind, maxd);
    } 
    else if(func_name=="Expo"){
        cout <<"Expo function is chosen"<<endl; 
        model = exponential_func(obs, mind, maxd);
    } 
    else{
        cout <<"func error"<<endl;
        exit(1);
    }
    return model;
}

void dofit(TString file, TString obs_var, TString min, TString max, TString dnn_var, TString dnn_cut, TString weightvar, TString procname, TString funcname, TString catname){
    TString path = "pdfs/"; 
  
    double mind = min.Atof();
    double maxd = max.Atof();

    TString cuttree = obs_var + " < " + max + " && " + obs_var + " > " + min + " && "+TString(dnn_cut);
           
    // Declare observable x
    RooRealVar* mgg = new RooRealVar(obs_var,obs_var,125,mind,maxd) ;
    RooRealVar* dnn = new RooRealVar(dnn_var,dnn_var,0.5,0,1) ;

    cout <<"before model"<<endl;
    RooAbsPdf* model = get_func(funcname,mgg, mind, maxd);
    model->Print("V");
    cout <<"after model"<<endl;

    TFile File(file);
    TTree* procTree = (TTree*)File.Get("tree");
    RooRealVar* evWeight = new RooRealVar(weightvar,weightvar,1,-1e10,1e10) ;

    RooArgSet obsAndWeight;
    obsAndWeight.add(*mgg);
    obsAndWeight.add(*dnn);
    obsAndWeight.add(*evWeight);

    RooDataSet* data = new RooDataSet("mc","mc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*procTree),Cut(cuttree)) ;
    //RooDataSet data("mc","mc",RooArgSet(obsAndWeight),RooFit::WeightVar(*evWeight),Import(*procTree),Cut(cuttree)) ;
    data->Print() ;
    //model->fitTo(data);
    
   
    cout <<"before nll"<<endl;
    RooNLLVar* nll = (RooNLLVar*)model->createNLL(*data);
    cout <<"after nll"<<endl;

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
        
    RooPlot* dtframe = mgg->frame(Range(mind,maxd,kTRUE),Title("mass"));
    data->plotOn(dtframe);
    model->plotOn(dtframe);
    //model->plotOn(dtframe, Components(*gauss1), LineStyle(kDashed), LineColor(kGreen));
    //model->plotOn(dtframe, Components(*gauss2), LineStyle(kDashed), LineColor(kRed));
    //model->plotOn(dtframe, Components(*gauss3), LineStyle(kDashed), LineColor(kOrange));
    model->paramOn(dtframe,Layout(0.55)) ;
    TCanvas* c1 = new TCanvas();
    dtframe->Draw();
    c1->SaveAs(path+procname+obs_var+"_"+funcname+catname+".png");
    c1->SaveAs(path+procname+obs_var+"_"+funcname+catname+".pdf");
    
    TString output_file = path+"wsinput."+funcname+procname+catname+".root";
    
    RooFitResult *result = minim.save("fitResult","Fit Results");
    result->Write();

    RooWorkspace *w = new RooWorkspace("pdf","workspace") ;    
    w->import(*model);
    w->Print();
    cout <<"write ws to "<<output_file<<endl;
    w->writeToFile(output_file);
}

void dimass_fit_to_ws(){
    
    TString path = "/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/notebook/ML/DNN_Trees/combine_sequential_DNN/post_ms/";
    TString ggHH_file = "sig_5.root";
    
    TString VH_file = "VHToGG_5.root";
    TString ttH_file = "ttHToGG_5.root";
    TString VBFH_file = "VBFToHH_5.root";
    TString ggH_file = "GluGluHtoGG_5.root";
    TString nonres_file = "nonres_5.root";
    TString data_file = "data_result.root"; 
    
   // category mass_sculpt_cut_sm = 1
    dofit(path+ggHH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","gghh","Gaussian","ggHHcat1"); //ggHH signal
    dofit(path+ttH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","tth","Gaussian","ggHHcat1"); //ttH bkg
    dofit(path+VH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","vh","Gaussian","ggHHcat1"); //VH bkg
    dofit(path+VBFH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","vbfh","Gaussian","ggHHcat1"); //VBFH bkg
    dofit(path+ggH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","ggh","Gaussian","ggHHcat1");
    //ggH bkg   
    dofit(path+nonres_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","nonresonant","Bern","ggHHcat1"); // nonres bkg 
    dofit(path+data_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight","data","Bern","ggHHcat1"); //data bkg
    
    dofit(path+ggHH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","gghh","CB","ggHHcat1"); //ggHH signal
    dofit(path+ttH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","tth","CB","ggHHcat1"); //ttH bkg
    dofit(path+VH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","vh","CB","ggHHcat1"); //VH bkg
    dofit(path+VBFH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","vbfh","CB","ggHHcat1"); //VBFH bkg
    dofit(path+ggH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==1","genweight_scale","ggh","CB","ggHHcat1");
    //ggH bkg   
    
    
    // category mass_sculpt_cut_sm = 0
    dofit(path+ggHH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","gghh","Gaussian","ggHHcat2"); //ggHH signal
    dofit(path+ttH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","tth","Gaussian","ggHHcat2"); //ttH bkg
    dofit(path+VH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","vh","Gaussian","ggHHcat2"); //VH bkg
    dofit(path+VBFH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","vbfh","Gaussian","ggHHcat2"); //VBFH bkg
    dofit(path+ggH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","ggh","Gaussian","ggHHcat2"); //ggH bkg
    dofit(path+data_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","nonresonant","Bern","ggHHcat2");// nonres bkg
    dofit(path+data_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight","data","Bern","ggHHcat2"); //data bkg
    
     dofit(path+ggHH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","gghh","CB","ggHHcat2"); //ggHH signal
    dofit(path+ttH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","tth","CB","ggHHcat2"); //ttH bkg
    dofit(path+VH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","vh","CB","ggHHcat2"); //VH bkg
    dofit(path+VBFH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","vbfh","CB","ggHHcat2"); //VBFH bkg
    dofit(path+ggH_file, "diphoton_mass", "110", "140", "mass_sculpt_cut_sm","mass_sculpt_cut_sm==0","genweight_scale","ggh","CB","ggHHcat2"); //ggH bkg

}
