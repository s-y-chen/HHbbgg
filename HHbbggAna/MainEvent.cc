//////////////////////////////////////////////////////////
// 
// Author : Nan Lu
//                   Caltech
// Date Created : March, 2021
//////////////////////////////////////////////////////////
#define MainEvent_cxx

#include "MainEvent.h"
#include <TH2.h>
#include <TStyle.h>
#include<iostream>

double MainEvent::DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double MainEvent::DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

bool MainEvent::jet_puid_sel(TString jet_puid, int puId, float pt){
  bool pass_jet_puid = false;
  if(jet_puid == "loose")
        pass_jet_puid = (pt > 50.) || (pt < 50. && puId >= 4);
  else if(jet_puid == "medium")
        pass_jet_puid = (pt > 50.) || (pt < 50. && puId >= 6);
  else if(jet_puid == "tight")
        pass_jet_puid = (pt > 50.) || (pt < 50. && puId >= 7);
  else if(jet_puid == "none")
        pass_jet_puid = true;
  return pass_jet_puid;
}