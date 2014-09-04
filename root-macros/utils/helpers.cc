#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <Riostream.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TDecompSVD.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include "log.h"
#include "extractor.h"

using namespace unilog;
using namespace massextractor;

TString extractor::getChannelLabel(TString ch) {
 
  TString label = "";
  if(ch =="ee") { label = "ee"; }
  if(ch =="mumu"){ label = "#mu#mu"; }
  if(ch =="emu"){ label = "e#mu"; }
  if(ch =="combined"){ label = "Dilepton Combined"; }

  return label;
}

TString extractor::getSampleLabel(TString systematic) {

  TString label = "";
  if(systematic.Contains("Nominal")) { label = "Nominal"; }
  else if(systematic.Contains("BTAG")) { label = "B-Tagging"; }
  else if(systematic.Contains("JER")) { label = "Jet Energy Resolution"; }
  else if(systematic.Contains("JES")) { label = "Jet Energy Scale"; }
  else if(systematic.Contains("PU")) { label = "Pile-Up"; }
  else if(systematic.Contains("TRIG")) { label = "Trigger Eff."; }
  else if(systematic.Contains("LEPT")) { label = "Lepton Eff."; }
  else if(systematic.Contains("BG")) { label = "Background"; }
  else if(systematic.Contains("DY")) { label = "Drell-Yan"; }
  else if(systematic.Contains("KIN")) { label = "Kinematic Reconstruction"; }
  else if(systematic.Contains("MATCH")) { label = "Matching"; }
  else if(systematic.Contains("SCALE")) { label = "Q$^2$ Scale"; }
  else if(systematic.Contains("HAD")) { label = "Model"; } //{ label = "Hadronization"; }
  else if(systematic.Contains("MASS")) { label = "Top Mass"; }
  else if(systematic.Contains("CR")) { label = "Color Reconnection"; }
  else if(systematic.Contains("UE")) { label = "Underlying Event"; }
  else if(systematic.Contains("PDF")) { label = "PDF"; }

  return label;
}

TString extractor::getSampleFromMass(TString sample, Double_t mass, bool nominal) {

  TString fullSampleName = "";

  // Variations down:
  if(mass < (nominalmass-5.5)) {
    if(nominal) { fullSampleName = "MASS_DOWN_6GEV"; }
    else { fullSampleName = sample + "_6NEG"; }
  }
  else if(mass < (nominalmass-2.5)) {
    if(nominal) { fullSampleName = "MASS_DOWN_3GEV"; }
    else { fullSampleName = sample + "_3NEG"; }
  }
  else if(mass < (nominalmass-0.5)) {
    if(nominal) { fullSampleName = "MASS_DOWN_1GEV"; }
    else { fullSampleName = sample + "_1NEG"; }
  }
  // Variations up:
  else if(mass > (nominalmass+5.5)) {
    if(nominal) { fullSampleName = "MASS_UP_1GEV"; }
    else { fullSampleName = sample + "_1POS"; }
  }
  else if(mass > (nominalmass+2.5)) {
    if(nominal) { fullSampleName = "MASS_UP_3GEV"; }
    else { fullSampleName = sample + "_3POS"; }
  }
  else if(mass > (nominalmass+0.5)) {
    if(nominal) { fullSampleName = "MASS_UP_6GEV"; }
    else { fullSampleName = sample + "_6POS"; }
  }
  // Nominal mass:
  else {
    if(nominal) { fullSampleName = "Nominal"; }
    else { fullSampleName = sample; }
  }

  return fullSampleName;
}

Double_t extractor::getMassFromSample(TString sample) {

  Double_t topmass = nominalmass;
  // The mass samples are marked with "GEV":
  if(sample.Contains("GEV")) {
    if(sample.Contains("UP")) {
      if(sample.Contains("1")) topmass += 1;
      else if(sample.Contains("3")) topmass += 3;
      else if(sample.Contains("6")) topmass += 6;
    }
    else if(sample.Contains("DOWN")) {
      if(sample.Contains("1")) topmass -= 1;
      else if(sample.Contains("3")) topmass -= 3;
      else if(sample.Contains("6")) topmass -= 6;
    }
  }
  else {
    if(sample.Contains("1POS")) topmass += 1;
    else if(sample.Contains("3POS")) topmass += 3;
    else if(sample.Contains("6POS")) topmass += 6;
    else if(sample.Contains("1NEG")) topmass -= 1;
    else if(sample.Contains("3NEG")) topmass -= 3;
    else if(sample.Contains("6NEG")) topmass -= 6;
  }

  return topmass;
}
