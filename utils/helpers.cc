#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sys/stat.h>
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
#include <TSystem.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TVirtualFitter.h>
#include <TMath.h>

#include "log.h"
#include "helpers.h"
#include "extractor.h"

using namespace unilog;
using namespace massextractor;

TFile * extractor::OpenFile(TString name, TString mode) {

  TString path;
  struct stat sb;

  // If we want to read an existing file, check first if it's there:
  if(mode == "read") {
    path = m_inputpath;
    ifstream inputFileStream(path + "/" + name);
    if(!inputFileStream.is_open()){
      LOG(logCRITICAL) << "File \"" << path << "/" << name << "\" does not exist!";
      throw 1;
    }
    inputFileStream.close();
  }
  else {
    path = m_outputpath;
    // If we want to write a file check that the folder exists:
    if (stat(path, &sb) != 0 || !S_ISDIR(sb.st_mode)) {
      gSystem->MakeDirectory(path);
    }
  }

  TFile * file = TFile::Open(path + "/" + name,mode);
  if(!file->IsOpen()) {
    LOG(logCRITICAL) << "Problem opening file \"" << path << "/" << name << "\"!";
    throw 1;    
  }
  
  return file;
}

std::vector<Double_t> massextractor::getBinningFromHistogram(TH1D * histo, Int_t startbin, Int_t nbins) {
  
  std::vector<Double_t> binning;

  if(nbins == 0) { nbins = histo->GetNbinsX(); }

  for (Int_t bin = startbin; bin <= nbins; bin++) { 
    binning.push_back(histo->GetBinLowEdge(bin));
  }
  binning.push_back(histo->GetBinLowEdge(nbins) + histo->GetBinWidth(nbins));
  
  return binning;
}

TString massextractor::getSampleLabel(TString systematic) {

  TString label = "";
  if(systematic.Contains("Nominal")) { label = "Nominal"; }
  else if(systematic.Contains("BTAG")) { label = "B-Tagging"; }
  else if(systematic.Contains("JER")) { label = "Jet Energy Resolution"; }
  else if(systematic.Contains("JES_MPF")) { label = "JES: MPF"; }
  else if(systematic.Contains("JES_INTERCAL")) { label = "JES: Intercalibration"; }
  else if(systematic.Contains("JES_UNCORR")) { label = "JES: Uncorrelated"; }
  else if(systematic.Contains("JES_BJES")) { label = "JES: bJES"; }
  else if(systematic.Contains("JES_FLAVOR")) { label = "JES: Flavor"; }
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

  if(systematic.Contains("PRED")) label += " (Theory)";

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
  // Nominal mass:
  else if(mass < (nominalmass+0.5)) {
    if(nominal) { fullSampleName = "Nominal"; }
    else { fullSampleName = sample; }
  }
  // Variations up:
  else if(mass < (nominalmass+1.5)) {
    if(nominal) { fullSampleName = "MASS_UP_1GEV"; }
    else { fullSampleName = sample + "_1POS"; }
  }
  else if(mass < (nominalmass+3.5)) {
    if(nominal) { fullSampleName = "MASS_UP_3GEV"; }
    else { fullSampleName = sample + "_3POS"; }
  }
  else if(mass < (nominalmass+6.5)) {
    if(nominal) { fullSampleName = "MASS_UP_6GEV"; }
    else { fullSampleName = sample + "_6POS"; }
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

template<class t>
bool extractor::isApprox(t a, t b, double eps) {
  if (fabs(a - b) < eps) { return true; }
  else { return false; }
}

float extractor::getTtbarXsec(float topmass, float energy, float* scaleerr, float * pdferr) {
    /*
     * all numbers following arxiv 1303.6254
     *
     */
    float mref=173.3;
    float referencexsec=0;
    float deltam=topmass-mref;


    float a1=0,a2=0;

    if(isApprox(energy,8.f,0.01)){
        a1=-1.1125;
        a2=0.070778;
        referencexsec=245.8;
        if(scaleerr)
            *scaleerr=0.034;
        if(pdferr)
            *pdferr=0.026;
    }
    else if(isApprox(energy,7.f,0.01)){
        a1=-1.24243;
        a2=0.890776;
        referencexsec=172.0;
        if(scaleerr)
            *scaleerr=0.034;
        if(pdferr)
            *pdferr=0.028;
    }

    float reldm=mref/(mref+deltam);

    float out= referencexsec* (reldm*reldm*reldm*reldm) * (1+ a1*(deltam)/mref + a2*(deltam/mref)*(deltam/mref));

    return out;
}

std::vector<std::string> &massextractor::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> massextractor::split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

Double_t massextractor::systStatErr(Double_t nominalStatErr, Double_t systStatErr) {
  return TMath::Abs(systStatErr - nominalStatErr);
}
