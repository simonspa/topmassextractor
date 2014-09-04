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

Double_t extractorYield::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << data;
    return data;
  }
  else {
    // Calculate the signal fraction from reconstructed events and TT background:
    Double_t fsignal = reco/(reco+ttbgr);

    // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
    Double_t bgr_other = bgr - ttbgr;

    // Calculate signal by subtracting backround from data, multiplied by signal fraction.
    Double_t signal = (data - bgr_other)*fsignal;

    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

  return signal;
  }
}


Double_t extractorYield::getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Scale the reco according to the different TTBar Cross sections (mass dependent):
  Double_t corr_reco = reco*getTtbarXsec(mass)/getTtbarXsec(nominalmass);
  LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    // Return the corrected reco event count plus the background:
    return (corr_reco + bgr);
  }
  else {
    // Return the reco event count corrected by the ttbar Xsec at given mass:
    return corr_reco;
  }
}

void extractorYieldBackground::prepareScaleFactor(TString systematic) {

  if(systematic.Contains("BG") || systematic.Contains("DY")) {

    if(systematic.Contains("UP")) { scaleFactor = 1.3; }
    if(systematic.Contains("DOWN")) { scaleFactor = 0.7; }
  }
}

Double_t extractorYieldBackground::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << data;
    return data;
  }
  else {
    // Calculate the signal fraction from reconstructed events and TT background:
    Double_t fsignal = reco/(reco+ttbgr);

    // Calculate the "others background" as difference between bgr and ttbgr (no xsec scaling)
    Double_t bgr_other = bgr - ttbgr;

    // Scale with DY or BG scale factor:
    bgr_other *= scaleFactor;

    // Calculate signal by subtracting backround from data, multiplied by signal fraction.
    Double_t signal = (data - bgr_other)*fsignal;

    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

    return signal;
  }
}


Double_t extractorYieldOtherSamples::getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  //variation for every bin:
  reco -= deltaRec.at(bin-1);

  // Call parent class signal calculation function:
  return extractorYield::getReco(bin, mass, reco, bgr, ttbgr);
}

Double_t extractorYieldOtherSamples::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  // variation for every mass sample in the current bin:
  reco -= deltaRec.at(bin-1);
  bgr -= deltaBgr.at(bin-1);
  ttbgr -= deltaTtbgr.at(bin-1);

  // Call parent class signal calculation function:
  return extractorYield::getSignal(bin, mass, data, reco, bgr, ttbgr);
}

void extractorYieldOtherSamples::calcDifferenceToNominal(TString nominal, TString systematic) {

  // Input files:
  TString nfilename, sfilename;
  bool splitDifference = true;

  if(systematic.Contains("HAD")) {
    nfilename = "preunfolded/POWHEG/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/MCATNLO/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else if(systematic.Contains("CR")) {
    nfilename = "preunfolded/PERUGIA11/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/PERUGIA11NoCR/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else if(systematic.Contains("UE")) {
    nfilename = "preunfolded/PERUGIA11/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/PERUGIA11mpiHi/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else {
    splitDifference = false;
    nfilename = "preunfolded/" + nominal + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/" + systematic + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }

  TFile * nominalfile = TFile::Open(nfilename,"read");
  TFile * systematicfile = TFile::Open(sfilename,"read");

  if(!nominalfile->IsOpen()) {
    LOG(logINFO) << "Failed to access file " << nfilename;
    throw 1;
  }
  LOG(logDEBUG2) << "Opened " << nfilename;

  if(!systematicfile->IsOpen()) {
    LOG(logINFO) << "Failed to access file " << sfilename;
    throw 1;
  }
  LOG(logDEBUG2) << "Opened " << sfilename;

  LOG(logDEBUG) << "Difference: " << nominal << " to " << systematic;

  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  TH1D * nominalReco = static_cast<TH1D*>(nominalfile->Get("aRecHist"));
  TH1D * varReco = static_cast<TH1D*>(systematicfile->Get("aRecHist"));

  TH1D * nominalBgr = static_cast<TH1D*>(nominalfile->Get("aBgrHist"));
  TH1D * varBgr = static_cast<TH1D*>(systematicfile->Get("aBgrHist"));

  TH1D * nominalTtbgr = static_cast<TH1D*>(nominalfile->Get("aTtBgrHist"));
  TH1D * varTtbgr = static_cast<TH1D*>(systematicfile->Get("aTtBgrHist"));

  for(Int_t bin = 1; bin <= nominalReco->GetNbinsX(); bin++) {
    Double_t rec = nominalReco->GetBinContent(bin) - varReco->GetBinContent(bin);
    Double_t bgr = nominalBgr->GetBinContent(bin) - varBgr->GetBinContent(bin);
    Double_t ttbgr = nominalTtbgr->GetBinContent(bin) - varTtbgr->GetBinContent(bin);

    if(splitDifference) { 
      rec = TMath::Abs(rec)/2;
      bgr = TMath::Abs(bgr)/2;
      ttbgr = TMath::Abs(ttbgr)/2;
      if(systematic.Contains("UP")) { rec *= -1; bgr *= -1; ttbgr *= -1; }
    }

    deltaRec.push_back(rec);
    deltaBgr.push_back(bgr);
    deltaTtbgr.push_back(ttbgr);

    LOG(logDEBUG3) << "Diff bin #" << bin << " reco: " << nominalReco->GetBinContent(bin) << " - " << varReco->GetBinContent(bin) << " = " << rec << " dbgr=" << bgr << " dttbgr=" << ttbgr;
  }

  delete nominalfile;
  delete systematicfile;
}
