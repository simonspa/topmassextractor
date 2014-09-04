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
    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << data << " (incl. bgr)";
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

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    // Return a full pseudo data set including scaled backgrounds:
    Double_t reco_data = getPseudoData(bin,mass,reco,bgr,ttbgr);
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " reco_data=" << reco_data;
    return reco_data;
  }
  else {
    // Scale the reco according to the different TTBar Cross sections (mass dependent):
    Double_t corr_reco = reco*getTtbarXsec(mass)/getTtbarXsec(nominalmass);
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

    // Return the reco event count corrected by the ttbar Xsec at given mass:
    return corr_reco;
  }
}

Double_t extractorYield::getPseudoData(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  Double_t xsecCorrection = getTtbarXsec(mass)/getTtbarXsec(nominalmass);

  // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
  Double_t bgr_other = bgr - ttbgr;
  Double_t pseudodata = (reco + ttbgr)*xsecCorrection + bgr_other;

  return pseudodata;
}

TFile * extractorYield::selectInputFile(TString sample) {
  // Input files for Total Yield mass extraction: preunfolded histograms:
  TString filename = "preunfolded/" + sample + "/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TFile * input = OpenFile(filename,"read");
  if(!input->IsOpen()) {
    LOG(logCRITICAL) << "Failed to access data file " << filename;
    throw 1;
  }
  LOG(logDEBUG) << "Successfully opened file " << filename;
  return input;
}

void extractorYield::setClosureSample(TString closure) {

  // Enable closure test:
  doClosure = true;

  // Input file:
  TString filename = "preunfolded/" + closure + "/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TFile * closureFile = OpenFile(filename,"read");

  // Histograms containing the background:
  TH1D * aRecHist = static_cast<TH1D*>(closureFile->Get("aRecHist"));
  TH1D * aTtBgrHist = static_cast<TH1D*>(closureFile->Get("aTtBgrHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(closureFile->Get("aBgrHist"));

  // Build a pseudo data set from the closure sample:
  pseudoData = new TH1D("pseudodata","pseudodata",aRecHist->GetNbinsX(),0,100);
  pseudoData->SetDirectory(0);
  Double_t mass = getMassFromSample(closure);

  LOG(logINFO) << "Running Closure test. Pseudo data taken from " << filename;
  LOG(logINFO) << "Pseudo data mass = " << mass;

  for(Int_t bin = 1; bin <= pseudoData->GetNbinsX(); bin++) {

    Double_t pdata = getPseudoData(bin,mass,aRecHist->GetBinContent(bin),aBgrHist->GetBinContent(bin),aTtBgrHist->GetBinContent(bin));

    LOG(logDEBUG) << "Closure: Bin #" << bin << " sig=" << aRecHist->GetBinContent(bin) << " pdat=" << pdata;
    // Write pseudo data with background:
    pseudoData->SetBinContent(bin,pdata);
  }

  delete closureFile;
}

TH1D * extractorYield::getSignalHistogram(Double_t mass, TFile * histos) {

  // Histogram containing data:
  TH1D * aDataHist;
  TString type = "";
  if(!doClosure) {
    aDataHist = static_cast<TH1D*>(histos->Get("aDataHist"));
    type = "data_";
  }
  else {
    aDataHist = static_cast<TH1D*>(pseudoData);
    type = "pseudodata_";
  }

  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist"));
  // Histograms containing the background:
  TH1D * aTtBgrHist = static_cast<TH1D*>(histos->Get("aTtBgrHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(histos->Get("aBgrHist"));

  // Create a new histogram with the same binning:
  Int_t nbins = aDataHist->GetNbinsX();
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Data hist has " << nbins << " bins, using " << (nbins-startbin+1);
  Double_t Xbins[nbins+2-startbin];
  for (Int_t bin = startbin; bin <= nbins; bin++) { Xbins[bin-startbin] = aRecHist->GetBinLowEdge(bin); }
  Xbins[nbins+1-startbin] = aRecHist->GetBinLowEdge(nbins) + aRecHist->GetBinWidth(nbins);
  TH1D * signalHist = new TH1D(type + m_channel + Form("_m%3.1f",mass),
			       type + m_channel + Form("_m%3.1f",mass),
			       nbins-startbin+1,
			       Xbins);

  // Store the binning for global use:
  bin_boundaries.clear();
  for(Int_t i = 0; i < nbins+2-startbin; i++) { bin_boundaries.push_back(Xbins[i]); }

  // Iterate over all bins:
  for(Int_t bin = startbin; bin <= nbins; bin++) {

    // Get signal corrected by ttbar background:
    Double_t signal = getSignal(bin, mass, aDataHist->GetBinContent(bin),
				aRecHist->GetBinContent(bin),
				aBgrHist->GetBinContent(bin),
				aTtBgrHist->GetBinContent(bin));
    
    // Write background subtrated signal:
    signalHist->SetBinContent(bin-startbin+1,signal);
    signalHist->SetBinError(bin-startbin+1,TMath::Sqrt(signal));
  }

  if((flags & FLAG_NORMALIZE_YIELD) != 0) {
    signalHist->Scale(1./signalHist->Integral("width"));
    LOG(logDEBUG) << "Normalized Data hist.";
  }

  for(Int_t b = 1; b <= signalHist->GetNbinsX(); ++b) {
    LOG(logDEBUG2) << "#" << b << ": " << signalHist->GetBinContent(b) << " +-" << signalHist->GetBinError(b);
  }

  // Return signal-only histogram:
  return signalHist;
}

TH1D * extractorYield::getSimulationHistogram(Double_t mass, TFile * histos) {

  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(histos->Get("aBgrHist"));
  TH1D * aTtBgrHist = static_cast<TH1D*>(histos->Get("aTtBgrHist"));

  // Create a new histogram with the same binning:
  Int_t nbins = aRecHist->GetNbinsX();
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Reco hist has " << nbins << " bins, using " << (nbins-startbin+1);
  Double_t Xbins[nbins+2-startbin];
  for (Int_t bin = startbin; bin <= nbins; bin++) Xbins[bin-startbin] = aRecHist->GetBinLowEdge(bin);
  Xbins[nbins+1-startbin] = aRecHist->GetBinLowEdge(nbins) + aRecHist->GetBinWidth(nbins);
  TH1D * simulationHist = new TH1D("reco_" + m_channel + Form("_m%3.1f",mass),
			       "reco_" + m_channel + Form("_m%3.1f",mass),
			       nbins-startbin+1,
			       Xbins);

  // Iterate over all bins:
  for(Int_t bin = startbin; bin <= nbins; bin++) {

    // Correct the Reco events for different TTBar Cross sections (mass dependent):
    Double_t corr_reco = getReco(bin, mass,aRecHist->GetBinContent(bin),aBgrHist->GetBinContent(bin),aTtBgrHist->GetBinContent(bin));

    // Write corrected Reco:
    simulationHist->SetBinContent(bin-startbin+1,corr_reco);
    simulationHist->SetBinError(bin-startbin+1,aRecHist->GetBinError(bin));
  }

  if((flags & FLAG_NORMALIZE_YIELD) != 0) {
    simulationHist->Scale(1./simulationHist->Integral("width"));
    LOG(logDEBUG) << "Normalized Reco hist.";
  }

  for(Int_t b = 1; b <= simulationHist->GetNbinsX(); ++b) {
    LOG(logDEBUG2) << "#" << b << ": " << simulationHist->GetBinContent(b) << " +-" << simulationHist->GetBinError(b);
  }
  
  // Return reco histogram:
  return simulationHist;
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

  // Subtract the difference in event count for the nominal mass bin and systematics variation for every bin:
  reco -= deltaRec.at(bin-1);
  bgr -= deltaBgr.at(bin-1);
  ttbgr -= deltaTtbgr.at(bin-1);

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
    nfilename = "preunfolded/POWHEG/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/MCATNLO/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else if(systematic.Contains("CR")) {
    nfilename = "preunfolded/PERUGIA11/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/PERUGIA11NoCR/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else if(systematic.Contains("UE")) {
    nfilename = "preunfolded/PERUGIA11/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/PERUGIA11mpiHi/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }
  else {
    splitDifference = false;
    nfilename = "preunfolded/" + nominal + "/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    sfilename = "preunfolded/" + systematic + "/" + m_channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  }

  TFile * nominalfile = OpenFile(nfilename,"read");
  TFile * systematicfile = OpenFile(sfilename,"read");

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

    LOG(logDEBUG3) << "Diff bin #" << bin << " reco: " << nominalReco->GetBinContent(bin) << " - " << varReco->GetBinContent(bin) << " = " << rec;
    LOG(logDEBUG3) << "Diff bin #" << bin << " bgr: " << nominalBgr->GetBinContent(bin) << " - " << varBgr->GetBinContent(bin) << " dbgr=" << bgr;
    LOG(logDEBUG3) << "Diff bin #" << bin << " ttbgr: " << nominalTtbgr->GetBinContent(bin) << " - " << varTtbgr->GetBinContent(bin) << " dbgr=" << ttbgr;
  }

  delete nominalfile;
  delete systematicfile;
}
