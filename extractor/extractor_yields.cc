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

Double_t extractorYield::getPseudoData(Int_t /*bin*/, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

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
  std::vector<Double_t> Xbins;
  for (Int_t bin = startbin; bin <= nbins; bin++) { Xbins.push_back(aRecHist->GetBinLowEdge(bin)); }
  Xbins.push_back(aRecHist->GetBinLowEdge(nbins) + aRecHist->GetBinWidth(nbins));

  TH1D * signalHist = new TH1D(type + m_channel + Form("_m%3.1f",mass),
			       type + m_channel + Form("_m%3.1f",mass),
			       Xbins.size()-1,
			       &Xbins.front());

  // Store the binning for global use:
  bin_boundaries = Xbins;

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
  std::vector<Double_t> Xbins;
  for (Int_t bin = startbin; bin <= nbins; bin++) { Xbins.push_back(aRecHist->GetBinLowEdge(bin)); }
  Xbins.push_back(aRecHist->GetBinLowEdge(nbins) + aRecHist->GetBinWidth(nbins));

  TH1D * simulationHist = new TH1D("reco_" + m_channel + Form("_m%3.1f",mass),
				   "reco_" + m_channel + Form("_m%3.1f",mass),
				   Xbins.size()-1,
				   &Xbins.front());

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
  bool splitDifference = true;
  bool systematicUp = false;
  if(systematic.Contains("UP")) { systematicUp = true; }

  if(systematic.Contains("HAD")) { nominal = "POWHEG"; systematic = "MCATNLO"; }
  else if(systematic.Contains("CR")) { nominal = "PERUGIA11"; systematic = "PERUGIA11NoCR"; }
  else if(systematic.Contains("UE")) { nominal = "PERUGIA11"; systematic = "PERUGIA11mpiHi"; }
  else { splitDifference = false; }

  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  deltaRec = calcSampleDifference(nominal,systematic,"aRecHist");
  deltaBgr = calcSampleDifference(nominal,systematic,"aBgrHist");
  deltaTtbgr = calcSampleDifference(nominal,systematic,"aTtBgrHist");

  for(size_t i = 0; i < deltaRec.size(); ++i) {
    if(splitDifference) { 
      deltaRec.at(i) = TMath::Abs(deltaRec.at(i))/2;
      deltaBgr.at(i) = TMath::Abs(deltaBgr.at(i))/2;
      deltaTtbgr.at(i) = TMath::Abs(deltaTtbgr.at(i))/2;
      if(systematicUp) { deltaRec.at(i) *= -1; deltaBgr.at(i) *= -1; deltaTtbgr.at(i) *= -1; }
    }
    LOG(logDEBUG3) << "Diff bin #" << i+1 << " reco:  = " << deltaRec.at(i);
    LOG(logDEBUG3) << "Diff bin #" << i+1 << " bgr:  = " << deltaBgr.at(i);
    LOG(logDEBUG3) << "Diff bin #" << i+1 << " ttbgr:  = " << deltaTtbgr.at(i);
  }
}
