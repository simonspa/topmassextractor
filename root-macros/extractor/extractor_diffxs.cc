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

void extractorDiffXSec::setClosureSample(TString /*closure*/) {
  LOG(logERROR) << "Not possible to request closure for DiffXSec extraction. Re-run unfolding with closure flag enabled in order to get closure pseudo data.";
  LOG(logINFO) << "Use \"setUnfoldingMass()\" to select different mass samples used for unfolding.";
}

void extractorDiffXSec::setUnfoldingMass(Double_t mass) {
  LOG(logINFO) << "Requested to use m_t=" << mass << " for unfolding.";
  unfoldingMass = mass;
}

void extractorDiffXSecScaled::prepareScaleFactors(TString ch, TString systematic) {

  Int_t sign = 0;
  if(systematic.Contains("UP")) { sign = 1; }
  else if(systematic.Contains("DOWN")) { sign = -1; }

  if(systematic.Contains("PDF")) {
    LOG(logDEBUG2) << "Preparing PDF uncertainty scale factors...";
    // PDF scaling for all five bins:
    if(ch == "combined") {
      scaleFactors.push_back(sign*0.0518874555267725);
      scaleFactors.push_back(sign*0.0307283998130669);
      scaleFactors.push_back(sign*0.0124526579196978);
      scaleFactors.push_back(sign*0.00149147332885075);
      scaleFactors.push_back(sign*0.00473063248751739);
    }
    else if(ch == "ee") {
      scaleFactors.push_back(sign*0.0520053995180302);
      scaleFactors.push_back(sign*0.0303351933038199);
      scaleFactors.push_back(sign*0.0126092037398061);
      scaleFactors.push_back(sign*0.00156712414239097);
      scaleFactors.push_back(sign*0.00486566872899394);
    }
    else if(ch == "emu") {
      scaleFactors.push_back(sign*0.0525837643853782);
      scaleFactors.push_back(sign*0.0303790574991724);
      scaleFactors.push_back(sign*0.0124525455719107);
      scaleFactors.push_back(sign*0.00151246828268544);
      scaleFactors.push_back(sign*0.00463153835628149);
    }
    else if(ch == "mumu") {
      scaleFactors.push_back(sign*0.04862795023471);
      scaleFactors.push_back(sign*0.0321488889850916);
      scaleFactors.push_back(sign*0.0123548702575585);
      scaleFactors.push_back(sign*0.00139357258960364);
      scaleFactors.push_back(sign*0.00496195593202099);
    }
  }
  std::stringstream sv;
  for(std::vector<Double_t>::iterator sf = scaleFactors.begin(); sf != scaleFactors.end(); ++sf) { sv << *sf << " "; }
  LOG(logDEBUG2) << "Scale factors for every bin: " << sv.str();
}

TH1D * extractorDiffXSec::getSignalHistogram(Double_t mass, TFile * histos) {

  // Histogram containing differential cross section from data:
  TH1D * aDiffXSecHist = static_cast<TH1D*>(histos->Get("unfoldedHistNorm"));

  // Create a new histogram with the same binning:
  Int_t nbins = aDiffXSecHist->GetNbinsX();
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Data hist has " << nbins << " bins, using " << (nbins-startbin+1);
  Double_t Xbins[nbins+2-startbin];
  for (Int_t bin = startbin; bin <= nbins; bin++) Xbins[bin-startbin] = aDiffXSecHist->GetBinLowEdge(bin);
  Xbins[nbins+1-startbin] = aDiffXSecHist->GetBinLowEdge(nbins) + aDiffXSecHist->GetBinWidth(nbins);
  TH1D * signalHist = new TH1D("diffxs_" + channel + Form("_m%3.1f",mass),
			       "diffxs_" + channel + Form("_m%3.1f",mass),
			       nbins-startbin+1,
			       Xbins);

  // Store the binning for global use:
  bin_boundaries.clear();
  for(Int_t i = 0; i < nbins+2-startbin; i++) { bin_boundaries.push_back(Xbins[i]); }

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    Double_t signal = getSignal(bin-startbin,mass,aDiffXSecHist->GetBinContent(bin));
    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << aDiffXSecHist->GetBinContent(bin) << " signal=" << signal;
    signalHist->SetBinContent(bin+1-startbin,signal);
    signalHist->SetBinError(bin+1-startbin,aDiffXSecHist->GetBinError(bin));
  }

  // Return DiffXSec signal histogram:
  return signalHist;
}

Double_t extractorDiffXSec::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return data;
}

Double_t extractorDiffXSecScaled::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return (1+scaleFactors.at(bin))*data;
}

TH1D * extractorDiffXSec::getSimulationHistogram(Double_t mass, TFile * histos) {

  std::vector<TString> filenames;
  TString filename, sample;

  if(mass < 167) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_166_massdown.root"; }
  else if(mass < 170) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_169_massdown.root"; }
  else if(mass < 172) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_massdown.root"; }
  else if(mass < 173) { sample = "Nominal"; filename = "_ttbarsignalplustau.root"; }
  else if(mass < 174) { sample = "MASS_UP"; filename = "_ttbarsignalplustau_massup.root"; }
  else if(mass < 176) { sample = "MASS_UP"; filename = "_ttbarsignalplustau_175_massup.root"; }
  else { sample = "MASS_UP"; filename = "_ttbarsignalplustau_178_massup.root"; }

  // Check if we need to add up several histograms:
  if(channel == "ee" || channel == "combined") { filenames.push_back("selectionRoot/" + sample + "/ee/ee" + filename); }
  if(channel == "emu" || channel == "combined") { filenames.push_back("selectionRoot/" + sample + "/emu/emu" + filename); }
  if(channel == "mumu" || channel == "combined") { filenames.push_back("selectionRoot/" + sample + "/mumu/mumu" + filename); }
  LOG(logDEBUG) << "Looking for mass " << mass << " files, found " << filenames.size() << " to be opened.";

  TH1D* aMcHist;
  TFile * input = TFile::Open(filenames.front(),"read");
  if(!input->IsOpen()) {
    LOG(logCRITICAL) << "Failed to access data file " << filenames.front();
    throw 1;
  }

  LOG(logDEBUG) << "Getting NLO curve from " << filenames.front();
  aMcHist = dynamic_cast<TH1D*>(input->Get("VisGenTTBar1stJetMass")->Clone());
  delete input;

  for(std::vector<TString>::iterator file = filenames.begin()+1; file != filenames.end(); ++file) {
    LOG(logDEBUG) << "Getting NLO curve from " << *file;
    TFile * input2 = TFile::Open(*file,"read");
    if(!input2->IsOpen()) {
      LOG(logCRITICAL) << "Failed to access data file " << *file;
      throw 1;
    }
    aMcHist->Add(dynamic_cast<TH1D*>(input2->Get("VisGenTTBar1stJetMass")));
    delete input2;
  }

  TH1D* aMcBinned;

  // Histogram containing differential cross section from data (just for the binning):
  TH1D * aDiffXSecHist = dynamic_cast<TH1D*>(histos->Get("unfoldedHistNorm"));
  Int_t nbins = aDiffXSecHist->GetNbinsX();
  Double_t Xbins[nbins+1];
  for (Int_t bin = 1; bin <= nbins; bin++) Xbins[bin-1] = aDiffXSecHist->GetBinLowEdge(bin);
  Xbins[nbins] = aDiffXSecHist->GetBinLowEdge(nbins) + aDiffXSecHist->GetBinWidth(nbins);

  // Globally scaling the MC statistical errors by getting the overall weight from Intergal() and GetEntries():
  aMcBinned = dynamic_cast<TH1D*>(aMcHist->Rebin(nbins,"madgraphplot",Xbins));

  for (Int_t bin=0; bin < nbins; bin++) {
    // Divide rebinned histogram's bin content by bin width factor (new/old):
    aMcBinned->SetBinError(bin+1,sqrt(aMcBinned->GetBinContent(bin+1))/((Xbins[bin+1]-Xbins[bin])/aMcHist->GetBinWidth(1)));
    aMcBinned->SetBinContent(bin+1,aMcBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/aMcHist->GetBinWidth(1)));
  }
  aMcBinned->Scale(1./aMcBinned->Integral("width"));
  
  // Create a new histogram with the reduced binning:
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Simulation hist has " << nbins << " bins, using " << (nbins-startbin+1);

  TH1D * simulationHist = new TH1D("reco_" + channel + Form("_m%3.1f",mass),
				   "reco_" + channel + Form("_m%3.1f",mass),
				   nbins-startbin+1,
				   &Xbins[startbin-1]);

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    simulationHist->SetBinContent(bin+1-startbin,aMcBinned->GetBinContent(bin));
    simulationHist->SetBinError(bin+1-startbin,aMcBinned->GetBinError(bin));
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << simulationHist->GetBinContent(bin+1-startbin) << " (err=" << simulationHist->GetBinError(bin+1-startbin) << ")";
  }

  LOG(logDEBUG) << "Returning Simulation histogram now.";
  return simulationHist;
}

TFile * extractorDiffXSec::selectInputFile(TString sample, TString ch) {

  // Overwrite the samples unfolded with different masses with just the nominal:
  if((flags & FLAG_UNFOLD_ALLMASSES) == 0 ) { 
    LOG(logDEBUG2) << "Overwriting: " << sample;

    // Mass samples from Nominal:
    if(sample.Contains("GEV") || sample.Contains("Nominal")) {
      sample = getSampleFromMass("Nominal",unfoldingMass,true);
    }

    // All other mass-varied samples:
    else if(sample.Contains("POS") || sample.Contains("NEG")) { 
      sample.Remove(sample.Length()-5);
      // Check if we have set a specific sample to use for unfolding:
      sample = getSampleFromMass(sample,unfoldingMass,false);
    }

    LOG(logDEBUG2) << "with: " << sample;
  }

  // Input files for Differential Cross section mass extraction: unfolded distributions
  TString filename = "UnfoldingResults/" + sample + "/" + ch + "/HypTTBar1stJetMassResults.root";
  TFile * input = TFile::Open(filename,"read");
  if(!input->IsOpen()) {
    LOG(logCRITICAL) << "Failed to access data file " << filename;
    throw 1;
  }
  LOG(logDEBUG) << "Successfully opened file " << filename;
  return input;
}
