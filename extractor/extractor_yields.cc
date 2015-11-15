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
#include "helpers.h"
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

Double_t extractorYield::getReco(Int_t bin, Double_t /*mass*/, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Other (i.e. non-ttbar) backgrounds:
  Double_t bgr_other = bgr - ttbgr;

  // Scale factor (normalizing to data):
  if(mcScalingFactor < 0.0001) {
    LOG(logCRITICAL) << "Data for scaling missing!";
    throw(1);
  }

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    // Return a full pseudo data set including scaled backgrounds:
    Double_t reco_data = (reco + ttbgr)*mcScalingFactor + bgr_other;
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " reco_data=" << reco_data;
    return reco_data;
  }
  else {
    // Scale the reco according to the different TTBar Cross sections (mass dependent):
    Double_t corr_reco = reco*mcScalingFactor;
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

    // Return the reco event count corrected by the ttbar Xsec at given mass:
    return corr_reco;
  }
}

Double_t extractorYieldScaled::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  // Return the signal after scaling with supplied scale factor:
  return (1+scaleFactors.at(bin-1))*data;
}

Double_t extractorYieldPrediction::getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  Int_t sign = 0;
  if(m_systematic.Contains("UP")) { sign = 1; }
  else if(m_systematic.Contains("DOWN")) { sign = -1; }

  if(m_systematic.Contains("MATCH")) {
    LOG(logDEBUG3) << (sign > 0 ? "Add" : "Subtract") << " Match Prediction: " 
		   << reco << (sign > 0 ? "+" : "-") 
		   << m_prediction_errors_rec.at(bin-1).first
		   << "=" << (reco + sign*m_prediction_errors_rec.at(bin-1).first);
    return extractorYield::getReco(bin, mass,
				   reco + sign*m_prediction_errors_rec.at(bin-1).first,
				   bgr + sign*m_prediction_errors_bgr.at(bin-1).first,
				   ttbgr + sign*m_prediction_errors_ttbgr.at(bin-1).first);
  }
  else if(m_systematic.Contains("SCALE")) {
    LOG(logDEBUG3) << (sign > 0 ? "Add" : "Subtract") << " Scale Prediction: " 
		   << reco << (sign > 0 ? "+" : "-") 
		   << m_prediction_errors_rec.at(bin-1).second
		   << "=" << (reco + sign*m_prediction_errors_rec.at(bin-1).second);
    return extractorYield::getReco(bin, mass,
				   reco + sign*m_prediction_errors_rec.at(bin-1).second,
				   bgr + sign*m_prediction_errors_bgr.at(bin-1).second,
				   ttbgr + sign*m_prediction_errors_ttbgr.at(bin-1).second);
  }
  else return extractorYield::getReco(bin, mass, reco, bgr, ttbgr);
}

Double_t extractorYield::getPseudoData(Int_t /*bin*/, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  Double_t xsecCorrection = getTtbarXsec(mass)/getTtbarXsec(nominalmass);

  // Calculate the "others background" as difference between bgr and ttbgr (no xsec scaling)
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
  TFile * closureFile = selectInputFile(closure);

  // Histograms containing the background:
  TH1D * aRecHist = static_cast<TH1D*>(closureFile->Get("aRecHist"));
  TH1D * aTtBgrHist = static_cast<TH1D*>(closureFile->Get("aTtBgrHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(closureFile->Get("aBgrHist"));

  // Build a pseudo data set from the closure sample:
  pseudoData = new TH1D("pseudodata","pseudodata",aRecHist->GetNbinsX(),0,100);
  pseudoData->SetDirectory(0);
  Double_t mass = getMassFromSample(closure);

  LOG(logINFO) << "Running Closure test. Pseudo data taken from closure sample " << closure;
  LOG(logINFO) << "Pseudo data mass = " << mass;

  for(Int_t bin = 1; bin <= pseudoData->GetNbinsX(); bin++) {

    Double_t pdata = getPseudoData(bin,mass,aRecHist->GetBinContent(bin),aBgrHist->GetBinContent(bin),aTtBgrHist->GetBinContent(bin));

    LOG(logDEBUG) << "Closure: Bin #" << bin << " sig=" << aRecHist->GetBinContent(bin) << " pdat=" << pdata;
    // Write pseudo data with background:
    pseudoData->SetBinContent(bin,pdata);
    pseudoData->SetBinError(bin,aRecHist->GetBinError(bin));
  }

  // Set the original number of events of this closure sample (no weights):
  pseudoData->SetEntries(aRecHist->GetEntries() + aBgrHist->GetEntries());

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

  // Store the total number of events (no weights applied) to later scale the statistical uncertainty. 
  // m_stat_ndata = aDataHist->GetEntries();
  // For this, always take the number of events from the data sample
  // - otherwise we end up with ridiculously large errors just because
  // the MC samples are larger.
  m_stat_ndata = dynamic_cast<TH1D*>(histos->Get("aDataHist"))->GetEntries();

  // Create a new histogram with the same binning:
  Int_t nbins = aDataHist->GetNbinsX();
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Data hist has " << nbins << " bins, using " << (nbins-startbin+1);
  std::vector<Double_t> Xbins = getBinningFromHistogram(aRecHist,startbin, nbins);

  TH1D * signalHist = new TH1D(type + m_channel + Form("_m%3.1f",mass),
			       type + m_channel + Form("_m%3.1f",mass),
			       Xbins.size()-1,
			       &Xbins.front());

  // Store the binning for global use:
  bin_boundaries = Xbins;

  // Calculate the scaling factor from the histogram integrals: X = (reco + ttbgr)/(data - bgr_other);
  mcScalingFactor = (aRecHist->Integral() + aTtBgrHist->Integral())/(aDataHist->Integral() - aBgrHist->Integral() + aTtBgrHist->Integral());
  LOG(logDEBUG) << "mcScalingFactor = " << mcScalingFactor;

  // Iterate over all bins:
  for(Int_t bin = startbin; bin <= nbins; bin++) {

    // Get signal corrected by ttbar background:
    Double_t signal = getSignal(bin, mass, aDataHist->GetBinContent(bin),
				aRecHist->GetBinContent(bin),
				aBgrHist->GetBinContent(bin),
				aTtBgrHist->GetBinContent(bin));
    
    // Write background subtrated signal:
    signalHist->SetBinContent(bin-startbin+1,signal);
    // Scale the error, so that the relative statistical error stays the same:
    signalHist->SetBinError(bin-startbin+1,aDataHist->GetBinError(bin)*signal/aDataHist->GetBinContent(bin));
  }

  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) {
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
  std::vector<Double_t> Xbins = getBinningFromHistogram(aRecHist,startbin,nbins);

  TH1D * simulationHist = new TH1D("reco_" + m_channel + Form("_m%3.1f",mass),
				   "reco_" + m_channel + Form("_m%3.1f",mass),
				   Xbins.size()-1,
				   &Xbins.front());

  TH1D * matchErrHist = new TH1D("", "", Xbins.size()-1, &Xbins.front());
  TH1D * scaleErrHist = new TH1D("", "", Xbins.size()-1, &Xbins.front());
  
  // Iterate over all bins:
  for(Int_t bin = startbin; bin <= nbins; bin++) {
    // Correct the Reco events for different TTBar Cross sections (mass dependent):
    Double_t corr_reco = getReco(bin, mass,aRecHist->GetBinContent(bin),aBgrHist->GetBinContent(bin),aTtBgrHist->GetBinContent(bin));
    simulationHist->SetBinContent(bin-startbin+1,corr_reco);

    // Scale the statistical error, so that the relative statistical error stays the same, no prediction errors yet:
    Double_t err_stat = aRecHist->GetBinError(bin)*corr_reco/aRecHist->GetBinContent(bin);
    simulationHist->SetBinError(bin-startbin+1,err_stat);

    if((flags & FLAG_NO_THEORYPREDICTION_ERRORS) == 0) {
      // Also populate the histograms used for the precition error calculation:
      matchErrHist->SetBinContent(bin-startbin+1,getReco(bin,mass,aRecHist->GetBinContent(bin)+m_prediction_errors_rec.at(bin-startbin).first,
							 aBgrHist->GetBinContent(bin)+m_prediction_errors_bgr.at(bin-startbin).first,
							 aTtBgrHist->GetBinContent(bin)+m_prediction_errors_ttbgr.at(bin-startbin).first));
      scaleErrHist->SetBinContent(bin-startbin+1,getReco(bin,mass,aRecHist->GetBinContent(bin)+m_prediction_errors_rec.at(bin-startbin).second,
							 aBgrHist->GetBinContent(bin)+m_prediction_errors_bgr.at(bin-startbin).second,
							 aTtBgrHist->GetBinContent(bin)+m_prediction_errors_ttbgr.at(bin-startbin).second));
    }
  }

  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) {
    simulationHist->Scale(1./simulationHist->Integral("width"));
    if((flags & FLAG_NO_THEORYPREDICTION_ERRORS) == 0) {
      matchErrHist->Scale(1./matchErrHist->Integral("width"));
      scaleErrHist->Scale(1./scaleErrHist->Integral("width"));
    }
    LOG(logDEBUG) << "Normalized Reco hist.";
  }

  // Recalculate the error including prediction uncertainties if requested:
  if((flags & FLAG_NO_THEORYPREDICTION_ERRORS) == 0) {
    for(Int_t bin = 1; bin <= simulationHist->GetNbinsX(); bin++) {
      Double_t err_match = simulationHist->GetBinContent(bin) - matchErrHist->GetBinContent(bin);
      Double_t err_scale = simulationHist->GetBinContent(bin) - scaleErrHist->GetBinContent(bin);
      Double_t err_stat = simulationHist->GetBinError(bin);
      Double_t error = TMath::Sqrt(err_stat*err_stat + err_match*err_match + err_scale*err_scale);
      simulationHist->SetBinError(bin-startbin+1,error);
    }
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

void extractorYieldScaled::prepareScaleFactors(TString systematic) {

  Int_t sign = 0;
  if(systematic.Contains("UP")) { sign = 1; }
  else if(systematic.Contains("DOWN")) { sign = -1; }

  if(systematic.Contains("PDF")) {
    LOG(logDEBUG2) << "Preparing PDF uncertainty scale factors...";
    scaleFactors = getPDFScaleFactors(sign,m_channel);
  }
  std::stringstream sv;
  for(std::vector<Double_t>::iterator sf = scaleFactors.begin(); sf != scaleFactors.end(); ++sf) { sv << *sf << " "; }
  LOG(logDEBUG2) << "Scale factors for every bin: " << sv.str();
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

Double_t extractorYieldBackground::getReco(Int_t bin, Double_t /*mass*/, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Other (i.e. non-ttbar) backgrounds:
  Double_t bgr_other = (bgr - ttbgr)*scaleFactor;

  // Scale factor (normalizing to data):
  if(mcScalingFactor < 0.0001) {
    LOG(logCRITICAL) << "Data for scaling missing!";
    throw(1);
  }

  if((flags & FLAG_DONT_SUBTRACT_BACKGROUND) != 0) {
    // Return a full pseudo data set including scaled backgrounds:
    Double_t reco_data = (reco + ttbgr)*mcScalingFactor + bgr_other;
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " reco_data=" << reco_data;
    return reco_data;
  }
  else {
    // Scale the reco according to the different TTBar Cross sections (mass dependent):
    Double_t corr_reco = reco*mcScalingFactor;
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

    // Return the reco event count corrected by the ttbar Xsec at given mass:
    return corr_reco;
  }
}

Double_t extractorYieldOtherSamples::getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Subtract the difference in event count for the nominal mass bin and systematics variation for every bin:
  LOG(logDEBUG3) << "Subtracting difference to Nominal: "
		 << reco << "-" << deltaRec.at(bin-1) << "=" << (reco-deltaRec.at(bin-1));
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

  // Get the total number of events (no weights applied) to later scale the statistical uncertainty:
  TFile * sysFile = selectInputFile(systematic);
  m_stat_nmc = static_cast<TH1D*>(sysFile->Get("aRecHist"))->GetEntries() + static_cast<TH1D*>(sysFile->Get("aBgrHist"))->GetEntries();
  delete sysFile;

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
    LOG(logDEBUG3) << "Diff bin #" << i+1 << " reco=" << deltaRec.at(i) << " bgr=" << deltaBgr.at(i) << " ttbgr=" << deltaTtbgr.at(i);
  }
}

std::vector<Double_t> extractorYield::calcSampleDifference(TString nominal, TString systematic, TString histogram) {

  std::vector<Double_t> difference;

  // Input files:
  TFile * nominalfile, * systematicfile;
  nominalfile = selectInputFile(nominal);
  systematicfile = selectInputFile(systematic);

  LOG(logDEBUG) << "Difference: " << nominal << " to " << systematic << ", hist " << histogram;

  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  TH1D * nominalHistogram = static_cast<TH1D*>(nominalfile->Get(histogram));
  TH1D * systHistogram = static_cast<TH1D*>(systematicfile->Get(histogram));

  // Calculate the difference:
  for(Int_t bin = 1; bin <= nominalHistogram->GetNbinsX(); bin++) {

    Double_t diff = static_cast<Double_t>(nominalHistogram->GetBinContent(bin)) - static_cast<Double_t>(systHistogram->GetBinContent(bin));
    difference.push_back(diff);
    LOG(logDEBUG3) << "Hist " << histogram << ": diff bin #" << bin << " " << nominalHistogram->GetBinContent(bin) << " - " << systHistogram->GetBinContent(bin) << " = " << diff;
  }

  delete nominalfile;
  delete systematicfile;
  return difference;
}

void extractorYield::getPredictionUncertainties() {

  m_prediction_errors_rec.clear();
  m_prediction_errors_bgr.clear();
  m_prediction_errors_ttbgr.clear();

  LOG(logDEBUG) << "Preparing theory prediction uncertainties.";

  std::vector<Double_t> aRecMatchUp = calcSampleDifference("Nominal","MATCH_UP","aRecHist");
  std::vector<Double_t> aBgrMatchUp = calcSampleDifference("Nominal","MATCH_UP","aBgrHist");
  std::vector<Double_t> aTtBgrMatchUp = calcSampleDifference("Nominal","MATCH_UP","aTtBgrHist");

  std::vector<Double_t> aRecScaleUp = calcSampleDifference("Nominal","SCALE_UP","aRecHist");
  std::vector<Double_t> aBgrScaleUp = calcSampleDifference("Nominal","SCALE_UP","aBgrHist");
  std::vector<Double_t> aTtBgrScaleUp = calcSampleDifference("Nominal","SCALE_UP","aTtBgrHist");

  std::vector<Double_t> aRecMatchDown = calcSampleDifference("Nominal","MATCH_DOWN","aRecHist");
  std::vector<Double_t> aBgrMatchDown = calcSampleDifference("Nominal","MATCH_DOWN","aBgrHist");
  std::vector<Double_t> aTtBgrMatchDown = calcSampleDifference("Nominal","MATCH_DOWN","aTtBgrHist");

  std::vector<Double_t> aRecScaleDown = calcSampleDifference("Nominal","SCALE_DOWN","aRecHist");
  std::vector<Double_t> aBgrScaleDown = calcSampleDifference("Nominal","SCALE_DOWN","aBgrHist");
  std::vector<Double_t> aTtBgrScaleDown = calcSampleDifference("Nominal","SCALE_DOWN","aTtBgrHist");

  // Take the theory prediction errors into account:
  for(size_t i = 0; i < aRecMatchUp.size(); ++i) {
    //Symmetrized errors for up/down:
    Double_t rec_match = (TMath::Abs(aRecMatchUp.at(i)) + TMath::Abs(aRecMatchDown.at(i)))/2;
    Double_t bgr_match = (TMath::Abs(aBgrMatchUp.at(i)) + TMath::Abs(aBgrMatchDown.at(i)))/2;
    Double_t ttbgr_match = (TMath::Abs(aTtBgrMatchUp.at(i)) + TMath::Abs(aTtBgrMatchDown.at(i)))/2;

    Double_t rec_scale = (TMath::Abs(aRecScaleUp.at(i)) + TMath::Abs(aRecScaleDown.at(i)))/2;
    Double_t bgr_scale = (TMath::Abs(aBgrScaleUp.at(i)) + TMath::Abs(aBgrScaleDown.at(i)))/2;
    Double_t ttbgr_scale = (TMath::Abs(aTtBgrScaleUp.at(i)) + TMath::Abs(aTtBgrScaleDown.at(i)))/2;

    m_prediction_errors_rec.push_back(std::make_pair(rec_match,rec_scale));
    m_prediction_errors_bgr.push_back(std::make_pair(bgr_match,bgr_scale));
    m_prediction_errors_ttbgr.push_back(std::make_pair(ttbgr_match,ttbgr_scale));
    LOG(logDEBUG) << "Pred.err. #" << i+1 << " reco = " << rec_match << " & " << rec_scale 
		   << ", bgr = " << bgr_match << " & " << bgr_scale
		  << ", ttbgr = " << ttbgr_match << " & " << ttbgr_scale;
  }
}

Double_t extractorYieldOtherSamples::getStatError(Double_t &statPos, Double_t &statNeg) {

  Double_t statisticsScaleFactor = m_stat_ndata/m_stat_nmc;
  LOG(logINFO) << "Events in samples: " << (doClosure ? "Pseudo Data" : "Data" ) << ": " << m_stat_ndata << ", "
	       << m_systematic << ": " << m_stat_nmc
	       << " - Statistics Scale Factor: " << statisticsScaleFactor;

  statPos = statErrorPos*statisticsScaleFactor;
  statNeg = statErrorNeg*statisticsScaleFactor;

  return (statPos+statNeg)/2;
}
