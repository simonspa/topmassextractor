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
#include "plotter.h"
#include "helpers.h"
#include "extractor.h"

using namespace unilog;
using namespace massextractor;

void extractorDiffXSec::setUnfoldingMass(Double_t mass) {
  LOG(logINFO) << "Requested to use m_t=" << mass << " for unfolding.";
  unfoldingMass = mass;
}

void extractorDiffXSecScaled::prepareScaleFactors(TString systematic) {

  Int_t sign = 0;
  if(systematic.Contains("UP")) { sign = 1; }
  else if(systematic.Contains("DOWN")) { sign = -1; }

  if(systematic.Contains("PDF")) {
    LOG(logDEBUG2) << "Preparing PDF uncertainty scale factors...";
    // PDF scaling for all five bins:
    if(m_channel == "combined") {
      scaleFactors.push_back(sign*0.0518874555267725);
      scaleFactors.push_back(sign*0.0307283998130669);
      scaleFactors.push_back(sign*0.0124526579196978);
      scaleFactors.push_back(sign*0.00149147332885075);
      scaleFactors.push_back(sign*0.00473063248751739);
    }
    else if(m_channel == "ee") {
      scaleFactors.push_back(sign*0.0520053995180302);
      scaleFactors.push_back(sign*0.0303351933038199);
      scaleFactors.push_back(sign*0.0126092037398061);
      scaleFactors.push_back(sign*0.00156712414239097);
      scaleFactors.push_back(sign*0.00486566872899394);
    }
    else if(m_channel == "emu") {
      scaleFactors.push_back(sign*0.0525837643853782);
      scaleFactors.push_back(sign*0.0303790574991724);
      scaleFactors.push_back(sign*0.0124525455719107);
      scaleFactors.push_back(sign*0.00151246828268544);
      scaleFactors.push_back(sign*0.00463153835628149);
    }
    else if(m_channel == "mumu") {
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
  std::vector<Double_t> Xbins = getBinningFromHistogram(aDiffXSecHist,startbin,nbins);

  TH1D * signalHist = new TH1D("diffxs_" + m_channel + Form("_m%3.1f",mass),
			       "diffxs_" + m_channel + Form("_m%3.1f",mass),
			       Xbins.size()-1,
			       &Xbins.front());

  // Store the binning for global use:
  bin_boundaries = Xbins;

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    Double_t signal = getSignal(bin-startbin,mass,aDiffXSecHist->GetBinContent(bin));
    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << aDiffXSecHist->GetBinContent(bin) << " signal=" << signal;

    signalHist->SetBinContent(bin+1-startbin,signal);
    // Scale the error, so that the relative statistical error stays the same:
    signalHist->SetBinError(bin+1-startbin,aDiffXSecHist->GetBinError(bin)*signal/aDiffXSecHist->GetBinContent(bin));
  }

  // Return DiffXSec signal histogram:
  return signalHist;
}

Double_t extractorDiffXSec::getSignal(Int_t /*bin*/, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return data;
}

Double_t extractorDiffXSecScaled::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return (1+scaleFactors.at(bin))*data;
}

TH1D * extractorDiffXSec::getSimulationHistogram(Double_t mass, TFile * histos) {

  std::vector<TString> channels;
  
  if(m_channel == "combined") { channels.push_back("ee"); channels.push_back("emu"); channels.push_back("mumu"); }
  else channels.push_back(m_channel);
  LOG(logDEBUG) << "Looking for mass " << mass << " files, found " << channels.size() << " files to be summed.";

  TH1D* aMcHist;
  TFile * input = selectInputFileTheory(channels.front(), getSampleFromMass(m_sample,mass,true));
  aMcHist = dynamic_cast<TH1D*>(input->Get("VisGenTTBar1stJetMass")->Clone());
  delete input;

  for(std::vector<TString>::iterator ch = channels.begin()+1; ch != channels.end(); ++ch) {
    TFile * input2 = selectInputFileTheory(*ch, getSampleFromMass(m_sample,mass,true));
    aMcHist->Add(dynamic_cast<TH1D*>(input2->Get("VisGenTTBar1stJetMass")));
    delete input2;
  }

  TH1D* aMcBinned;

  // Histogram containing differential cross section from data (just for the binning):
  TH1D * aDiffXSecHist = dynamic_cast<TH1D*>(histos->Get("unfoldedHistNorm"));
  Int_t nbins = aDiffXSecHist->GetNbinsX();
  std::vector<Double_t> Xbins = getBinningFromHistogram(aDiffXSecHist);

  // Globally scaling the MC statistical errors by getting the overall weight from Intergal() and GetEntries():
  aMcBinned = dynamic_cast<TH1D*>(aMcHist->Rebin(nbins,"madgraphplot",&Xbins.front()));

  for (Int_t bin=0; bin < nbins; bin++) {
    // Divide rebinned histogram's bin content by bin width factor (new/old):
    aMcBinned->SetBinError(bin+1,sqrt(aMcBinned->GetBinContent(bin+1))/((Xbins.at(bin+1)-Xbins.at(bin))/aMcHist->GetBinWidth(1)));
    aMcBinned->SetBinContent(bin+1,aMcBinned->GetBinContent(bin+1)/((Xbins.at(bin+1)-Xbins.at(bin))/aMcHist->GetBinWidth(1)));
  }
  aMcBinned->Scale(1./aMcBinned->Integral("width"));
  
  // Create a new histogram with the reduced binning:
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Simulation hist has " << nbins << " bins, using " << (nbins-startbin+1);

  TH1D * simulationHist = new TH1D("reco_" + m_channel + Form("_m%3.1f",mass),
				   "reco_" + m_channel + Form("_m%3.1f",mass),
				   Xbins.size()-startbin,
				   &Xbins.at(startbin-1));

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    simulationHist->SetBinContent(bin+1-startbin,aMcBinned->GetBinContent(bin));
    if((flags & FLAG_NO_THEORYPREDICITION_ERRORS) == 0) {
      Double_t err_match = aMcBinned->GetBinError(bin) + m_prediction_errors.at(bin-startbin).first;
      Double_t err_scale = aMcBinned->GetBinError(bin) + m_prediction_errors.at(bin-startbin).second;
      Double_t error = TMath::Sqrt(aMcBinned->GetBinError(bin)*aMcBinned->GetBinError(bin) + err_match*err_match + err_scale*err_scale);
      simulationHist->SetBinError(bin+1-startbin,error);
    }
    else {
      simulationHist->SetBinError(bin+1-startbin,aMcBinned->GetBinError(bin));
    }
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << simulationHist->GetBinContent(bin+1-startbin) << " (err=" << simulationHist->GetBinError(bin+1-startbin) << ")";
  }

  LOG(logDEBUG) << "Returning Simulation histogram now.";
  return simulationHist;
}

TFile * extractorDiffXSec::selectInputFileTheory(TString channel, TString sample) {

  TString filename;
  if(sample.Contains("MASS") || sample == "Nominal") {
    if(getMassFromSample(sample) < 167) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_166_massdown.root"; }
    else if(getMassFromSample(sample) < 170) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_169_massdown.root"; }
    else if(getMassFromSample(sample) < 172) { sample = "MASS_DOWN"; filename = "_ttbarsignalplustau_massdown.root"; }
    else if(getMassFromSample(sample) < 173) { sample = "Nominal"; filename = "_ttbarsignalplustau.root"; }
    else if(getMassFromSample(sample) < 174) { sample = "MASS_UP"; filename = "_ttbarsignalplustau_massup.root"; }
    else if(getMassFromSample(sample) < 176) { sample = "MASS_UP"; filename = "_ttbarsignalplustau_175_massup.root"; }
    else { sample = "MASS_UP"; filename = "_ttbarsignalplustau_178_massup.root"; }
  }
  else if(sample.Contains("MATCH_UP")) filename = "_ttbarsignalplustau_matchingup.root";
  else if(sample.Contains("MATCH_DOWN")) filename = "_ttbarsignalplustau_matchingdown.root";
  else if(sample.Contains("SCALE_UP")) filename = "_ttbarsignalplustau_scaleup.root";
  else if(sample.Contains("SCALE_DOWN")) filename = "_ttbarsignalplustau_scaledown.root";

  // Input files for Differential Cross section mass extraction: NLO curves
  TString path = "selectionRoot/" + sample + "/" + channel + "/" + channel + filename;
  LOG(logDEBUG) << "Getting NLO curve from " << path;
  TFile * input = OpenFile(path,"read");
  LOG(logDEBUG) << "Successfully opened file " << path;
  return input;
}

TFile * extractorDiffXSec::selectInputFile(TString sample) {

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
  TString filename = "UnfoldingResults/" + sample + "/" + m_channel + "/HypTTBar1stJetMassResults.root";
  TFile * input = OpenFile(filename,"read");
  LOG(logDEBUG) << "Successfully opened file " << filename;
  return input;
}

std::pair<TGraphErrors*,TF1*> extractorDiffXSec::getFittedChiSquare(std::vector<Double_t> masses, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc) {

  // Not using the covariance matrix method, but just summing the individual chi2 of the bins:
  if((flags & FLAG_DONT_USE_COVARIANCE) != 0) { return extractor::getFittedChiSquare(masses, data, mc); }

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TGraphErrors * chi2sum = new TGraphErrors();
  TString gname = "chi2_" + m_channel + "_sum";
  chi2sum->SetTitle(gname);

  // Cross-check: we have the same number of rho-S bins:
  if(data.size() != mc.size()) {
    LOG(logCRITICAL) << "Bin numbers don't match!";
    throw(1);
  }

  std::vector<TGraphErrors*> dataFits, mcFits;
  // Loop over all bins we have:
  for(size_t bin = 0; bin < data.size(); bin++) {
    TGraphErrors * dataFit = new TGraphErrors();
    TGraphErrors * mcFit = new TGraphErrors();

    // Get the likelihood for the two functions:
    createIntersectionChiSquare(data.at(bin),mc.at(bin),1+bin, dataFit, mcFit);
    dataFits.push_back(dataFit);
    mcFits.push_back(mcFit);
  }

  // Calculate the overall Chi2 using all bin correlations from covariance matrix:
  // Get the inverse covariance matrix:
  TMatrixD * invCov = getInverseCovMatrix(m_sample);
  // Now we have to fits to MC and data for all bins, let's calculate the Chi2.
  // Scanning over the fits:
  for(Int_t i = 0; i < dataFits.at(0)->GetN(); i++) {
    // Prepare the vector containing the (dat - mc) difference for all bins:
    TVectorD v(dataFits.size());
    Double_t x_dat;
    std::stringstream sv;
    for(size_t bin = 0; bin < dataFits.size(); bin++) {
      Double_t y_dat, y_mc;
      dataFits.at(bin)->GetPoint(i,x_dat,y_dat);
      mcFits.at(bin)->GetPoint(i,x_dat,y_mc);
      v[bin] = y_dat - y_mc;
      sv << v[bin] << " ";
    }
    // Do the matrix multiplication with the inverse covariance matrix:
    TVectorD v2 = v;
    v2 *= *invCov;
    Double_t chi2 = v2*v;
    LOG(logDEBUG4) << "(" << i << ") Chi2 = V^T*COV^-1*V = (" << x_dat << "/" << chi2 << ") "
		   << "(V:" << v.GetNoElements() << "x1, {" << sv.str() << "} "
		   << "COV: " << invCov->GetNrows() << "x" << invCov->GetNcols() << ")";
    chi2sum->SetPoint(i,x_dat,chi2);
  }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + m_channel + "_sum";
    c = new TCanvas(cname,cname);
    c->cd();
    chi2sum->SetMarkerStyle(20);
    chi2sum->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2sum->GetYaxis()->SetTitle("#chi^{2}");
    chi2sum->Draw("AP");
    DrawDecayChLabel(getChannelLabel(m_channel));
    DrawCMSLabels();
    chi2sum->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + cname + ".pdf"); }
  }

  // Fit the graph
  chi2sum->Fit("pol2","Q","",170,174.5);
  
  finalChiSquare = std::make_pair(chi2sum,chi2sum->GetFunction("pol2"));
  return finalChiSquare;
}

TMatrixD * extractorDiffXSec::getInverseCovMatrix(TString sample) {

  // Input files for Differential Cross section mass extraction: unfolded distributions
  TString filename = "SVD/" + sample + "/Unfolding_" + m_channel + "_TtBar_Mass_HypTTBar1stJetMass.root";
  TString histogramname = "SVD_" + m_channel + "_TtBar_Mass_HypTTBar1stJetMass_" + sample + "_STATCOV";

  TFile * input = OpenFile(filename,"read");
  LOG(logDEBUG) << "Successfully opened covariance matrix file " << filename;

  // Histogram containing differential cross section from data:
  TH2D * statCovNorm = static_cast<TH2D*>(input->Get(histogramname));
  if(!statCovNorm) { LOG(logCRITICAL) << "Failed to get histogram " << histogramname; }
  else { LOG(logDEBUG) << "Fetched " << histogramname; }

  // Cut away over and underflow bins:
  LOG(logDEBUG2) << "Creating a " << statCovNorm->GetNbinsX()-2 << "-by-" << statCovNorm->GetNbinsY()-2 << " COV matrix:";

  // Read the matrix from file:
  TMatrixD * cov = new TMatrixD(statCovNorm->GetNbinsX()-2,statCovNorm->GetNbinsY()-2);
  for(Int_t y = 0; y < cov->GetNcols(); y++) {
    std::stringstream st;
    for(Int_t x = 0; x < cov->GetNrows(); x++) {
      (*cov)(x,y) = statCovNorm->GetBinContent(x+2,y+2);
      st << std::setw(15) << std::setprecision(5) << (*cov)(x,y);
    }
    LOG(logDEBUG3) << st.str();
  }

  // Invert the covariance matrix:
  cov->Invert();

  LOG(logDEBUG2) << "Inverted COV matrix:";
  // Just printing is for debugging:
  for(Int_t y = 0; y < cov->GetNcols(); y++) {
    std::stringstream st;
    for(Int_t x = 0; x < cov->GetNrows(); x++) {
      st << std::setw(15) << std::setprecision(5) << (*cov)(x,y); 
    }
    LOG(logDEBUG3) << st.str();
  }

  return cov;
}

void extractorDiffXSec::getPredictionUncertainties() {

  m_prediction_errors.clear();
  LOG(logDEBUG2) << "Preparing theory prediction uncertainties for " << m_sample;

  std::vector<Double_t> aDiffXSecMatchUp = calcSampleDifference(m_sample,"MATCH_UP","VisGenTTBar1stJetMass",true);
  std::vector<Double_t> aDiffXSecScaleUp = calcSampleDifference(m_sample,"SCALE_UP","VisGenTTBar1stJetMass",true);

  std::vector<Double_t> aDiffXSecMatchDown = calcSampleDifference(m_sample,"MATCH_DOWN","VisGenTTBar1stJetMass",true);
  std::vector<Double_t> aDiffXSecScaleDown = calcSampleDifference(m_sample,"SCALE_DOWN","VisGenTTBar1stJetMass",true);

  // Take the theory prediction errors into account:
  for(size_t i = 0; i < aDiffXSecMatchUp.size(); ++i) {
    Double_t dxsec_match = (TMath::Abs(aDiffXSecMatchUp.at(i)) + TMath::Abs(aDiffXSecMatchDown.at(i)))/2;
    Double_t dxsec_scale = (TMath::Abs(aDiffXSecScaleUp.at(i)) + TMath::Abs(aDiffXSecScaleDown.at(i)))/2;

    m_prediction_errors.push_back(std::make_pair(dxsec_match,dxsec_scale));
    LOG(logDEBUG) << "Pred. err #" << i+1 << " diffxs = " << dxsec_match << " & " << dxsec_scale;
  }
}
