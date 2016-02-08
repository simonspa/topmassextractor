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
    scaleFactors = getPDFScaleFactors(sign,m_channel);
  }
  std::stringstream sv;
  for(std::vector<Double_t>::iterator sf = scaleFactors.begin(); sf != scaleFactors.end(); ++sf) { sv << *sf << " "; }
  LOG(logDEBUG2) << "Scale factors for every bin: " << sv.str();
}

TH1D * extractorDiffXSec::getSignalHistogram(Double_t mass, TFile * histos) {

  // Histogram containing differential cross section from data:
  TH1D * aDiffXSecHist;
  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { aDiffXSecHist = static_cast<TH1D*>(histos->Get("unfoldedHistNorm")); }
  else { aDiffXSecHist = static_cast<TH1D*>(histos->Get("unfoldedHist")); }
  

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
    signalHist->SetBinContent(bin+1-startbin,signal);
    // Scale the error, so that the relative statistical error stays the same:
    Double_t staterr = aDiffXSecHist->GetBinError(bin)*signal/aDiffXSecHist->GetBinContent(bin);
    signalHist->SetBinError(bin+1-startbin,staterr);

    LOG(logDEBUG3) << "Bin #" << bin << ": data=" << aDiffXSecHist->GetBinContent(bin) << "+-" << aDiffXSecHist->GetBinError(bin) 
		   << " signal=" << signal << "+-" << staterr;
  }

  // Store the total number of events in the signal histogram:
  m_nevents = signalHist->Integral();
  LOG(logDEBUG2) << "Signal histogram integral: " << m_nevents;

  // Return DiffXSec signal histogram:
  return signalHist;
}

Double_t extractorDiffXSec::getSignal(Int_t /*bin*/, Double_t /*mass*/, Double_t data, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return data;
}

Double_t extractorDiffXSec::getReco(Int_t /*bin*/, Double_t /*mass*/, Double_t reco, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return reco;
}

Double_t extractorDiffXSecScaled::getReco(Int_t bin, Double_t /*mass*/, Double_t reco, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  return (1+scaleFactors.at(bin))*reco;
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
    // Since Sumw2() has been called for the histograms, Add() does recalc the errors correctly:
    aMcHist->Add(dynamic_cast<TH1D*>(input2->Get("VisGenTTBar1stJetMass")));
    delete input2;
  }

  TH1D* aMcBinned;

  // Histogram containing differential cross section from data (just for the binning):
  TH1D * aDiffXSecHist;
  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { aDiffXSecHist = dynamic_cast<TH1D*>(histos->Get("unfoldedHistNorm")); }
  else { aDiffXSecHist = dynamic_cast<TH1D*>(histos->Get("unfoldedHist")); }

  Int_t nbins = aDiffXSecHist->GetNbinsX();
  std::vector<Double_t> Xbins = getBinningFromHistogram(aDiffXSecHist);

  // Globally scaling the MC statistical errors by getting the overall weight from Integral() and GetEntries():
  aMcBinned = dynamic_cast<TH1D*>(aMcHist->Rebin(nbins,"madgraphplot",&Xbins.front()));

  for (Int_t bin=0; bin < nbins; bin++) {
    // Divide rebinned histogram's bin content and bin error by bin width factor (new/old):
    aMcBinned->SetBinError(bin+1,aMcBinned->GetBinError(bin+1)/((Xbins.at(bin+1)-Xbins.at(bin))/aMcHist->GetBinWidth(1)));
    aMcBinned->SetBinContent(bin+1,aMcBinned->GetBinContent(bin+1)/((Xbins.at(bin+1)-Xbins.at(bin))/aMcHist->GetBinWidth(1)));
  }
  
  // Either normalize the histogram to 1:
  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { aMcBinned->Scale(1./aMcBinned->Integral("width")); }
  // or scale to data:
  else { aMcBinned->Scale(aDiffXSecHist->Integral("width")/aMcBinned->Integral("width")); }
  
  // Create a new histogram with the reduced binning:
  Int_t startbin = 1;
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0) { startbin = nbins; }
  LOG(logDEBUG) << "Simulation hist has " << nbins << " bins, using " << (nbins-startbin+1);

  TH1D * simulationHist = new TH1D("reco_" + m_channel + Form("_m%3.1f",mass),
				   "reco_" + m_channel + Form("_m%3.1f",mass),
				   Xbins.size()-startbin,
				   &Xbins.at(startbin-1));

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    simulationHist->SetBinContent(bin+1-startbin,getReco(bin-startbin,mass,aMcBinned->GetBinContent(bin)));
    if((flags & FLAG_NO_THEORYPREDICTION_ERRORS) == 0) {
      Double_t err_match = aMcBinned->GetBinError(bin) + m_prediction_errors.at(bin-startbin).first;
      Double_t err_scale = aMcBinned->GetBinError(bin) + m_prediction_errors.at(bin-startbin).second;
      Double_t error = TMath::Sqrt(aMcBinned->GetBinError(bin)*aMcBinned->GetBinError(bin) + err_match*err_match + err_scale*err_scale);
      simulationHist->SetBinError(bin+1-startbin,error);
    }
    else {
      simulationHist->SetBinError(bin+1-startbin,aMcBinned->GetBinError(bin));
    }
    LOG(logDEBUG3) << "Bin #" << bin << ": reco=" << simulationHist->GetBinContent(bin+1-startbin) << "+-" << simulationHist->GetBinError(bin+1-startbin);
  }

  LOG(logDEBUG) << "Returning Simulation histogram now.";
  return simulationHist;
}

TFile * extractorDiffXSec::selectInputFileTheory(TString channel, TString sample) {

  TString filename;

  if((flags & FLAG_USE_NLO) != 0) {
    // Input files for Differential Cross section mass extraction: NLO curves
    if(sample.Contains("MASS") || sample == "Nominal") {
      if(getMassFromSample(sample) < 164) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m163p5.root"; }
      else if(getMassFromSample(sample) < 167) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m166p5.root"; }
      else if(getMassFromSample(sample) < 170) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m169p5.root"; }
      else if(getMassFromSample(sample) < 172) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m171p5.root"; }
      else if(getMassFromSample(sample) < 173) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m172p5.root"; }
      else if(getMassFromSample(sample) < 174) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m173p5.root"; }
      else if(getMassFromSample(sample) < 176) { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m175p5.root"; }
      else { sample = "POWHEG"; filename = "_ttbarsignalplustau_powhegbox_m178p5.root"; }
    }
    else if(sample.Contains("MATCH_UP")) filename = "_ttbarsignalplustau_matchingup.root";
    else if(sample.Contains("MATCH_DOWN")) filename = "_ttbarsignalplustau_matchingdown.root";
    else if(sample.Contains("SCALE_UP")) filename = "_ttbarsignalplustau_scaleup.root";
    else if(sample.Contains("SCALE_DOWN")) filename = "_ttbarsignalplustau_scaledown.root";
    else if(sample.Contains("POWHEG")) { sample = "POWHEG_lo"; filename = "_ttbarsignalplustau_powheg.root"; }
    else if(sample.Contains("MCATNLO")) filename = "_ttbarsignalplustau_mcatnlo.root";
    else if(sample.Contains("PERUGIA11mpiHi")) filename = "_ttbarsignalplustau_Perugia11mpiHi.root";
    else if(sample.Contains("PERUGIA11TeV")) filename = "_ttbarsignalplustau_Perugia11TeV.root";
    else if(sample.Contains("PERUGIA11NoCR")) filename = "_ttbarsignalplustau_Perugia11NoCR.root";
    else if(sample.Contains("PERUGIA11")) filename = "_ttbarsignalplustau_Perugia11.root";
  } else {
    // Input files for Differential Cross section mass extraction: MadGraph LO curves
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
  }

  TString path = "selectionRoot/" + sample + "/" + channel + "/" + channel + filename;
  TFile * input = OpenFile(path,"read");
  LOG(logDEBUG) << "Successfully opened file " << path;
  return input;
}

TFile * extractorDiffXSecGenLevelPrediction::selectInputFileTheory(TString channel, TString sample) {

  TString filename;
  TString samplename = m_systematic;

  if((flags & FLAG_USE_NLO) != 0) {
    samplename.Prepend("POWHEG_");
    // Input files for Differential Cross section mass extraction: NLO curves
    if(getMassFromSample(sample) < 164) { filename = "_ttbarsignalplustau_powhegbox_m163p5.root"; }
    else if(getMassFromSample(sample) < 167) { filename = "_ttbarsignalplustau_powhegbox_m166p5.root"; }
    else if(getMassFromSample(sample) < 170) { filename = "_ttbarsignalplustau_powhegbox_m169p5.root"; }
    else if(getMassFromSample(sample) < 172) { filename = "_ttbarsignalplustau_powhegbox_m171p5.root"; }
    else if(getMassFromSample(sample) < 173) { filename = "_ttbarsignalplustau_powhegbox_m172p5.root"; }
    else if(getMassFromSample(sample) < 174) { filename = "_ttbarsignalplustau_powhegbox_m173p5.root"; }
    else if(getMassFromSample(sample) < 176) { filename = "_ttbarsignalplustau_powhegbox_m175p5.root"; }
    else { filename = "_ttbarsignalplustau_powhegbox_m178p5.root"; }
  } else {
    LOG(logCRITICAL) << "Currently no per-mass samples are available for systematic variations at gen level!";
    throw;
  }

  TString path = "selectionRoot/" + samplename + "/" + channel + "/" + channel + filename;
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
  TGraphErrors * chi2sum_plotting = new TGraphErrors();
  TString gname = "chi2_" + m_channel + "_sum";
  chi2sum->SetTitle(gname);
  chi2sum_plotting->SetTitle(gname);

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

  // Bin number (starting from 1) to be dropped in order to
  // satistfy NDOF reduction due to the usage of normalized distributions
  // 0 returns all bins.
  Int_t drop_bin = 0;
  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { drop_bin = 1; }

  // Calculate the overall Chi2 using all bin correlations from covariance matrix:
  // Get the covariance matrix:
  TMatrixD * cov = getCovMatrix(m_sample);

  // Now we have to fits to MC and data for all bins, let's calculate the Chi2.
  // Scanning over the fits:
  Int_t i_plotting = 0;
  for(Int_t i = 0; i < dataFits.at(0)->GetN(); i++) {
    // Prepare the vector containing the (dat - mc) difference for all bins:
    TVectorD v(dataFits.size() - (drop_bin > 0 ? 1 : 0));
    std::vector<Double_t> mcstats;

    Double_t x_dat;
    std::stringstream sv;
    size_t output_bin = 0;
    for(size_t bin = 0; bin < dataFits.size(); bin++) {

      // Store the fit error for the current scan point:
      mcstats.push_back(mcFits.at(bin)->GetErrorY(i));

      // Drop bin for NDOF satisfaction:
      if(static_cast<Int_t>(bin) == drop_bin-1) { continue; }

      Double_t y_dat, y_mc;
      dataFits.at(bin)->GetPoint(i,x_dat,y_dat);
      mcFits.at(bin)->GetPoint(i,x_dat,y_mc);
      v[output_bin] = y_dat - y_mc;
      sv << v[output_bin] << " ";
      output_bin++;
    }

    // Invert the Covariance matrix with the MC statistical errors from the current scan point added:
    TMatrixD * invCov = invertCovMatrix(cov, mcstats, drop_bin);

    // Do the matrix multiplication with the inverse covariance matrix:
    TVectorD v2 = v;
    v2 *= *invCov;
    Double_t chi2 = v2*v;
    LOG(logDEBUG4) << "(" << i << ") Chi2 = V^T*COV^-1*V = (" << x_dat << "/" << chi2 << ") "
		   << "(V:" << v.GetNoElements() << "x1, {" << sv.str() << "} "
		   << "COV: " << invCov->GetNrows() << "x" << invCov->GetNcols() << ")";
    chi2sum->SetPoint(i,x_dat,chi2);

    // Check if this is within the range we want to plot:
    if(PLOT_LOWER_LIMIT_DIFFXS <= x_dat && x_dat <= PLOT_UPPER_LIMIT_DIFFXS && i%(dataFits.at(0)->GetN()/500) == 0) {
      chi2sum_plotting->SetPoint(i_plotting,x_dat,chi2);
      i_plotting++;
    }
    delete invCov;
  }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + m_channel + "_sum";
    c = new TCanvas(cname,cname);
    c->cd();
    chi2sum_plotting->SetMarkerStyle(20);
    chi2sum_plotting->SetLineWidth(3);
    chi2sum_plotting->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2sum_plotting->GetYaxis()->SetTitle("#chi^{2}");
    chi2sum_plotting->Draw("AL");
    DrawDecayChLabel(m_channel,((flags & FLAG_DONT_PLOT_CHANNELLABELS) == 0),-1);
    DrawCMSLabels();
    rescaleGraph(chi2sum_plotting);
    chi2sum_plotting->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + cname + ".pdf"); }
  }

  // Fit the graph
  chi2sum->Fit("pol2","Q","",170,174.5);
  
  finalChiSquare = std::make_pair(chi2sum,chi2sum->GetFunction("pol2"));
  return finalChiSquare;
}

TMatrixD * extractorDiffXSec::readMatrix(TString sample, TString channel) {

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
  TString filename = "SVD/" + sample + "/Unfolding_" + channel + "_TtBar_Mass_HypTTBar1stJetMass.root";
  TString histogramname;
  if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) == 0) { histogramname = "SVD_" + channel + "_TtBar_Mass_HypTTBar1stJetMass_" + sample + "_STATCOV"; }
  else { histogramname = "SVD_" + channel + "_TtBar_Mass_HypTTBar1stJetMass_" + sample + "_STATCOVNORM"; }

  TFile * input = OpenFile(filename,"read");
  LOG(logDEBUG) << "Successfully opened covariance matrix file " << filename;

  // Histogram containing differential cross section from data:
  TH2D * statCovNorm = static_cast<TH2D*>(input->Get(histogramname));
  if(!statCovNorm) { LOG(logCRITICAL) << "Failed to get histogram " << histogramname; }
  else { LOG(logDEBUG) << "Fetched " << histogramname; }

  // Cut away over and underflow bins:
  LOG(logDEBUG2) << "Creating a " << statCovNorm->GetNbinsX()-2 << "-by-" << statCovNorm->GetNbinsY()-2 << " COV matrix";

  // Build the matrix:
  TMatrixD * cov = new TMatrixD(statCovNorm->GetNbinsX()-2,statCovNorm->GetNbinsY()-2);
  for(Int_t y = 0; y < cov->GetNcols(); y++) {
    std::stringstream st;
    for(Int_t x = 0; x < cov->GetNrows(); x++) {
      (*cov)(x,y) = statCovNorm->GetBinContent(x+2,y+2)/(statCovNorm->GetXaxis()->GetBinWidth(x+2)*statCovNorm->GetYaxis()->GetBinWidth(y+2));
      st << std::setw(15) << std::setprecision(5) << (*cov)(x,y);
    }
    LOG(logDEBUG3) << st.str();
  }

  LOG(logDEBUG2) << "COV Determinant: " << std::setprecision(5) << cov->Determinant();
  LOG(logDEBUG2) << "Sqrt of COV Diagonals (stat. errors per bin): ";
  std::stringstream sigma;
  for(Int_t y = 0; y < cov->GetNcols(); y++) {
    for(Int_t x = 0; x < cov->GetNrows(); x++) {
      if(y == x) sigma << std::setw(15) << std::setprecision(5) << TMath::Sqrt((*cov)(x,y));
    }
  }
  LOG(logDEBUG2) << sigma.str();

  // Close the COV input file and set pointer back to output file:
  if(input->IsOpen()) {
    input->Close();
    if(m_root_output != NULL && m_root_output->IsOpen()) m_root_output->cd();
  }

  return cov;
}

TMatrixD * extractorDiffXSec::getCovMatrix(TString sample) {

  TMatrixD * cov;
  // Single channel - just read the information:
  if(m_channel != "combined") cov = readMatrix(sample,m_channel);
  // Combined channel - add up the three contributions:
  else {
    TMatrixD * cov_ee = readMatrix(sample,"ee");
    TMatrixD * cov_emu = readMatrix(sample,"emu");
    TMatrixD * cov_mumu = readMatrix(sample,"mumu");

    // New matrix containing the sum:
    cov = new TMatrixD(cov_ee->GetNrows(),cov_ee->GetNcols());
    for(Int_t y = 0; y < cov_ee->GetNcols(); y++) {
      std::stringstream st;
      for(Int_t x = 0; x < cov_ee->GetNrows(); x++) {
	// If the COV is not normalized, we can just linearly add the three matrices:
	if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) == 0) {
	  (*cov)(x,y) = (*cov_ee)(x,y) + (*cov_emu)(x,y) + (*cov_mumu)(x,y);
	}
	// If the COVs are normalized we have to inverse the contributions for the sum:
	else {
	  (*cov)(x,y) = 1/(1/(*cov_ee)(x,y) + 1/(*cov_emu)(x,y) + 1/(*cov_mumu)(x,y));
	}
	st << std::setw(15) << std::setprecision(5) << (*cov)(x,y);
      }
      LOG(logDEBUG3) << st.str();
    }

    LOG(logDEBUG2) << "COV Determinant: " << std::setprecision(5) << cov->Determinant();
    LOG(logDEBUG2) << "Sqrt of COV Diagonals (stat. errors per bin): ";
    std::stringstream sigma;
    for(Int_t y = 0; y < cov->GetNcols(); y++) {
      for(Int_t x = 0; x < cov->GetNrows(); x++) {
	if(y == x) sigma << std::setw(15) << std::setprecision(5) << TMath::Sqrt((*cov)(x,y));
      }
    }
    LOG(logDEBUG2) << sigma.str();
  }

  return cov;
}

TMatrixD * extractorDiffXSec::invertCovMatrix(TMatrixD * inputcov, std::vector<Double_t> mcstats, Int_t drop_bin) {

  TMatrixD * cov = new TMatrixD(*inputcov);

  if(static_cast<size_t>(cov->GetNcols()) != mcstats.size()) {
    LOG(logCRITICAL) << "Bin numbers don't match!";
    throw(1);
  }

  // Add MC statistical uncertainty to diagonal elements:
  LOG(logDEBUG3) << "Added MC contributions to diagonal elements";
  std::stringstream mc;
  for(Int_t y = 0; y < cov->GetNcols(); y++) {
    std::stringstream st;
    for(Int_t x = 0; x < cov->GetNrows(); x++) {
      if(x == y) {
	// Add squared mcstat error:
	(*cov)(x,y) = (*cov)(x,y) + mcstats.at(x)*mcstats.at(x);
	st << std::setw(15) << std::setprecision(5) << (*cov)(x,y);
	mc << std::setw(15) << std::setprecision(5) << mcstats.at(x);
      }
      else st << std::setw(15) << std::setprecision(5) << " ";
    }
    LOG(logDEBUG4) << st.str();
  }
  LOG(logDEBUG3) << "MC stat errors: " << mc.str();
 
  // Invert the covariance matrix:
  cov->Invert();
  LOG(logDEBUG3) << "Inverted COV matrix";

  TMatrixD * covFinal;
  // Check if we need to drop one bin. This might be necessary due to normalization which
  // reduces the NDOF by one.
  if(drop_bin > 0 && drop_bin <= cov->GetNrows()) {
    covFinal = new TMatrixD(cov->GetNrows()-1, cov->GetNcols()-1);
    LOG(logDEBUG4) << "Dropping bin " << drop_bin << " from convariance matrix.";

    // Now remove the additional bin:
    Int_t cov_x = 0, cov_y = 0;
    for(Int_t y = 0; y < covFinal->GetNcols(); y++) {
      // Check if we need to drop this entry:
      if(y == drop_bin-1) { cov_y++; }

      for(Int_t x = 0; x < covFinal->GetNrows(); x++) {
	// Check if we need to drop this entry:
	if(x == drop_bin-1) { cov_x++; }

	(*covFinal)(x,y) = (*cov)(cov_x,cov_y);
	cov_x++;
      }
      cov_x = 0; cov_y++;
    }
    LOG(logDEBUG3) << "Reduced inverted COV matrix to " << covFinal->GetNrows() << "-by-" << covFinal->GetNcols();
  }
  else if(drop_bin > 0) {
    LOG(logCRITICAL) << "Failed to drop COV bin" << drop_bin << ". Continue with full COV matrix.";
    covFinal = cov;
  }
  else {
    covFinal = cov;
  }

  // Just printing is for debugging:
  for(Int_t y = 0; y < covFinal->GetNcols(); y++) {
    std::stringstream st;
    for(Int_t x = 0; x < covFinal->GetNrows(); x++) {
      st << std::setw(15) << std::setprecision(5) << (*covFinal)(x,y); 
    }
    LOG(logDEBUG4) << st.str();
  }

  // Return the matrix:
  return covFinal;
}


std::vector<Double_t> extractorDiffXSec::calcSampleDifference(TString nominal, TString systematic, TString histogram) {

  std::vector<Double_t> difference;

  std::vector<TString> channels;
  if(m_channel == "combined") { channels.push_back("ee"); channels.push_back("emu"); channels.push_back("mumu"); }
  else channels.push_back(m_channel);
  LOG(logDEBUG) << "Found " << channels.size() << " files to be summed.";

  // Input files:
  TFile * nominalfile, * systematicfile, *binningfile;
  binningfile = selectInputFile(nominal);
  LOG(logDEBUG) << "Difference: " << nominal << " to " << systematic << ", hist " << histogram;

  nominalfile = selectInputFileTheory(channels.front(),nominal);
  systematicfile = selectInputFileTheory(channels.front(),systematic);
  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  TH1D * nominalHistogram = static_cast<TH1D*>(nominalfile->Get(histogram)->Clone());
  TH1D * systHistogram = static_cast<TH1D*>(systematicfile->Get(histogram)->Clone());

  for(std::vector<TString>::iterator ch = channels.begin()+1; ch != channels.end(); ++ch) {
    nominalfile = selectInputFileTheory(*ch, nominal);
    systematicfile = selectInputFileTheory(*ch, systematic);
    nominalHistogram->Add(dynamic_cast<TH1D*>(nominalfile->Get(histogram)));
    systHistogram->Add(dynamic_cast<TH1D*>(systematicfile->Get(histogram)));
  }

  TH1D * nominalHistBinned, * systHistBinned;

  // Rebin both nominal and systematic histograms:
  TH1D * aDiffXSecHist = dynamic_cast<TH1D*>(binningfile->Get("unfoldedHistNorm"));
  Int_t nbins = aDiffXSecHist->GetNbinsX();
  std::vector<Double_t> Xbins = getBinningFromHistogram(aDiffXSecHist);

  // Globally scaling the MC statistical errors by getting the overall weight from Intergal() and GetEntries():
  nominalHistBinned = dynamic_cast<TH1D*>(nominalHistogram->Rebin(nbins,"madgraphplot",&Xbins.front()));
  systHistBinned = dynamic_cast<TH1D*>(systHistogram->Rebin(nbins,"madgraphplot",&Xbins.front()));

  for (Int_t bin=0; bin < nbins; bin++) {
    // Divide rebinned histogram's bin content by bin width factor (new/old):
    nominalHistBinned->SetBinError(bin+1,sqrt(nominalHistBinned->GetBinContent(bin+1))/((Xbins.at(bin+1)-Xbins.at(bin))/nominalHistogram->GetBinWidth(1)));
    systHistBinned->SetBinError(bin+1,sqrt(systHistBinned->GetBinContent(bin+1))/((Xbins.at(bin+1)-Xbins.at(bin))/systHistogram->GetBinWidth(1)));

    nominalHistBinned->SetBinContent(bin+1,nominalHistBinned->GetBinContent(bin+1)/((Xbins.at(bin+1)-Xbins.at(bin))/nominalHistogram->GetBinWidth(1)));
    systHistBinned->SetBinContent(bin+1,systHistBinned->GetBinContent(bin+1)/((Xbins.at(bin+1)-Xbins.at(bin))/systHistogram->GetBinWidth(1)));
  }
  nominalHistBinned->Scale(1./nominalHistBinned->Integral("width"));
  systHistBinned->Scale(1./systHistBinned->Integral("width"));

  // Calculate the difference:
  for(Int_t bin = 1; bin <= nominalHistBinned->GetNbinsX(); bin++) {
    Double_t diff = static_cast<Double_t>(nominalHistBinned->GetBinContent(bin)) - static_cast<Double_t>(systHistBinned->GetBinContent(bin));
    difference.push_back(diff);
    LOG(logDEBUG3) << "Hist " << histogram << ": diff bin #" << bin << " " << nominalHistBinned->GetBinContent(bin) << " - " << systHistBinned->GetBinContent(bin) << " = " << diff;
  }

  delete nominalfile;
  delete systematicfile;
  delete binningfile;

  return difference;
}

void extractorDiffXSec::getPredictionUncertainties() {

  m_prediction_errors.clear();
  LOG(logDEBUG2) << "Preparing theory prediction uncertainties.";

  std::vector<Double_t> aDiffXSecMatchUp = calcSampleDifference("Nominal","MATCH_UP","VisGenTTBar1stJetMass");
  std::vector<Double_t> aDiffXSecScaleUp = calcSampleDifference("Nominal","SCALE_UP","VisGenTTBar1stJetMass");

  std::vector<Double_t> aDiffXSecMatchDown = calcSampleDifference("Nominal","MATCH_DOWN","VisGenTTBar1stJetMass");
  std::vector<Double_t> aDiffXSecScaleDown = calcSampleDifference("Nominal","SCALE_DOWN","VisGenTTBar1stJetMass");

  // Take the theory prediction errors into account:
  for(size_t i = 0; i < aDiffXSecMatchUp.size(); ++i) {
    Double_t dxsec_match = (TMath::Abs(aDiffXSecMatchUp.at(i)) + TMath::Abs(aDiffXSecMatchDown.at(i)))/2;
    Double_t dxsec_scale = (TMath::Abs(aDiffXSecScaleUp.at(i)) + TMath::Abs(aDiffXSecScaleDown.at(i)))/2;

    m_prediction_errors.push_back(std::make_pair(dxsec_match,dxsec_scale));
    LOG(logDEBUG) << "Pred. err #" << i+1 << " diffxs = " << dxsec_match << " & " << dxsec_scale;
  }
}

Double_t extractorDiffXSecPrediction::getReco(Int_t bin, Double_t /*mass*/, Double_t reco, Double_t /*bgr*/, Double_t /*ttbgr*/) {
  // Shift the values up/down by the calculated difference:
  LOG(logDEBUG3) << "Bin #" << bin << ": reco_in=" << reco << ", reco_out=" << (reco + m_shiftFactors.at(bin-1));
  return (reco + m_shiftFactors.at(bin));
}

void extractorDiffXSecPrediction::prepareShiftFactors(TString systematic) {

  m_shiftFactors.clear();

  if(systematic.Contains("HAD")) {
    LOG(logDEBUG2) << "Preparing sample difference between POWHEG and MCATNLO.";
    m_shiftFactors = calcSampleDifference("POWHEG","MCATNLO","VisGenTTBar1stJetMass");
  }
  else if(systematic.Contains("CR")) {
    LOG(logDEBUG2) << "Preparing sample difference between PERUGIA11 and PERUGIA11NoCR.";
    m_shiftFactors = calcSampleDifference("PERUGIA11","PERUGIA11NoCR","VisGenTTBar1stJetMass");
  }
  else if(systematic.Contains("UE1")) {
    LOG(logDEBUG2) << "Preparing sample difference between PERUGIA11 and PERUGIA11mpiHi.";
    m_shiftFactors = calcSampleDifference("PERUGIA11","PERUGIA11mpiHi","VisGenTTBar1stJetMass");
  }
  else if(systematic.Contains("UE2")) {
    LOG(logDEBUG2) << "Preparing sample difference between PERUGIA11 and PERUGIA11TeV.";
    m_shiftFactors = calcSampleDifference("PERUGIA11","PERUGIA11TeV","VisGenTTBar1stJetMass");
  }
  else if(systematic.Contains("MATCH") || systematic.Contains("SCALE")) {
    LOG(logDEBUG2) << "Preparing theory prediction uncertainty factors for " << systematic;
    m_shiftFactors = calcSampleDifference("Nominal",systematic,"VisGenTTBar1stJetMass");
  }
  else {
    LOG(logWARNING) << "Systematic " << systematic << " can't be used as theory prediction error.";
    return;
  }

  std::stringstream sv;
  for(std::vector<Double_t>::iterator sf = m_shiftFactors.begin(); sf != m_shiftFactors.end(); ++sf) { sv << *sf << " "; }
  LOG(logDEBUG2) << "Factors for every bin: " << sv.str();
}

