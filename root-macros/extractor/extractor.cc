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

std::vector<TGraphErrors* > extractor::splitBins(TString type, std::vector<Double_t> masses, std::vector<TH1D*> histograms) {

  std::vector<TGraphErrors* > separated_bins;

  // For every bin, prepare a new TGraph containing all mass points:
  for(Int_t bin = 1; bin <= histograms.at(0)->GetNbinsX(); bin++) {

    LOG(logDEBUG2) << "Filling bin " << bin << " mass points ";
 
    // Prepare new TGraphErrors for this bin:
    TGraphErrors * thisbin = new TGraphErrors();
    thisbin->SetTitle(type + "_" + m_channel + Form("_masses_bin%i",bin));

    for(std::vector<TH1D*>::iterator hist = histograms.begin(); hist != histograms.end(); ++hist) {
      // Set the corresponding bin in output vector:
      thisbin->SetPoint(hist-histograms.begin(), masses.at(hist-histograms.begin()), (*hist)->GetBinContent(bin));
      thisbin->SetPointError(hist-histograms.begin(),0,(*hist)->GetBinError(bin));
    }
    separated_bins.push_back(thisbin);
  }
  return separated_bins;
}

void extractor::shiftGraph(TGraphErrors* ingraph, Double_t xshift, Double_t yshift) {

  size_t npoints = ingraph->GetN();
  for(size_t i = 0; i < npoints; i++) {
    // Shift by x and y:
    Double_t inX, inY;
    ingraph->GetPoint(i,inX,inY);
    ingraph->SetPoint(i,inX-xshift,inY-yshift);
    //ingraph->SetPointError(i,ingraph->GetErrorX(i),ingraph->GetErrorY(i));
  }
}

TGraphErrors * extractor::createIntersectionChiSquare(TGraphErrors* first, TGraphErrors* second, Int_t bin, TGraphErrors * fit1, TGraphErrors * fit2) {

  TGraphErrors * chi2_graph = new TGraphErrors();
  TString gname = "chi2_" + m_channel + Form("_bin%i",bin);
  chi2_graph->SetTitle(gname);

  size_t n = first->GetN();
  Double_t xmin1,xmax1,ymin1,ymax1;
  first->GetPoint(0,xmin1,ymin1);
  first->GetPoint(n-1,xmax1,ymax1);
  Double_t xmin2,xmax2,ymin2,ymax2;
  second->GetPoint(0,xmin2,ymin2);
  second->GetPoint(n-1,xmax2,ymax2);

  Double_t xmin = std::min(xmin1*0.95,xmin2*0.95);
  Double_t xmax = std::min(xmax1*1.05,xmax2*1.05);
  Double_t ymeana = (ymax1+ymin1)/2;
  Double_t ymeanb = (ymax2+ymin2)/2;

  Double_t xshift = (xmax+xmin)/2;
  Double_t yshift = (ymeana+ymeanb)/2;


  if((flags & FLAG_DONT_SHIFT_GRAPHS) == 0) {
    // Get shifted graphs:
    LOG(logDEBUG) << "Shift all coordinates by [" << xshift << "/" << yshift << "]";
    shiftGraph(first,xshift,yshift);
    shiftGraph(second,xshift,yshift);
    xmin -= xshift; xmax -= xshift;
  }
  // Else use input graphs directly.

  Double_t interval = (xmax-xmin)/(static_cast<Double_t>(granularity));
  std::vector<Double_t> scanPoints(granularity+1,0);
  std::vector<Double_t> confIntervalData(scanPoints.size(),0);  
  std::vector<Double_t> confIntervalMC(scanPoints.size(),0);  

  // Prepare scan point vector:
  for(size_t i = 0; i < scanPoints.size(); i++) {
    scanPoints.at(i) = xmin + static_cast<Double_t>(i)*interval;
  }

  LOG(logDEBUG2) << "Prepared for fitting, " << scanPoints.size() << " scan points in [" << xmin << "," << xmax << "]";

  second->Fit("pol2","Q","",xmin,xmax);
  TF1 * secondFit = second->GetFunction("pol2");
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(scanPoints.size(),1,&scanPoints.at(0),&confIntervalMC.at(0),confidenceLevel);
  TGraphErrors *secondconf = new TGraphErrors(second->GetN());
  // Add additional point, just for drawing:
  secondconf->SetPoint(0,second->GetX()[0] - 1, 0);
  for (Int_t i = 1; i <= second->GetN(); i++) { secondconf->SetPoint(i, second->GetX()[i-1], 0); }
  secondconf->SetPoint(second->GetN(),second->GetX()[second->GetN()-1] + 1, 0);
  // Compute the confidence intervals at the x points of the created graph
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(secondconf);

  first->Fit("pol2","Q","",xmin,xmax);
  TF1 * firstFit = first->GetFunction("pol2");
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(scanPoints.size(),1,&scanPoints.at(0),&confIntervalData.at(0),confidenceLevel);

  // Pick fixed error for data from middle point:
  Double_t fixedError = first->GetErrorY(first->GetN()/2);

  for(size_t i = 0; i < scanPoints.size(); i++) {
    Double_t a = firstFit->EvalPar(&scanPoints.at(i));
    Double_t b = secondFit->EvalPar(&scanPoints.at(i));

    // For data, do not use fitted confidence interval since all points are correlated. Use fixed error:
    Double_t awidth = fixedError;
    // For MC, use confidence interval of the fit - all points are uncorrelated:
    Double_t bwidth = confIntervalMC.at(i);
    Double_t chi2 = chiSquare(b,awidth*awidth+bwidth*bwidth,a);
    
    LOG(logDEBUG4) << "Scan " << i << "@" << scanPoints.at(i) << ": a=" << a << "(" << awidth << "/" << confIntervalData.at(i) << ") b=" << b << "(" << bwidth << ") chi2=" << chi2;
    chi2_graph->SetPoint(i, scanPoints.at(i), chi2);

    fit1->SetPoint(i, scanPoints.at(i), a);
    fit2->SetPoint(i, scanPoints.at(i), b);
  }

  if((flags & FLAG_DONT_SHIFT_GRAPHS) == 0) { 
    // Shift all graphs back to initial position:
    shiftGraph(chi2_graph,-1*xshift,0);
    shiftGraph(secondconf,-1*xshift,-1*yshift);
    shiftGraph(second,-1*xshift,-1*yshift);
    shiftGraph(first,-1*xshift,-1*yshift);
    shiftGraph(fit1,-1*xshift,-1*yshift);
    shiftGraph(fit2,-1*xshift,-1*yshift);
  }

  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TCanvas* c = 0;

    // Naming:
    TString cname = "chi2_" + m_channel + Form("_bin%i",bin);
    TString c2name = "inputs_" + m_channel + Form("_bin%i",bin);
    TString secondname = "input_mc_" + m_channel + Form("_bin%i",bin);
    TString firstname = "input_dat_" + m_channel + Form("_bin%i",bin);

    c = new TCanvas(cname,cname);
    c->cd();
    chi2_graph->SetMarkerStyle(20);
    chi2_graph->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2_graph->GetYaxis()->SetTitle("#chi^{2}");
    chi2_graph->Draw("AP");
    DrawDecayChLabel(getChannelLabel(),bin);
    DrawCMSLabels();
    chi2_graph->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + cname + ".pdf"); }


    c = new TCanvas(c2name,c2name);
    c->cd();
    TLegend *leg = new TLegend();
    setLegendStyle(leg);

    second->SetTitle("");
    second->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    second->GetYaxis()->SetTitle(getQuantity());

    secondconf->SetTitle("");
    secondconf->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    secondconf->GetYaxis()->SetTitle(getQuantity());
    secondconf->GetXaxis()->SetLimits(second->GetX()[0]-1,second->GetX()[second->GetN()-1] + 1);

    setStyle(second,"madgraph");
    setStyleAndFillLegend(secondconf,"madgraph",leg);
    setStyleAndFillLegend(first,"data",leg);

    secondconf->Draw("A E3");
    second->Draw("SAME P");
    second->Write(secondname);

    first->SetPoint(0,first->GetX()[0] - 1,first->GetY()[0]);
    first->SetPoint(first->GetN()-1,first->GetX()[first->GetN()-1] + 1,first->GetY()[first->GetN()-1]);
    first->Draw("SAME L E3");
    first->Write(firstname);

    DrawDecayChLabel(getChannelLabel(),bin);
    DrawCMSLabels();
    leg->Draw();

    c->Write(c2name);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + c2name + ".pdf"); }
  }

  return chi2_graph;
}

std::pair<TGraphErrors*,TF1*> extractor::getFittedChiSquare(std::vector<Double_t> /*masses*/, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc) {

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TGraphErrors * chi2sum = new TGraphErrors();
  TString gname = "chi2_" + m_channel + "_sum";
  chi2sum->SetTitle(gname);

  // Cross-check: we have the same number of rho-S bins:
  if(data.size() != mc.size()) {
    LOG(logCRITICAL) << "Bin numbers don't match!";
    throw(1);
  }

  // Loop over all bins we have:
  for(size_t bin = 0; bin < data.size(); bin++) {

    // Get the likelihood for the two functions:
    TGraphErrors* chi2 = createIntersectionChiSquare(data.at(bin),mc.at(bin),1+bin, NULL, NULL);

    // Discard insignificant bins if requested:
    Double_t maxVal = TMath::MaxElement(chi2->GetN(),chi2->GetY());
    if(maxVal < chi2significance) {
      LOG(logWARNING) << "Channel " << m_channel << " bin " << bin+1 << " has low significance: max(chi2) = " << maxVal << " < " << chi2significance;
    }
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
    DrawDecayChLabel(getChannelLabel());
    DrawCMSLabels();
    chi2sum->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + cname + ".pdf"); }
  }

  // Fit the graph
  chi2sum->Fit("pol2","Q","",170,174.5);
  
  finalChiSquare = std::make_pair(chi2sum,chi2sum->GetFunction("pol2"));
  return finalChiSquare;
}

Double_t extractor::chiSquare(const Double_t center, const Double_t widthsquared, const Double_t eval) {
  return static_cast<Double_t>(static_cast<Double_t>(eval - center)*static_cast<Double_t>(eval - center))/static_cast<Double_t>(widthsquared);
}

std::pair<TGraphErrors*,TF1*> extractor::getChiSquare(std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc) {

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TString name = "chi2_" + m_channel;

  TGraphErrors * chisquare = new TGraphErrors();
  chisquare->SetTitle(name);

  for(UInt_t point = 0; point < masses.size(); ++point) {
    Double_t chi2 = 0;
    // Iterate over all bins:
    for(Int_t bin = 1; bin <= data.at(point)->GetNbinsX(); bin++) {
      chi2 += chiSquare(mc.at(point)->GetBinContent(bin),data.at(point)->GetBinError(bin)*data.at(point)->GetBinError(bin),data.at(point)->GetBinContent(bin));
    }
    chisquare->SetPoint(point, masses.at(point), chi2);
  }

  chisquare->Fit("pol2","Q");

  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TCanvas* c = new TCanvas(name+"_c",name+"_c");
    c->cd();
    chisquare->SetMarkerStyle(20);
    chisquare->GetXaxis()->SetTitle("m_{t} [GeV]");
    chisquare->GetYaxis()->SetTitle("#chi^{2}");
    chisquare->Draw("AP");
    DrawDecayChLabel(getChannelLabel());
    DrawCMSLabels();
    chisquare->Write(name);
    c->Write(name + "_c");
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + name + ".pdf"); }
  }

  finalChiSquare = std::make_pair(chisquare,chisquare->GetFunction("pol2"));
  return finalChiSquare;
}

Double_t extractor::getMinimum(std::pair<TGraphErrors*,TF1*> finalChiSquare) {
  
  Double_t chi2min, x_chi2min, x_left, x_right;

  if((flags & FLAG_RETURN_FITMIN) != 0) {
    // Return the fit function's minimum:
    TF1 * fit = finalChiSquare.second;
    chi2min = fit->GetMinimum(0,330);
    x_chi2min = fit->GetMinimumX(0,330);

    // Statictical error: vary chi2 by +-1 (going up the curve left and right by dChi2 = 1):
    x_left = fit->GetX(chi2min+1,0,x_chi2min);
    x_right = fit->GetX(chi2min+1, x_chi2min,330);
    
    LOG(logDEBUG) << "Minimized by fitting central region w/ pol2. Interval [" << x_left << " - " << x_chi2min << " - " << x_right << "]";
  }
  else {
    // Search the minimum of the TGraph:
    TGraphErrors * graph = finalChiSquare.first;
    double* y = graph->GetY();
    double* x = graph->GetX();
    Int_t locmin = TMath::LocMin(graph->GetN(),y);
    chi2min = y[locmin];
    x_chi2min = x[locmin];
   
    Double_t scanGranularity = 0.005; // GeV
    Double_t scanDistance = 3; // GeV
    // Get left and right bound by scanning the graph and returning the value where Chi2=Chi2Min+1:
    for(x_left = x_chi2min; x_left > x_chi2min - scanDistance; x_left -= scanGranularity) { if(graph->Eval(x_left) >= chi2min+1) break; }
    for(x_right = x_chi2min; x_right < x_chi2min + scanDistance; x_right += scanGranularity) { if(graph->Eval(x_right) >= chi2min+1) break; }

    LOG(logDEBUG) << "Minimized by finding minimum graph point. Interval [" << x_left << " - " << x_chi2min << " - " << x_right << "]";
  }

  // Just take the difference:
  x_left = x_chi2min - x_left;
  x_right = x_right - x_chi2min;

  LOG(logDEBUG) << "Minimum Chi2 is " << chi2min << " at " << x_chi2min << " +" << x_right << "-" << x_left;
  
  // Store the averaged statistical error:
  statErrorPos = x_right;
  statErrorNeg = x_left;
  return x_chi2min;
}

Double_t extractor::getStatError() {

  // Just return the statistical error calculated from chi2:
  return (statErrorPos+statErrorNeg)/2;
}

Double_t extractor::getStatError(Double_t &statPos, Double_t &statNeg) {

  // Just return the statistical error calculated from chi2:
  statPos = statErrorPos;
  statNeg = statErrorNeg;
  return (statErrorPos+statErrorNeg)/2;
}

void extractor::getControlPlots(std::vector<TH1D*> histograms) {

  if((flags & FLAG_STORE_HISTOGRAMS) == 0) { return; }

  TString canvastitle = histograms.at(0)->GetName();
  TCanvas* c = new TCanvas(canvastitle,canvastitle);
  c->cd();

  Int_t colors[7] = {kRed+1,kMagenta+2,kBlue+1,kGreen+2,kBlack,kRed-1,kGreen};
  for(std::vector<TH1D*>::iterator h = histograms.begin(); h != histograms.end(); ++h) {
    (*h)->SetFillStyle(0);
    (*h)->SetLineWidth(2);
    (*h)->SetLineColor(colors[h-histograms.begin()]);

    if(h - histograms.begin() > 0) { (*h)->Draw("HIST SAME"); }
    else {
      (*h)->GetXaxis()->SetTitle("#rho_{s}");
      (*h)->GetYaxis()->SetTitle("Events");
      (*h)->Draw("HIST");
    }
  }

  DrawDecayChLabel(getChannelLabel());
  DrawCMSLabels();
  c->Write(canvastitle);
  if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + canvastitle + ".pdf"); }

  return;
}

Double_t extractor::getTopMass() {

  // Just return mass if extraction has been executed already:
  if(extractedMass > 0) { return extractedMass; }

  std::vector<TH1D*> data_hists;
  std::vector<TH1D*> mc_hists;
  std::vector<Double_t> masses;

  for(std::vector<TString>::iterator sample = samples.begin(); sample != samples.end(); ++sample) {

    TFile * datafile = selectInputFile(*sample);
    Double_t topmass = getMassFromSample(*sample);

    LOG(logDEBUG) << "Top Mass for Sample " << (*sample) << " m_t=" << topmass;

    // Subtract the estimated background from the data:
    TH1D * data = getSignalHistogram(topmass,datafile);
    data->SetDirectory(0);
    // Also get the MC sample:
    TH1D * mc = getSimulationHistogram(topmass,datafile);
    mc->SetDirectory(0);

    masses.push_back(topmass);
    data_hists.push_back(data);
    mc_hists.push_back(mc);

    delete datafile;
  }

  // Store the output histograms into this file:
  TFile * output;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    output = TFile::Open(basepath + "/" + getRootFilename(),"update");
    gDirectory->pwd();
  }

  LOG(logDEBUG2) << "Requesting control plots...";
  getControlPlots(data_hists);
  getControlPlots(mc_hists);

  LOG(logDEBUG2) << "Calculating Chi2 value...";
  std::pair<TGraphErrors*,TF1*> fit;
  if((flags & FLAG_CHISQUARE_FROM_FITS) == 0) { fit = getChiSquare(masses,data_hists,mc_hists); }
  else {
    LOG(logDEBUG2) << "Splitting histograms into bins...";
    std::vector<TGraphErrors*> data_graphs = splitBins("data",masses,data_hists);
    std::vector<TGraphErrors*> mc_graphs = splitBins("mc",masses,mc_hists);
    fit = getFittedChiSquare(masses,data_graphs,mc_graphs);
  }
  
  LOG(logDEBUG2) << "Minimizing global Chi2 distribution...";
  extractedMass = getMinimum(fit);

  //if(output->IsOpen()) { delete output; }
  return extractedMass;
}

extractor::extractor(TString ch, TString sample, uint32_t steeringFlags) : statErrorPos(0), statErrorNeg(0), extractedMass(0), m_channel(ch), m_sample(sample), samples(), bin_boundaries(), flags(steeringFlags), doClosure(false) {

  // Do not add histograms to the directory listing:
  TH1::AddDirectory(kFALSE);

  // Set the histogram styles:
  setHHStyle(*gStyle);

  // This is our nominal mass variation sample:
  if(sample == "Nominal") {
    samples.push_back("MASS_DOWN_6GEV");
    samples.push_back("MASS_DOWN_3GEV");
    samples.push_back("MASS_DOWN_1GEV");
    samples.push_back("Nominal");
    samples.push_back("MASS_UP_1GEV");
    samples.push_back("MASS_UP_3GEV");
    samples.push_back("MASS_UP_6GEV");
  }
  // We are looking at a systematic uncertainty sample here:
  else {
    samples.push_back(sample+"_6NEG");
    samples.push_back(sample+"_3NEG");
    samples.push_back(sample+"_1NEG");
    samples.push_back(sample);
    samples.push_back(sample+"_1POS");
    samples.push_back(sample+"_3POS");
    samples.push_back(sample+"_6POS");
  }

  if((flags & FLAG_NORMALIZE_YIELD) != 0
     && (flags & FLAG_LASTBIN_EXTRACTION) != 0) {
    LOG(logERROR) << "Normalization of a single bin doesn't make any sense. Dropping "
		  << "NORMALIZE_YIELD in favor for extracting from last bin only.";
    flags &= ~FLAG_NORMALIZE_YIELD;
  }

  std::stringstream s;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0 ) { s << "STORE_HISTOGRAMS "; }
  if((flags & FLAG_NORMALIZE_YIELD) != 0 ) { s << "NORMALIZE_YIELD "; }
  if((flags & FLAG_LASTBIN_EXTRACTION) != 0 ) { s << "LASTBIN_EXTRACTION "; }
  if((flags & FLAG_UNFOLD_ALLMASSES) != 0 ) { s << "UNFOLD_ALLMASSES "; }
  if((flags & FLAG_CHISQUARE_FROM_FITS) != 0 ) { s << "CHISQUARE_FROM_FITS "; }
  if((flags & FLAG_DONT_SHIFT_GRAPHS) != 0 ) { s << "DONT_SHIFT_GRAPHS "; }
  if((flags & FLAG_STORE_PDFS) != 0 ) { s << "STORE_PDFS "; }
  if((flags & FLAG_RETURN_FITMIN) != 0 ) { s << "RETURN_FITMIN "; }
  if((flags & FLAG_DONT_USE_COVARIANCE) != 0 ) { s << "DONT_USE_COVARIANCE "; }

  LOG(logDEBUG) << "Initialized. Flags shipped: " << s.str();
}


