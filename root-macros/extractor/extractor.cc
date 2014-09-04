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

Double_t extractor::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

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


Double_t extractor::getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) {

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

TH1D * extractor::getSignalHistogram(Double_t mass, TFile * histos) {

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
  TH1D * signalHist = new TH1D(type + channel + Form("_m%3.1f",mass),
			       type + channel + Form("_m%3.1f",mass),
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
  }

  if((flags & FLAG_NORMALIZE_YIELD) != 0) {
    signalHist->Sumw2();
    signalHist->Scale(1./signalHist->Integral());
    LOG(logDEBUG) << "Normalized Data hist.";
  }

  // Return signal-only histogram:
  return signalHist;
}

TH1D * extractor::getSimulationHistogram(Double_t mass, TFile * histos) {

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
  TH1D * simulationHist = new TH1D("reco_" + channel + Form("_m%3.1f",mass),
			       "reco_" + channel + Form("_m%3.1f",mass),
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
    simulationHist->Sumw2();
    simulationHist->Scale(1./simulationHist->Integral());
    LOG(logDEBUG) << "Normalized Reco hist.";
  }

  // Return reco histogram:
  return simulationHist;
}

std::vector<TGraphErrors* > extractor::splitBins(TString type, std::vector<Double_t> masses, std::vector<TH1D*> histograms) {

  std::vector<TGraphErrors* > separated_bins;

  // For every bin, prepare a new TGraph containing all mass points:
  for(Int_t bin = 1; bin <= histograms.at(0)->GetNbinsX(); bin++) {

    LOG(logDEBUG2) << "Filling bin " << bin << " mass points ";
 
    // Prepare new TGraphErrors for this bin:
    TGraphErrors * thisbin = new TGraphErrors();
    thisbin->SetTitle(type + "_" + channel + Form("_masses_bin%i",bin));

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
  TString gname = "chi2_" + channel + Form("_bin%i",bin);
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
    TString cname = "chi2_" + channel + Form("_bin%i",bin);
    TString c2name = "inputs_" + channel + Form("_bin%i",bin);
    TString secondname = "input_mc_" + channel + Form("_bin%i",bin);
    TString firstname = "input_dat_" + channel + Form("_bin%i",bin);

    c = new TCanvas(cname,cname);
    c->cd();
    chi2_graph->SetMarkerStyle(20);
    chi2_graph->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2_graph->GetYaxis()->SetTitle("#chi^{2}");
    chi2_graph->Draw("AP");
    DrawDecayChLabel(getChannelLabel(channel),bin);
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

    DrawDecayChLabel(getChannelLabel(channel),bin);
    DrawCMSLabels();
    leg->Draw();

    c->Write(c2name);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(basepath + "/" + c2name + ".pdf"); }
  }

  return chi2_graph;
}

std::pair<TGraphErrors*,TF1*> extractor::getFittedChiSquare(TString ch, std::vector<Double_t> /*masses*/, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc) {

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TGraphErrors * chi2sum = new TGraphErrors();
  TString gname = "chi2_" + ch + "_sum";
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
    TGraphErrors* chi2 = createIntersectionChiSquare(data.at(bin),mc.at(bin),1+bin, dataFit, mcFit);

    dataFits.push_back(dataFit);
    mcFits.push_back(mcFit);

    // Discard insignificant bins if requested:
    Double_t maxVal = TMath::MaxElement(chi2->GetN(),chi2->GetY());
    if(maxVal < chi2significance) {
      LOG(logWARNING) << "Channel " << ch << " bin " << bin+1 << " has low significance: max(chi2) = " << maxVal << " < " << chi2significance;
    }

    // If we don't use or have the Covariance matrix, just use statistical errors:
    if((flags & FLAG_DONT_USE_COVARIANCE) != 0) {
      // Sum them all:
      for(Int_t i = 0; i < chi2->GetN(); i++) {
	Double_t xsum,ysum,x,y;
	chi2->GetPoint(i,x,y);
	chi2sum->GetPoint(i,xsum,ysum);
	LOG(logDEBUG4) << "Adding (" << x << "/" << y << ") to (" << xsum << "/" << ysum << ")";
	chi2sum->SetPoint(i,x,y+ysum);
      }
    }
  }

  // Calculate the overall Chi2 using all bin correlations from covariance matrix:
  if((flags & FLAG_DONT_USE_COVARIANCE) == 0) {
    // Get the inverse covariance matrix:
    TMatrixD * invCov = getInverseCovMatrix(ch,m_sample);
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
  }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + ch + "_sum";
    c = new TCanvas(cname,cname);
    c->cd();
    chi2sum->SetMarkerStyle(20);
    chi2sum->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2sum->GetYaxis()->SetTitle("#chi^{2}");
    chi2sum->Draw("AP");
    DrawDecayChLabel(getChannelLabel(channel));
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

std::pair<TGraphErrors*,TF1*> extractor::getChiSquare(TString ch, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc) {

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TString name = "chi2_" + ch;

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
    DrawDecayChLabel(getChannelLabel(ch));
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

TFile * extractor::selectInputFile(TString sample, TString ch) {
  // Input files for Total Yield mass extraction: preunfolded histograms:
  TString filename = "preunfolded/" + sample + "/" + ch + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TFile * input = TFile::Open(filename,"read");
  if(!input->IsOpen()) {
    LOG(logCRITICAL) << "Failed to access data file " << filename;
    throw 1;
  }
  LOG(logDEBUG) << "Successfully opened file " << filename;
  return input;
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

  DrawDecayChLabel(getChannelLabel(channel));
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

    TFile * datafile = selectInputFile(*sample,channel);
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
  if((flags & FLAG_CHISQUARE_FROM_FITS) == 0) { fit = getChiSquare(channel,masses,data_hists,mc_hists); }
  else {
    LOG(logDEBUG2) << "Splitting histograms into bins...";
    std::vector<TGraphErrors*> data_graphs = splitBins("data",masses,data_hists);
    std::vector<TGraphErrors*> mc_graphs = splitBins("mc",masses,mc_hists);
    fit = getFittedChiSquare(channel,masses,data_graphs,mc_graphs);
  }
  
  LOG(logDEBUG2) << "Minimizing global Chi2 distribution...";
  extractedMass = getMinimum(fit);

  //if(output->IsOpen()) { delete output; }
  return extractedMass;
}

void extractor::setClosureSample(TString closure) {

  // Enable closure test:
  doClosure = true;

  // Input file:
  TString filename = "preunfolded/" + closure + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TFile * closureFile = TFile::Open(filename,"read");

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
    
    Double_t xsecCorrection = getTtbarXsec(mass)/getTtbarXsec(nominalmass);

    // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
    Double_t bgr_other = aBgrHist->GetBinContent(bin) - aTtBgrHist->GetBinContent(bin);

    Double_t pdata = (aRecHist->GetBinContent(bin) + aTtBgrHist->GetBinContent(bin))*xsecCorrection + bgr_other;

    LOG(logDEBUG) << "Closure: Bin #" << bin << " sig=" << aRecHist->GetBinContent(bin) << " pdat=" << pdata;
    // Write pseudo data with background:
    pseudoData->SetBinContent(bin,pdata);
  }

  delete closureFile;
}

extractor::extractor(TString ch, TString sample, uint32_t steeringFlags) : statErrorPos(0), statErrorNeg(0), extractedMass(0), channel(ch), m_sample(sample), samples(), bin_boundaries(), flags(steeringFlags), doClosure(false) {

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

void extractorBackground::prepareScaleFactor(TString systematic) {

  if(systematic.Contains("BG") || systematic.Contains("DY")) {

    if(systematic.Contains("UP")) { scaleFactor = 1.3; }
    if(systematic.Contains("DOWN")) { scaleFactor = 0.7; }
  }
}

TMatrixD * extractor::getInverseCovMatrix(TString ch, TString sample) {

  // Input files for Differential Cross section mass extraction: unfolded distributions
  TString filename = "SVD/" + sample + "/Unfolding_" + ch + "_TtBar_Mass_HypTTBar1stJetMass.root";
  TString histogramname = "SVD_" + ch + "_TtBar_Mass_HypTTBar1stJetMass_" + sample + "_STATCOV";

  TFile * input = TFile::Open(filename,"read");
  if(!input->IsOpen()) {
    LOG(logCRITICAL) << "Failed to access covariance matrix file " << filename;
    throw 1;
  }
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

