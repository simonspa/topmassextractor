#include <iostream>
#include <vector>
#include <string>
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
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
#include "extract.h"

Double_t nominalmass = 172.5;
Double_t lumi = 19712;

Int_t granularity = 500;
Double_t confidenceLevel = 0.95;

using namespace unilog;


Double_t extractorOtherSamples::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  // variation for every mass sample in the current bin:
  reco -= deltaRec.at(bin-1);
  bgr -= deltaBgr.at(bin-1);
  ttbgr -= deltaTtbgr.at(bin-1);

  // Call parent class signal calculation function:
  return extractor::getSignal(bin, mass, data, reco, bgr, ttbgr);
}

Double_t extractor::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Calculate the signal fraction from reconstructed events and TT background:
  Double_t fsignal = reco/(reco+ttbgr);

  // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
  Double_t bgr_other = bgr - ttbgr;

  // Calculate signal by subtracting backround from data, multiplied by signal fraction.
  Double_t signal = (data - bgr_other)*fsignal;

  LOG(logDEBUG2) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

  return signal;
}

Double_t extractorBackground::getSignal(Int_t bin, Double_t /*mass*/, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Calculate the signal fraction from reconstructed events and TT background:
  Double_t fsignal = reco/(reco+ttbgr);

  // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
  Double_t bgr_other = bgr - ttbgr;

  // Scale with DY or BG scale factor:
  bgr_other *= scaleFactor;

  // Calculate signal by subtracting backround from data, multiplied by signal fraction.
  Double_t signal = (data - bgr_other)*fsignal;

  LOG(logDEBUG2) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

  return signal;
}

Double_t extractor::getReco(Int_t bin, Double_t mass, Double_t reco) {

  // Scale the reco accoring to the different TTBar Cross sections (mass dependent):
  Double_t corr_reco = reco*getTtbarXsec(mass)/getTtbarXsec(nominalmass);
  LOG(logDEBUG2) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

  // Return the reco event count corrected by the ttbar Xsec at given mass:
  return corr_reco;
}

Double_t extractorOtherSamples::getReco(Int_t bin, Double_t mass, Double_t reco) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  //variation for every bin:
  reco -= deltaRec.at(bin-1);

  // Call parent class signal calculation function:
  return extractor::getReco(bin, mass, reco);
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
  for (Int_t bin = startbin; bin <= nbins; bin++) Xbins[bin-startbin] = aRecHist->GetBinLowEdge(bin);
  Xbins[nbins+1-startbin] = aRecHist->GetBinLowEdge(nbins) + aRecHist->GetBinWidth(nbins);
  TH1D * signalHist = new TH1D(type + channel + Form("_m%3.1f",mass),
			       type + channel + Form("_m%3.1f",mass),
			       nbins-startbin+1,
			       Xbins);

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
    Double_t corr_reco = getReco(bin, mass,aRecHist->GetBinContent(bin));

    // Write corrected Reco:
    simulationHist->SetBinContent(bin-startbin+1,corr_reco);
  }

  if((flags & FLAG_NORMALIZE_YIELD) != 0) {
    simulationHist->Sumw2();
    simulationHist->Scale(1./simulationHist->Integral());
    LOG(logDEBUG) << "Normalized Reco hist.";
  }

  // Return reco histogram:
  return simulationHist;
}

std::vector<TH1D* > extractor::splitBins(TString type, std::vector<TH1D*> histograms) {

  std::vector<TH1D* > separated_bins;
  Int_t nbins = histograms.at(0)->GetNbinsX();

  // For every bin, prepare a new histogram containing all mass points:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    LOG(logDEBUG) << "Filling bin " << bin << " mass points ";
 
    // Prepare new vector for this bin:
    TString hname = type + "_" + channel + Form("_masses_bin%i",bin);
    TH1D* thisbin = new TH1D(hname,hname,histograms.size(),0,histograms.size());

    Int_t newbin = 1;
    for(std::vector<TH1D*>::iterator hist = histograms.begin(); hist != histograms.end(); ++hist) {
      // Set the corresponding bin in output vector:
      thisbin->SetBinContent(newbin,(*hist)->GetBinContent(bin));
      thisbin->SetBinError(newbin,(*hist)->GetBinError(bin));
      newbin++;
    }

    separated_bins.push_back(thisbin);
  }
  
  return separated_bins;
}

std::pair<TGraphErrors*,TGraphErrors*> extractor::fitMassBins(TString ch, Int_t bin, std::vector<Double_t> masses, TH1D* data, TH1D* mc) {

  std::pair<TGraphErrors*,TGraphErrors*> allfits;

  TString mname, dname, cname;
  mname.Form("mc_%i_",bin);
  dname.Form("dat_%i_",bin);
  cname.Form("dat_mc_%i_",bin);

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    c = new TCanvas(cname+ch,cname+ch);
    c->cd();
  }

  TGraphErrors * graph_mc = new TGraphErrors();
  graph_mc->SetTitle(mname+ch);

  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph_mc->SetPoint(point, masses.at(point), mc->GetBinContent(point+1));
    graph_mc->SetPointError(point,0,mc->GetBinError(point+1));
  }

  TLegend *leg = new TLegend();
  setLegendStyle(leg);

  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    graph_mc->SetTitle("");
    graph_mc->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    graph_mc->GetYaxis()->SetTitle(getQuantity());

    setStyleAndFillLegend(graph_mc,"madgraph",leg);
    graph_mc->Draw("A P");
    graph_mc->Write(mname+ch);
  }

  TGraphErrors * graph = new TGraphErrors();
  graph->SetTitle(dname+ch);

  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph->SetPoint(point, masses.at(point), data->GetBinContent(point+1));
    graph->SetPointError(point,0,data->GetBinError(point+1));
  }

  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    setStyleAndFillLegend(graph,"data",leg);

    graph->Draw("SAME P E1");
    DrawDecayChLabel(getChannelLabel(ch));
    DrawCMSLabels();

    // Also draw legend:
    leg->Draw();

    graph->Write(dname+ch);
    c->Write();
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(cname + ch + ".pdf"); }
  }

  allfits = std::make_pair(graph,graph_mc);
  return allfits;
}

TGraphErrors * extractor::getShiftedGraph(TGraphErrors* ingraph, Double_t xshift, Double_t yshift) {

  TGraphErrors * outgraph = new TGraphErrors();
  
  size_t npoints = ingraph->GetN();
  for(size_t i = 0; i < npoints; i++) {
    // Shift by x and y:
    Double_t inX, inY;
    ingraph->GetPoint(i,inX,inY);
    outgraph->SetPoint(i,inX-xshift,inY-yshift);
    outgraph->SetPointError(i,ingraph->GetErrorX(i),ingraph->GetErrorY(i));
  }

  return outgraph;
}

TGraphErrors * extractor::createIntersectionChiSquare(std::pair<TGraphErrors*,TGraphErrors*> fits, Int_t bin) {

  TGraphErrors * chi2_graph = new TGraphErrors();
  TGraphErrors * chi2_tmp = new TGraphErrors();
  
  chi2_graph->SetTitle("chi2_" + channel + Form("_bin%i",bin));

  size_t n = fits.first->GetN();
  Double_t xmin1,xmax1,ymin1,ymax1;
  fits.first->GetPoint(0,xmin1,ymin1);
  fits.first->GetPoint(n-1,xmax1,ymax1);
  Double_t xmin2,xmax2,ymin2,ymax2;
  fits.second->GetPoint(0,xmin2,ymin2);
  fits.second->GetPoint(n-1,xmax2,ymax2);

  Double_t xmin = std::min(xmin1*0.95,xmin2*0.95);
  Double_t xmax = std::min(xmax1*1.05,xmax2*1.05);
  Double_t ymeana = (ymax1+ymin1)/2;
  Double_t ymeanb = (ymax2+ymin2)/2;

  Double_t xshift = (xmax+xmin)/2;
  Double_t yshift = (ymeana+ymeanb)/2;

  // Shift!
  TGraphErrors * first, * second;
  if((flags & FLAG_DONT_SHIFT_GRAPHS) != 0) {
    // Use input graphs directly:
    first = fits.first;
    second = fits.second;
  }
  else {
    // Get shifted graphs:
    LOG(logDEBUG) << "Shift all coordinates by [" << xshift << "/" << yshift << "]";
    first = getShiftedGraph(fits.first,xshift,yshift);
    second = getShiftedGraph(fits.second,xshift,yshift);
    xmin -= xshift; xmax -= xshift;
  }

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
    // Double_t awidth = confIntervalData.at(i);
    // For MC, use confidence interval of the fit - all points are uncorrelated:
    Double_t bwidth = confIntervalMC.at(i);
    Double_t chi2 = chiSquare(b,awidth*awidth+bwidth*bwidth,a);
    
    LOG(logDEBUG3) << "Scan " << i << "@" << scanPoints.at(i) << ": a=" << a << "(" << awidth << "/" << confIntervalData.at(i) << ") b=" << b << "(" << bwidth << ") chi2=" << chi2;
    chi2_tmp->SetPoint(i, scanPoints.at(i), chi2);
  }

  if((flags & FLAG_DONT_SHIFT_GRAPHS) == 0) { chi2_graph = getShiftedGraph(chi2_tmp,-1*xshift,0); }
  else { chi2_graph = chi2_tmp; }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + channel + Form("_bin%i",bin);
    c = new TCanvas(cname,cname);
    c->cd();
    chi2_graph->SetMarkerStyle(20);
    chi2_graph->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2_graph->GetYaxis()->SetTitle("#chi^{2}");
    chi2_graph->Draw("AP");
    DrawDecayChLabel(getChannelLabel(channel));
    DrawCMSLabels();
    chi2_graph->Write();
    c->Write();
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(cname + ".pdf"); }
  }

  return chi2_graph;
}

std::pair<TGraphErrors*,TF1*> extractor::getFittedChiSquare(TString ch, std::vector<Double_t> masses, std::vector<std::pair<TGraphErrors*,TGraphErrors*> > fits) {

  std::pair<TGraphErrors*,TF1*> finalChiSquare;
  TGraphErrors * chi2sum = new TGraphErrors();
  chi2sum->SetTitle("chi2_" + channel + "_sum");

  // Loop over all bins we have:
  for(std::vector<std::pair<TGraphErrors*,TGraphErrors*> >::iterator binfits = fits.begin(); binfits != fits.end(); binfits++) {
    // Get the likelihood for the two functions:
    TGraphErrors* chi2 = createIntersectionChiSquare(*binfits,1+binfits-fits.begin());

    // Sum them all:
    for(size_t i = 0; i < chi2->GetN(); i++) {
      Double_t xsum,ysum,x,y;
      chi2->GetPoint(i,x,y);
      chi2sum->GetPoint(i,xsum,ysum);
      LOG(logDEBUG3) << "Adding (" << x << "/" << y << ") to (" << xsum << "/" << ysum << ")";
      chi2sum->SetPoint(i,x,y+ysum);
    }
  }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + channel + "_sum";
    c = new TCanvas(cname,cname);
    c->cd();
    chi2sum->SetMarkerStyle(20);
    chi2sum->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2sum->GetYaxis()->SetTitle("#chi^{2}");
    chi2sum->Draw("AP");
    DrawDecayChLabel(getChannelLabel(channel));
    DrawCMSLabels();
    chi2sum->Write();
    c->Write();
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(cname + ".pdf"); }
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

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    c = new TCanvas(name+"_c",name+"_c");
    c->cd();
  }

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
    chisquare->SetMarkerStyle(20);
    chisquare->GetXaxis()->SetTitle("m_{t} [GeV]");
    chisquare->GetYaxis()->SetTitle("#chi^{2}");
    chisquare->Draw("AP");
    DrawDecayChLabel(getChannelLabel(ch));
    DrawCMSLabels();
    chisquare->Write(name);
    c->Write();
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(name + ".pdf"); }
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
   
    Double_t scanGranularity = 0.01; // GeV
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
  c->Write();
  if((flags & FLAG_STORE_PDFS) != 0) { c->Print(canvastitle + ".pdf"); }

  return;
}

Double_t extractor::getTopMass() {

  // Do not add histograms to the directory listing:
  TH1::AddDirectory(kFALSE);

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

  std::vector<TH1D*> separated_data = splitBins("data",data_hists);
  std::vector<TH1D*> separated_mc = splitBins("mc",mc_hists);

  TFile * output;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    output = TFile::Open("MassFitRates.root","update");
    gDirectory->pwd();
  }

  getControlPlots(data_hists);
  getControlPlots(mc_hists);

  std::vector<std::pair<TGraphErrors*,TGraphErrors*> > fits;
  for(UInt_t bin = 0; bin < separated_data.size(); ++bin) {
    fits.push_back(fitMassBins(channel,bin+1,masses,separated_data.at(bin),separated_mc.at(bin)));
  }

  std::pair<TGraphErrors*,TF1*> fit;
  if((flags & FLAG_CHISQUARE_FROM_FITS) == 0) { fit = getChiSquare(channel,masses,data_hists,mc_hists); }
  else { fit = getFittedChiSquare(channel,masses,fits); }
  
  extractedMass = getMinimum(fit);

  //if(output->IsOpen()) { delete output; }
  return extractedMass;
}

void extractorOtherSamples::calcDifferenceToNominal(TString nominal, TString systematic) {

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

    LOG(logDEBUG2) << "Diff bin #" << bin << " reco: " << nominalReco->GetBinContent(bin) << " - " << varReco->GetBinContent(bin) << " = " << rec << " dbgr=" << bgr << " dttbgr=" << ttbgr;
  }

  delete nominalfile;
  delete systematicfile;
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

void extractorDiffXSec::setClosureSample(TString /*closure*/) {
  LOG(logERROR) << "Can't set closure sample - use unfolding step to prepare closure data!";
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

extractor::extractor(TString ch, TString sample, uint32_t steeringFlags) : statErrorPos(0), statErrorNeg(0), extractedMass(0), channel(ch), samples(), flags(steeringFlags), doClosure(false) {

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

  LOG(logINFO) << "Flags shipped: " << s.str();
  LOG(logINFO) << "Initialized.";
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

void extractor::setHHStyle(TStyle& HHStyle)
{
    const int fontstyle=42;
    HHStyle.SetPalette(1);
        
    // ==============
    //  Canvas
    // ==============
            
    HHStyle.SetCanvasBorderMode(0);
    HHStyle.SetCanvasColor(kWhite);
    HHStyle.SetCanvasDefH(600); //Height of canvas
    HHStyle.SetCanvasDefW(600); //Width of canvas
    HHStyle.SetCanvasDefX(0);   //Position on screen
    HHStyle.SetCanvasDefY(0);
            
    // ==============
    //  Pad
    // ==============
            
    HHStyle.SetPadBorderMode(0);
    // HHStyle.SetPadBorderSize(Width_t size = 1);
    HHStyle.SetPadColor(kWhite);
    HHStyle.SetPadGridX(false);
    HHStyle.SetPadGridY(false);
    HHStyle.SetGridColor(kGray);
    HHStyle.SetGridStyle(3);
    HHStyle.SetGridWidth(1);
            
    // ==============
    //  Frame
    // ==============
            
    HHStyle.SetFrameBorderMode(0);
    HHStyle.SetFrameBorderSize(1);
    HHStyle.SetFrameFillColor(0);
    HHStyle.SetFrameFillStyle(0);
    HHStyle.SetFrameLineColor(1);
    HHStyle.SetFrameLineStyle(1);
    HHStyle.SetFrameLineWidth(1);
            
    // ==============
    //  Histo
    // ==============

    HHStyle.SetErrorX(0.0);
    HHStyle.SetEndErrorSize(8);
            
    // HHStyle.SetHistFillColor(1);
    // HHStyle.SetHistFillStyle(0);
    // HHStyle.SetHistLineColor(1);
    HHStyle.SetHistLineStyle(0);
    HHStyle.SetHistLineWidth(1);
    // HHStyle.SetLegoInnerR(Float_t rad = 0.5);
    // HHStyle.SetNumberContours(Int_t number = 20);

    // HHStyle.SetErrorMarker(20);
            
    HHStyle.SetMarkerStyle(20);
            
    // ==============
    //  Fit/function
    // ==============
            
    HHStyle.SetOptFit(0);
    HHStyle.SetFitFormat("5.4g");
    HHStyle.SetFuncColor(2);
    HHStyle.SetFuncStyle(1);
    HHStyle.SetFuncWidth(1);
            
    // ==============
    //  Date
    // ============== 
            
    HHStyle.SetOptDate(0);
    // HHStyle.SetDateX(Float_t x = 0.01);
    // HHStyle.SetDateY(Float_t y = 0.01);
            
    // =====================
    //  Statistics Box
    // =====================
            
    HHStyle.SetOptFile(0);
    HHStyle.SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
    HHStyle.SetStatColor(kWhite);
    HHStyle.SetStatFont(fontstyle);
    HHStyle.SetStatFontSize(0.025);
    HHStyle.SetStatTextColor(1);
    HHStyle.SetStatFormat("6.4g");
    HHStyle.SetStatBorderSize(1);
    HHStyle.SetStatH(0.1);
    HHStyle.SetStatW(0.15);
    // HHStyle.SetStatStyle(Style_t style = 1001);
    // HHStyle.SetStatX(Float_t x = 0);
    // HHStyle.SetStatY(Float_t y = 0);
            
    // ==============
    //  Margins
    // ==============

    HHStyle.SetPadTopMargin(0.1);
    HHStyle.SetPadBottomMargin(0.15);
    HHStyle.SetPadLeftMargin(0.20);
    HHStyle.SetPadRightMargin(0.05);
            
    // ==============
    //  Global Title
    // ==============
            
    HHStyle.SetOptTitle(0);
    HHStyle.SetTitleFont(fontstyle);
    HHStyle.SetTitleColor(1);
    HHStyle.SetTitleTextColor(1);
    HHStyle.SetTitleFillColor(10);
    HHStyle.SetTitleFontSize(0.05);
    // HHStyle.SetTitleH(0); // Set the height of the title box
    // HHStyle.SetTitleW(0); // Set the width of the title box
    // HHStyle.SetTitleX(0); // Set the position of the title box
    // HHStyle.SetTitleY(0.985); // Set the position of the title box
    // HHStyle.SetTitleStyle(Style_t style = 1001);
    // HHStyle.SetTitleBorderSize(2);
            
    // ==============
    //  Axis titles
    // ==============
            
    HHStyle.SetTitleColor(1, "XYZ");
    HHStyle.SetTitleFont(fontstyle, "XYZ");
    HHStyle.SetTitleSize(0.05, "XYZ");
    // HHStyle.SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // HHStyle.SetTitleYSize(Float_t size = 0.02);
    HHStyle.SetTitleXOffset(1.0);
    HHStyle.SetTitleYOffset(1.7);
    // HHStyle.SetTitleOffset(1.1, "Y"); // Another way to set the Offset
            
    // ==============
    //  Axis Label
    // ==============
            
    //HHStyle.SetLabelColor(1, "XYZ");
    HHStyle.SetLabelFont(fontstyle, "XYZ");
    HHStyle.SetLabelOffset(0.007, "XYZ");
    HHStyle.SetLabelSize(0.04, "XYZ");
            
    // ==============
    //  Axis
    // ==============
            
    HHStyle.SetAxisColor(1, "XYZ");
    HHStyle.SetStripDecimals(kTRUE);
    HHStyle.SetTickLength(0.03, "XYZ");
    HHStyle.SetNdivisions(510, "XYZ");
    HHStyle.SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    HHStyle.SetPadTickY(1);
            
    // Change for log plots:
    HHStyle.SetOptLogx(0);
    HHStyle.SetOptLogy(0);
    HHStyle.SetOptLogz(0);
            
    // ==============
    //  Text
    // ==============
            
    HHStyle.SetTextAlign(11);
    HHStyle.SetTextAngle(0);
    HHStyle.SetTextColor(1);
    HHStyle.SetTextFont(fontstyle);
    HHStyle.SetTextSize(0.05);
            
    // =====================
    //  Postscript options:
    // =====================
            
    HHStyle.SetPaperSize(20.,20.);
    // HHStyle.SetLineScalePS(Float_t scale = 3);
    // HHStyle.SetLineStyleString(Int_t i, const char* text);
    // HHStyle.SetHeaderPS(const char* header);
    // HHStyle.SetTitlePS(const char* pstitle);
            
    // HHStyle.SetBarOffset(Float_t baroff = 0.5);
    // HHStyle.SetBarWidth(Float_t barwidth = 0.5);
    // HHStyle.SetPaintTextFormat(const char* format = "g");
    // HHStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    // HHStyle.SetTimeOffset(Double_t toffset);
    // HHStyle.SetHistMinimumZero(kTRUE);
}


// Draw label for Decay Channel in upper left corner of plot
void extractor::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void extractor::DrawCMSLabels(int cmsprelim, double energy, double textSize) {

    const char *text;
    if(cmsprelim ==2 ) {//Private work for PhDs students
        text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else if (cmsprelim==1) {//CMS preliminary label
        text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else {//CMS label
        text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    }
    
    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label->SetY2NDC(1.0);
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}

TString extractor::getChannelLabel(TString ch) {
 
  TString label = "";
  if(ch =="ee") { label = "ee"; }
  if(ch =="mumu"){ label = "#mu#mu"; }
  if(ch =="emu"){ label = "e#mu"; }
  if(ch =="combined"){ label = "Dilepton Combined"; }

  return label;
}

TString extractor::getSampleLabel(TString systematic) {

  TString label = "";
  if(systematic.Contains("Nominal") || systematic.Contains("MASS")) { label = "Nominal"; }
  else if(systematic.Contains("BTAG")) { label = "B-Tagging"; }
  else if(systematic.Contains("JER")) { label = "Jet Energy Resolution"; }
  else if(systematic.Contains("JES")) { label = "Jet Energy Scale"; }
  else if(systematic.Contains("PU")) { label = "Pile-Up"; }
  else if(systematic.Contains("TRIG")) { label = "Trigger"; }
  else if(systematic.Contains("LEPT")) { label = "Lepton"; }
  else if(systematic.Contains("BG")) { label = "Background"; }
  else if(systematic.Contains("DY")) { label = "Drell-Yan"; }
  else if(systematic.Contains("KIN")) { label = "Kinematic Reconstruction"; }
  else if(systematic.Contains("MATCH")) { label = "Matching"; }
  else if(systematic.Contains("SCALE")) { label = "Scale"; }
  else if(systematic.Contains("HAD")) { label = "Model"; } //{ label = "Hadronization"; }
  else if(systematic.Contains("MASS")) { label = "Top Mass"; }
  else if(systematic.Contains("CR")) { label = "Color Reconnection"; }
  else if(systematic.Contains("UE")) { label = "Underlying Event"; }

  return label;
}

void extractor::setStyle(TGraphErrors *hist)
{
  hist->SetLineWidth(1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.7);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelOffset(0.007);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
  //hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}"+" #left[GeV^{-1}#right]"); 
}


void extractor::setStyleAndFillLegend(TGraphErrors* hist, TString name, TLegend *leg) {

  setStyle(hist);

  if(name == "data"){
    hist->SetLineWidth(0);
    if(leg && doClosure) leg->AddEntry(hist, "Pseudo Data",  "p");
    else if(leg) leg->AddEntry(hist, "Data",  "p");
  }

  /*  if(name != "data") {
    hist->SetMarkerStyle(1);
    hist->SetMarkerSize(0);
    }*/
  
  if(name == "madgraph") {
    hist->SetMarkerColor(kRed+1);
    hist->SetLineColor(kRed+1);
    if(leg) leg->AddEntry(hist, "MadGraph+Pythia",  "p");
  }
}

void extractor::setLegendStyle(TLegend *leg)
{
    double x1 = 0.560, y1 = 0.755;
    double height = 0.075, width = 0.275;

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}

void extractorBackground::prepareScaleFactor(TString systematic) {

  if(systematic.Contains("BG") || systematic.Contains("DY")) {

    if(systematic.Contains("UP")) { scaleFactor = 1.3; }
    if(systematic.Contains("DOWN")) { scaleFactor = 0.7; }
  }
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

  for(Int_t bin = startbin; bin <= nbins; bin++) {
    LOG(logDEBUG2) << "Bin #" << bin << ": data=" << aDiffXSecHist->GetBinContent(bin);
    signalHist->SetBinContent(bin+1-startbin,aDiffXSecHist->GetBinContent(bin));
    signalHist->SetBinError(bin+1-startbin,aDiffXSecHist->GetBinError(bin));
  }

  // Return DiffXSec signal histogram:
  return signalHist;
}

TH1D * extractorDiffXSec::getSimulationHistogram(Double_t mass, TFile * histos) {

  std::vector<TString> filenames;
  TString generator = "MADGRAPH", filename;

  TString sample;

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

  aMcHist->Scale(1./aMcHist->Integral("width"));
  aMcBinned = dynamic_cast<TH1D*>(aMcHist->Rebin(nbins,"madgraphplot",Xbins));

  for (Int_t bin=0; bin < nbins; bin++) {
    // Condense matrices to arrays for plotting
    aMcBinned->SetBinContent(bin+1,aMcBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/aMcHist->GetBinWidth(1)));
  }
  aMcBinned->Sumw2();
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
    LOG(logDEBUG2) << "Bin #" << bin << ": reco=" << aMcBinned->GetBinContent(bin);
    simulationHist->SetBinContent(bin+1-startbin,aMcBinned->GetBinContent(bin));
    simulationHist->SetBinError(bin+1-startbin,aMcBinned->GetBinError(bin));
  }

  LOG(logDEBUG) << "Returning Simulation histogram now.";
  return simulationHist;
}

TFile * extractorDiffXSec::selectInputFile(TString sample, TString ch) {

  // Overwrite the samples unfolded with different masses with just the nominal:
  if((flags & FLAG_UNFOLD_ALLMASSES) == 0 ) { 
    LOG(logDEBUG2) << "Overwriting: " << sample;
    // Mass samples from Nominal:
    if(sample.Contains("GEV")) { sample = "Nominal"; }
    // All other mass-varied samples:
    else if(sample.Contains("POS") || sample.Contains("NEG")) { sample.Remove(sample.Length()-5); }
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



void extract_yield(std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags) {

  // #######################################
  // ###              YIELD              ###
  // #######################################

  // Do the same (or similar things) for the systematic uncertainties:
  std::vector<TString> syst_on_nominal;
  syst_on_nominal.push_back("MATCH_UP");
  syst_on_nominal.push_back("MATCH_DOWN");
  syst_on_nominal.push_back("SCALE_UP");
  syst_on_nominal.push_back("SCALE_DOWN");
  syst_on_nominal.push_back("HAD_UP");
  syst_on_nominal.push_back("HAD_DOWN");
  syst_on_nominal.push_back("CR_UP");
  syst_on_nominal.push_back("CR_DOWN");
  syst_on_nominal.push_back("UE_UP");
  syst_on_nominal.push_back("UE_DOWN");
  
  std::vector<TString> syst_bg;
  syst_bg.push_back("BG_UP");
  syst_bg.push_back("BG_DOWN");

  std::vector<TString> systematics;
  systematics.push_back("JES_UP"); systematics.push_back("JES_DOWN");
  systematics.push_back("JER_UP"); systematics.push_back("JER_DOWN");
  systematics.push_back("PU_UP"); systematics.push_back("PU_DOWN");
  systematics.push_back("TRIG_UP"); systematics.push_back("TRIG_DOWN");
  systematics.push_back("KIN_UP"); systematics.push_back("KIN_DOWN");
  systematics.push_back("LEPT_UP"); systematics.push_back("LEPT_DOWN");
  systematics.push_back("BTAG_UP"); systematics.push_back("BTAG_DOWN");
  systematics.push_back("BTAG_LJET_UP"); systematics.push_back("BTAG_LJET_DOWN");
  systematics.push_back("BTAG_PT_UP"); systematics.push_back("BTAG_PT_DOWN");
  systematics.push_back("BTAG_ETA_UP"); systematics.push_back("BTAG_ETA_DOWN");
  systematics.push_back("BTAG_LJET_PT_UP"); systematics.push_back("BTAG_LJET_PT_DOWN");
  systematics.push_back("BTAG_LJET_ETA_UP"); systematics.push_back("BTAG_LJET_ETA_DOWN");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream SystOutputFile("MassFitRatesSystematics_" + *ch + ".txt", std::ofstream::trunc);
    SystOutputFile << "Top Mass, Channel: " << *ch << endl;
    SystOutputFile << "Systematic & Syst. error on m_t & [GeV] \\\\" << endl;
    SystOutputFile << "\\hline" << std::endl;


    extractor * mass_samples = new extractor(*ch,"Nominal", flags | FLAG_STORE_HISTOGRAMS);
    if(closure) mass_samples->setClosureSample(closure_sample);

    Double_t topmass = mass_samples->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_samples->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;
    delete mass_samples;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;

    // Systematic Variations with own samples:
    for(std::vector<TString>::iterator syst = syst_on_nominal.begin(); syst != syst_on_nominal.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorOtherSamples * matchscale_samples = new extractorOtherSamples(*ch,"Nominal", flags,(*syst));
      if(closure) matchscale_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = matchscale_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("UP")) SystOutputFile << matchscale_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(3) << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      delete matchscale_samples;
    }

    // Getting DrellYan/Background variations:
    for(std::vector<TString>::iterator syst = syst_bg.begin(); syst != syst_bg.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorBackground * bg_samples = new extractorBackground(*ch,"Nominal",flags,(*syst));
      if(closure) bg_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = bg_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;
      
      if(syst->Contains("UP")) SystOutputFile << bg_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(3) << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
    }

    // Systematic Variations produced by varying nominal samples:
    Double_t btag_syst_pos = 0, btag_syst_neg = 0;
    for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractor * variation_samples = new extractor(*ch,*syst,flags);
      if(closure) variation_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = variation_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("BTAG")) {
	if(delta > 0) btag_syst_pos += delta*delta;
	else btag_syst_neg += delta*delta;
      }
      else {
	if(syst->Contains("UP")) SystOutputFile << variation_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}_{";
	else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      }
    }
    SystOutputFile << mass_samples->getSampleLabel("BTAG") << " & $^{+" << setprecision(3) <<  sqrt(btag_syst_pos) << "}_{"
		   << setprecision(3) << sqrt(btag_syst_neg) << "}$ \\\\" << endl;

    total_syst_pos = sqrt(total_syst_pos);
    total_syst_neg = sqrt(total_syst_neg);
    LOG(logRESULT) << "Channel " << *ch << ": m_t = " << setprecision(6) << topmass << setprecision(3) << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV";

    SystOutputFile << "Channel " << *ch << ": m_t = " << setprecision(6) << topmass << setprecision(3) << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV" << endl;
    SystOutputFile.close();
  }
  return;
}


void extract_diffxsec(std::vector<TString> channels, uint32_t flags) {

  std::vector<TString> systematics;
  systematics.push_back("MATCH_UP"); systematics.push_back("MATCH_DOWN");
  systematics.push_back("SCALE_UP"); systematics.push_back("SCALE_DOWN");
  // If we only take one mass for unfolding, we need to add this as systematic error:
  if((flags & FLAG_UNFOLD_ALLMASSES) == 0) {
    systematics.push_back("MASS_UP"); systematics.push_back("MASS_DOWN");
  }
  systematics.push_back("BG_UP"); systematics.push_back("BG_DOWN");
  systematics.push_back("JES_UP"); systematics.push_back("JES_DOWN");
  systematics.push_back("JER_UP"); systematics.push_back("JER_DOWN");
  systematics.push_back("PU_UP"); systematics.push_back("PU_DOWN");
  systematics.push_back("TRIG_UP"); systematics.push_back("TRIG_DOWN");
  systematics.push_back("KIN_UP"); systematics.push_back("KIN_DOWN");
  systematics.push_back("LEPT_UP"); systematics.push_back("LEPT_DOWN");
  systematics.push_back("BTAG_UP"); systematics.push_back("BTAG_DOWN");
  systematics.push_back("BTAG_LJET_UP"); systematics.push_back("BTAG_LJET_DOWN");
  systematics.push_back("BTAG_PT_UP"); systematics.push_back("BTAG_PT_DOWN");
  systematics.push_back("BTAG_ETA_UP"); systematics.push_back("BTAG_ETA_DOWN");
  systematics.push_back("BTAG_LJET_PT_UP"); systematics.push_back("BTAG_LJET_PT_DOWN");
  systematics.push_back("BTAG_LJET_ETA_UP"); systematics.push_back("BTAG_LJET_ETA_DOWN");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream DiffSystOutputFile("MassFitDiffXSecSystematics_" + *ch + ".txt", std::ofstream::trunc);
    DiffSystOutputFile << "Top Mass, Channel: " << *ch << endl;
    DiffSystOutputFile << "Systematic & Syst. error on m_t & [GeV] \\\\" << endl;
    DiffSystOutputFile << "\\hline" << std::endl;


    extractorDiffXSec * mass_diffxs = new extractorDiffXSec(*ch,"Nominal", flags | FLAG_STORE_HISTOGRAMS);
    Double_t topmass = mass_diffxs->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_diffxs->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;

    // Systematic Variations with own samples:
    extractorDiffXSec * var, * var2;
    Double_t diff = 0;

    LOG(logDEBUG) << "Getting HAD variation...";
    var = new extractorDiffXSec(*ch,"POWHEG", flags);
    var2 = new extractorDiffXSec(*ch,"MCATNLO", flags);
    diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("HAD_UP") << " & $\\pm " << setprecision(3) << diff << "$ \\\\" << endl;
    LOG(logINFO) << "HAD - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2;

    LOG(logDEBUG) << "Getting CR variation...";
    var = new extractorDiffXSec(*ch,"PERUGIA11NoCR", flags);
    var2 = new extractorDiffXSec(*ch,"PERUGIA11", flags);
    diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("CR_UP") << " & $\\pm " << setprecision(3) << diff << "$ \\\\" << endl;
    LOG(logINFO) << "CR - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2;

    LOG(logDEBUG) << "Getting UE variation...";
    extractorDiffXSec * var3;
    var = new extractorDiffXSec(*ch,"PERUGIA11mpiHi", flags);
    var2 = new extractorDiffXSec(*ch,"PERUGIA11TeV", flags);
    var3 = new extractorDiffXSec(*ch,"PERUGIA11", flags);
    diff = (TMath::Abs(var->getTopMass() - var3->getTopMass()) + TMath::Abs(var2->getTopMass() - var3->getTopMass()))/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("UE_UP") << " & $\\pm " << setprecision(3) << diff << "$ \\\\" << endl;
    LOG(logINFO) << "UE - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2;

    // Systematic Variations produced by varying nominal samples:
    Double_t btag_syst_pos = 0, btag_syst_neg = 0;
    for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorDiffXSec * variation_diffxs = new extractorDiffXSec(*ch,*syst,flags);

      Double_t topmass_variation = variation_diffxs->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("BTAG")) {
	if(delta > 0) btag_syst_pos += delta*delta;
	else btag_syst_neg += delta*delta;
      }
      else {
	if(syst->Contains("UP")) DiffSystOutputFile << variation_diffxs->getSampleLabel((*syst)) << " & $^{" << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}_{";
	else DiffSystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      }
    }
    DiffSystOutputFile << mass_diffxs->getSampleLabel("BTAG") << " & $^{+" << setprecision(3) <<  sqrt(btag_syst_pos) << "}_{"
		       << setprecision(3) << sqrt(btag_syst_neg) << "}$ \\\\" << endl;

    total_syst_pos = sqrt(total_syst_pos);
    total_syst_neg = sqrt(total_syst_neg);
    LOG(logRESULT) << "Channel " << *ch << ": m_t = " << setprecision(6) << topmass << setprecision(3) << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV";

    DiffSystOutputFile << "Channel " << *ch << ": m_t = " << setprecision(6) << topmass << setprecision(3) << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV" << endl;
    DiffSystOutputFile.close();
    delete mass_diffxs;
  }
  return;
}

void extract() {
  Log::ReportingLevel() = Log::FromString("DEBUG3");

  bool closure = true;
  TString closure_sample = "Nominal";

  const uint32_t flags = FLAG_NORMALIZE_YIELD;

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");

  extract_diffxsec(channels,flags);
  //extract_yield(channels,closure,closure_sample,flags);

  return;
}
