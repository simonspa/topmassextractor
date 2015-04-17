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
#include "plotter.h"
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

void extractor::getTrueConfidenceIntervals(TGraphErrors * gr, Double_t cl) {
  getTrueConfidenceIntervals(gr->GetN(),1,gr->GetX(), gr->GetEY(), cl);
  for (Int_t i=0; i<gr->GetN(); i++)
    gr->SetPoint(i, gr->GetX()[i], ((TF1*)(TVirtualFitter::GetFitter())->GetUserFunc())->Eval(gr->GetX()[i]));
}

void extractor::getTrueConfidenceIntervals(Int_t n, Int_t ndim, const Double_t* x, Double_t* ci, Double_t cl) {

  // Get the confidence intervals from root:
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(n,ndim,x,ci,cl);
  // Get the fit function:
  TF1 *f = dynamic_cast<TF1*>((TVirtualFitter::GetFitter())->GetUserFunc());
  // Divide the entries by Chi2/NDF of the fit again:
  for(Int_t i = 0; i < n; i++) {
    ci[i] = ci[i]/TMath::Sqrt(f->GetChisquare()/f->GetNDF());
  }

}

TGraphErrors * extractor::createIntersectionChiSquare(TGraphErrors* data, TGraphErrors* mc_stat, Int_t bin, TGraphErrors * fit1, TGraphErrors * fit2) {

  TGraphErrors * chi2_graph = new TGraphErrors();
  TGraphErrors * chi2_graph_plotting = new TGraphErrors();
  TString gname = "chi2_" + m_channel + Form("_bin%i",bin);
  //chi2_graph->SetTitle(gname);
  chi2_graph_plotting->SetTitle(gname);

  size_t n = data->GetN();
  Double_t xmin1,xmax1,ymin1,ymax1;
  data->GetPoint(0,xmin1,ymin1);
  data->GetPoint(n-1,xmax1,ymax1);
  Double_t xmin2,xmax2,ymin2,ymax2;
  mc_stat->GetPoint(0,xmin2,ymin2);
  mc_stat->GetPoint(n-1,xmax2,ymax2);

  Double_t xmin = std::min(xmin1*0.995,xmin2*0.995);
  Double_t xmax = std::max(xmax1*1.005,xmax2*1.005);
  Double_t xmin_plotting = std::min(xmin1,xmin2);
  Double_t xmax_plotting = std::max(xmax1,xmax2);
  Double_t ymeana = (ymax1+ymin1)/2;
  Double_t ymeanb = (ymax2+ymin2)/2;

  Double_t xshift = (xmax+xmin)/2;
  Double_t yshift = (ymeana+ymeanb)/2;

  LOG(logDEBUG2) << "xmin/xmax: " << xmin << "/" << xmax << " (plotting: " << xmin_plotting << "/" << xmax_plotting << ")";

  if((flags & FLAG_DONT_SHIFT_GRAPHS) == 0) {
    // Get shifted graphs:
    LOG(logDEBUG) << "Shift all coordinates by [" << xshift << "/" << yshift << "]";
    shiftGraph(data,xshift,yshift);
    shiftGraph(mc_stat,xshift,yshift);
    xmin -= xshift; xmax -= xshift;
    xmin_plotting -= xshift; xmax_plotting -= xshift;
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
  mc_stat->Fit("pol2","F EX0 S Q","",xmin,xmax);
  TF1 * mc_statFit = mc_stat->GetFunction("pol2");
  getTrueConfidenceIntervals(scanPoints.size(),1,&scanPoints.at(0),&confIntervalMC.at(0),confidenceLevel);

  TGraphErrors *mc_statconf = new TGraphErrors(mc_stat->GetN());
  // Add additional point, just for drawing:
  mc_statconf->SetPoint(0,mc_stat->GetX()[0] - 1, 0);
  for (Int_t i = 1; i <= mc_stat->GetN(); i++) { mc_statconf->SetPoint(i, mc_stat->GetX()[i-1], 0); }
  mc_statconf->SetPoint(mc_stat->GetN(),mc_stat->GetX()[mc_stat->GetN()-1] + 1, 0);
  // Compute the confidence intervals at the x points of the created graph
  getTrueConfidenceIntervals(mc_statconf,confidenceLevel);

  data->Fit("pol2","F Q","",xmin,xmax);
  TF1 * dataFit = data->GetFunction("pol2");
  getTrueConfidenceIntervals(scanPoints.size(),1,&scanPoints.at(0),&confIntervalData.at(0),confidenceLevel);

  // Pick fixed error for data from middle point:
  Double_t fixedError = data->GetErrorY(data->GetN()/2);

  size_t i_plotting = 0;
  for(size_t i = 0; i < scanPoints.size(); i++) {
    Double_t a = dataFit->EvalPar(&scanPoints.at(i));
    Double_t b = mc_statFit->EvalPar(&scanPoints.at(i));

    // For data, do not use fitted confidence interval since all points are correlated. Use fixed error:
    Double_t awidth = fixedError;
    // For MC, use confidence interval of the fit - all points are uncorrelated:
    Double_t bwidth = confIntervalMC.at(i);

    Double_t chi2 = chiSquare(b,bwidth*bwidth,awidth*awidth,a);
    
    LOG(logDEBUG4) << "Scan " << i << "@" << scanPoints.at(i) << ": a=" << a << "(" << awidth << "/" << confIntervalData.at(i) << ") b=" << b << "(" << bwidth << ") chi2=" << chi2;
    chi2_graph->SetPoint(i, scanPoints.at(i), chi2);

    // Not all points go into the final plot, restrict to range between mass points:
    if(xmin_plotting <= scanPoints.at(i) && scanPoints.at(i) <= xmax_plotting && i%(scanPoints.size()/500) == 0) {
      chi2_graph_plotting->SetPoint(i_plotting, scanPoints.at(i), chi2);
      i_plotting++;
    }

    if(fit1) {
      fit1->SetPoint(i, scanPoints.at(i), a);
      fit1->SetPointError(i, 0, awidth);
    }
    if(fit2) {
      fit2->SetPoint(i, scanPoints.at(i), b);
      fit2->SetPointError(i, 0, bwidth);
    }
  }

  if((flags & FLAG_DONT_SHIFT_GRAPHS) == 0) { 
    // Shift all graphs back to initial position:
    shiftGraph(chi2_graph,-1*xshift,0);
    shiftGraph(chi2_graph_plotting,-1*xshift,0);
    shiftGraph(mc_statconf,-1*xshift,-1*yshift);
    shiftGraph(mc_stat,-1*xshift,-1*yshift);
    shiftGraph(data,-1*xshift,-1*yshift);
    if(fit1) shiftGraph(fit1,-1*xshift,-1*yshift);
    if(fit2) shiftGraph(fit2,-1*xshift,-1*yshift);
  }

  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TCanvas* c = 0;

    // Naming:
    TString cname = "chi2_" + m_channel + Form("_bin%i",bin);
    TString c2name = "inputs_" + m_channel + Form("_bin%i",bin);
    TString mc_statname = "input_mc_" + m_channel + Form("_bin%i",bin);
    TString dataname = "input_dat_" + m_channel + Form("_bin%i",bin);

    c = new TCanvas(cname,cname);
    c->cd();
    chi2_graph_plotting->SetMarkerStyle(20);
    chi2_graph_plotting->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2_graph_plotting->GetYaxis()->SetTitle("#chi^{2}");
    chi2_graph_plotting->Draw("AP");
    DrawDecayChLabel(m_channel,bin,bin_boundaries);
    DrawCMSLabels();
    rescaleGraph(chi2_graph_plotting,1.25);
    chi2_graph_plotting->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + cname + ".pdf"); }


    c = new TCanvas(c2name,c2name);
    c->cd();
    TLegend *leg = new TLegend();
    setLegendStyle(leg);

    mc_stat->SetTitle("");
    mc_stat->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    mc_stat->GetYaxis()->SetTitle(getQuantity());

    mc_statconf->SetTitle("");
    mc_statconf->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    mc_statconf->GetYaxis()->SetTitle(getQuantity());
    mc_statconf->GetXaxis()->SetLimits(mc_stat->GetX()[0]-1,mc_stat->GetX()[mc_stat->GetN()-1] + 1);

    setStyle(mc_stat,"madgraph");
    setStyleAndFillLegend(mc_statconf,"madgraph",leg);
    setStyleAndFillLegend(data,"data",leg, doClosure);

    mc_statconf->Draw("A E3");
    mc_stat->Draw("SAME P");
    rescaleGraph(mc_statconf);
    mc_stat->Write(mc_statname);

    data->SetPoint(0,data->GetX()[0] - 1,data->GetY()[0]);
    data->SetPoint(data->GetN()-1,data->GetX()[data->GetN()-1] + 1,data->GetY()[data->GetN()-1]);
    data->Draw("SAME L E3");
    data->Write(dataname);

    DrawDecayChLabel(m_channel,bin, bin_boundaries);
    DrawCMSLabels();
    leg->Draw();

    c->Write(c2name);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + c2name + ".pdf"); }
  }

  return chi2_graph;
}

std::pair<TGraphErrors*,TF1*> extractor::getFittedChiSquare(std::vector<Double_t> masses, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc) {

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

  // Loop over all bins we have:
  for(size_t bin = 0; bin < data.size(); bin++) {

    // Get the likelihood for the two functions:
    TGraphErrors* chi2 = createIntersectionChiSquare(data.at(bin),mc.at(bin),1+bin);

    // Discard insignificant bins if requested:
    Double_t maxVal = TMath::MaxElement(chi2->GetN(),chi2->GetY());
    if(maxVal < chi2significance) {
      LOG(logWARNING) << "Channel " << m_channel << " bin " << bin+1 << " has low significance: max(chi2) = " << maxVal << " < " << chi2significance;
    }

    // Sum them all:
    Int_t i_plotting = 0;
    for(Int_t i = 0; i < chi2->GetN(); i++) {
      Double_t xsum,ysum,x,y;
      chi2->GetPoint(i,x,y);
      chi2sum->GetPoint(i,xsum,ysum);
      LOG(logDEBUG4) << "Adding (" << x << "/" << y << ") to (" << xsum << "/" << ysum << ")";
      chi2sum->SetPoint(i,x,y+ysum);
      
      // Check if this is within the range we want to plot:
      if(masses.front() <= x && x <= masses.back() && i%(chi2->GetN()/500) == 0) {
	chi2sum_plotting->SetPoint(i_plotting,x,y+ysum);
	i_plotting++;
      }
    }
  }

  TCanvas* c = 0;
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    TString cname = "chi2_" + m_channel + "_sum";
    c = new TCanvas(cname,cname);
    c->cd();
    chi2sum_plotting->SetMarkerStyle(20);
    chi2sum_plotting->GetXaxis()->SetTitle("m_{t} [GeV]");
    chi2sum_plotting->GetYaxis()->SetTitle("#chi^{2}");
    chi2sum_plotting->Draw("AP");
    DrawDecayChLabel(m_channel);
    DrawCMSLabels();
    rescaleGraph(chi2sum_plotting);
    chi2sum_plotting->Write(gname);
    chi2sum->Write(gname);
    c->Write(cname);
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + cname + ".pdf"); }
  }

  // Fit the graph
  chi2sum->Fit("pol2","Q","",170,174.5);
  
  finalChiSquare = std::make_pair(chi2sum,chi2sum->GetFunction("pol2"));
  return finalChiSquare;
}

Double_t extractor::chiSquare(const Double_t center, const Double_t center_widthsquared, const Double_t eval_widthsquared, const Double_t eval) {
  // Do not take the "center" statictical error into account if this flag is set:
  if((flags & FLAG_IGNORE_MC_STATERR) != 0) {
    return static_cast<Double_t>(static_cast<Double_t>(eval - center)*static_cast<Double_t>(eval - center))/static_cast<Double_t>(eval_widthsquared);
  }
  // Do not take the "eval" statistical error into account if this flag is set:
  else if((flags & FLAG_INFINITE_DATA_STATISTICS) != 0) {
    return static_cast<Double_t>(static_cast<Double_t>(eval - center)*static_cast<Double_t>(eval - center))/static_cast<Double_t>(center_widthsquared);
  }
  else {
    return static_cast<Double_t>(static_cast<Double_t>(eval - center)*static_cast<Double_t>(eval - center))/static_cast<Double_t>(center_widthsquared + eval_widthsquared);
  }
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
      chi2 += chiSquare(mc.at(point)->GetBinContent(bin),0,data.at(point)->GetBinError(bin)*data.at(point)->GetBinError(bin),data.at(point)->GetBinContent(bin));
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
    DrawDecayChLabel(m_channel);
    DrawCMSLabels();
    chisquare->Write(name);
    c->Write(name + "_c");
    if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + name + ".pdf"); }
  }

  finalChiSquare = std::make_pair(chisquare,chisquare->GetFunction("pol2"));
  return finalChiSquare;
}

Double_t extractor::getMinimum(std::pair<TGraphErrors*,TF1*> finalChiSquare) {
  
  Double_t chi2min, x_chi2min = 0.0, x_left, x_right;

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

    bool firstmaximum = false, minimized = false;
    // This scanning approach requires the graph to first have a maximum and then the minimum we are looking for (excludes global minima at the edges):
    for(Int_t i = 1; i < graph->GetN(); i++) {
      if(y[i] < y[i-1] && !firstmaximum) {
	LOG(logDEBUG2) << "Found first maximum at " << x[i-1];
	firstmaximum = true;
      }
      if(y[i] > y[i-1] && firstmaximum) {
	LOG(logDEBUG2) << "Found center minimum at " << x[i-1];
	chi2min = y[i-1];
	x_chi2min = x[i-1];
	minimized = true;
	break;
      }
    }

    // We could not find a minimum in this graph!
    if(!minimized) {
      LOG(logCRITICAL) << "No suitable minimum found in this distribution!";
      // FIXME for now, just keep on running...
    }

    Double_t scanGranularity = 0.005; // GeV
    Double_t scanDistance = 5; // GeV
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

  DrawDecayChLabel(m_channel);
  DrawCMSLabels();
  c->Write(canvastitle);
  if((flags & FLAG_STORE_PDFS) != 0) { c->Print(m_outputpath + "/" + canvastitle + ".pdf"); }

  return;
}

Double_t extractor::getTopMass() {

  // Just return mass if extraction has been executed already:
  if(extractedMass > 0) { return extractedMass; }

  LOG(logDEBUG) << "Will write all output data into: " << m_outputpath;
  std::vector<TH1D*> data_hists;
  std::vector<TH1D*> mc_hists;
  std::vector<Double_t> masses;

  // Calculate the theory prediction errors to nominal sample if requested:
  if(m_requestPredictionErrors) { getPredictionUncertainties(); }

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
  if((flags & FLAG_STORE_HISTOGRAMS) != 0) {
    m_root_output = OpenFile(getRootFilename(),"update");
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

  if(m_root_output != NULL && m_root_output->IsOpen()) { m_root_output->Close(); }
  return extractedMass;
}

extractor::extractor(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags) : 
  statErrorPos(0),
  statErrorNeg(0),
  extractedMass(0),
  m_inputpath(inputpath),
  m_outputpath(outputpath),
  m_channel(ch),
  m_sample(sample),
  samples(),
  bin_boundaries(),
  m_isSystematicVariation(false),
  m_requestPredictionErrors(false),
  flags(steeringFlags),
  doClosure(false),
  m_root_output(NULL) {

  // Do not add histograms to the directory listing:
  TH1::AddDirectory(kFALSE);

  // Set the histogram styles:
  setHHStyle(*gStyle);

  // Supress Root info and warnings:
  gErrorIgnoreLevel=kError;

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
    // This is a systematic variation run:
    m_isSystematicVariation = true;

    samples.push_back(sample+"_6NEG");
    samples.push_back(sample+"_3NEG");
    samples.push_back(sample+"_1NEG");
    samples.push_back(sample);
    samples.push_back(sample+"_1POS");
    samples.push_back(sample+"_3POS");
    samples.push_back(sample+"_6POS");
  }

  // Settle the request to calculate theory prediction errors if we need it:
  if((flags & FLAG_NO_THEORYPREDICTION_ERRORS) == 0) { 
    m_requestPredictionErrors = true;
  }

  // Adapt the output path for files if we look at systematics.
  // If sample is "Nominal" but still a systematic, the child class should take
  // care of setting the output path.
  if(m_isSystematicVariation && m_sample != "Nominal") {
    m_outputpath += "/" + m_sample;
  }

  if((flags&FLAG_INFINITE_DATA_STATISTICS) != 0) { 
    LOG(logWARNING) << "This run will assume infinite statistics for the data sample!";
  }

  LOG(logDEBUG) << "Reading input files from \"" << m_inputpath << "\".";
  LOG(logDEBUG) << "Initialized. Flags shipped: " << listFlags(flags);
}
