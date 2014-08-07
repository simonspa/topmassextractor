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
#include "log.h"
#include "extract.h"

Double_t nominalmass = 172.5;
Double_t lumi = 19712;

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

Double_t extractor::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Calculate the signal fraction from reconstructed events and TT background:
  Double_t fsignal = reco/(reco+ttbgr);

  // Calculate the "others backgorund" as difference between bgr and ttbgr (no xsec scaling)
  Double_t bgr_other = bgr - ttbgr;

  // Calculate signal by subtracting backround from data, multiplied by signal fraction.
  Double_t signal = (data - bgr_other)*fsignal;

  LOG(logDEBUG2) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

  return signal;
}

Double_t extractorBackground::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

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
  if(!doClosure) aDataHist = static_cast<TH1D*>(histos->Get("aDataHist")->Clone());
  else aDataHist = static_cast<TH1D*>(pseudoData->Clone());

  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist"));
  // Histograms containing the background:
  TH1D * aTtBgrHist = static_cast<TH1D*>(histos->Get("aTtBgrHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(histos->Get("aBgrHist"));

  Int_t nbins = aDataHist->GetNbinsX();
  LOG(logDEBUG) << "Data hist has " << nbins << " bins.";

  // Iterate over all bins:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    // Get signal corrected by ttbar background:
    Double_t signal = getSignal(bin, mass, aDataHist->GetBinContent(bin),
				aRecHist->GetBinContent(bin),
				aBgrHist->GetBinContent(bin),
				aTtBgrHist->GetBinContent(bin));
    
    // Write background subtrated signal:
    aDataHist->SetBinContent(bin,signal);
  }

  // Return signal-only histogram:
  return aDataHist;
}

TH1D * extractor::getSimulationHistogram(Double_t mass, TFile * histos) {

  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist")->Clone());

  Int_t nbins = aRecHist->GetNbinsX();
  LOG(logDEBUG) << "Reco hist has " << nbins << " bins.";

  // Iterate over all bins:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    // Correct the Reco events for different TTBar Cross sections (mass dependent):
    Double_t corr_reco = getReco(bin, mass,aRecHist->GetBinContent(bin));

    // Write corrected Reco:
    aRecHist->SetBinContent(bin,corr_reco);
  }

  // Return reco histogram:
  return aRecHist;
}

std::vector<std::vector<Double_t> > extractor::splitBins(std::vector<TH1D*> histograms) {

  std::vector<std::vector<Double_t> > separated_bins;
  Int_t nbins = histograms.at(0)->GetNbinsX();

  // For every bin, prepare a new histogram containing all mass points:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    LOG(logDEBUG) << "Filling bin " << bin << " mass points ";
 
    // Prepare new vector for this bin:
    std::vector<Double_t> thisbin;

    for(std::vector<TH1D*>::iterator hist = histograms.begin(); hist != histograms.end(); ++hist) {
      // Set the corresponding bin in output vector:
      thisbin.push_back((*hist)->GetBinContent(bin));
    }

    separated_bins.push_back(thisbin);
  }
  
  return separated_bins;
}

std::vector<TF1*> extractor::fitMassBins(TString channel, Int_t bin, std::vector<Double_t> masses, std::vector<Double_t> data, std::vector<Double_t> mc) {

  std::vector<TF1*> allfits;

  TString mname, dname, cname;
  mname.Form("mc_%i_",bin);
  dname.Form("dat_%i_",bin);
  cname.Form("dat_mc_%i_",bin);

  TCanvas* c;
  if(storeHistograms) {
    c = new TCanvas(cname+channel,cname+channel);
    c->cd();
  }

  TGraphErrors * graph_mc = new TGraphErrors();
  graph_mc->SetTitle(mname+channel);

  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph_mc->SetPoint(point, masses.at(point), mc.at(point));
    //graph_mc->SetPointError(point,0,sqrt(mc.at(point)));
  }

  //graph_mc->Fit("pol2","Q");

  TF1 * mcfit = graph_mc->GetFunction("pol2");
  allfits.push_back(mcfit);

  TLegend *leg = new TLegend();
  setLegendStyle(leg);

  if(storeHistograms) {
    graph_mc->SetTitle("");
    graph_mc->GetXaxis()->SetTitle("m_{t} #left[GeV#right]");
    graph_mc->GetYaxis()->SetTitle("Events");

    setStyleAndFillLegend(graph_mc,"madgraph",leg);
    graph_mc->Draw("A P E1");
    graph_mc->Write(mname+channel);
  }

  TGraphErrors * graph = new TGraphErrors();
  graph->SetTitle(dname+channel);

  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph->SetPoint(point, masses.at(point), data.at(point));
    graph->SetPointError(point,0,sqrt(data.at(point)));
  }

  //graph->Fit("pol2","Q");

  TF1 * datfit = graph->GetFunction("pol2");
  allfits.push_back(datfit);

  if(storeHistograms) {
    setStyleAndFillLegend(graph,"data",leg);

    graph->Draw("SAME P E1");
    DrawDecayChLabel(getChannelLabel(channel));
    DrawCMSLabels();

    // Also draw legend:
    leg->Draw();

    graph->Write(dname+channel);
    c->Write();
  }

  return allfits;
}

TF1 * extractor::getChiSquare(TString channel, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc) {

  TString name = "chi2_";

  TCanvas* c;
  if(storeHistograms) {
    c = new TCanvas(name+channel+"_c",name+channel+"_c");
    c->cd();
  }

  TGraphErrors * chisquare = new TGraphErrors();
  chisquare->SetTitle(name+channel);

  for(UInt_t point = 0; point < masses.size(); ++point) {

    Double_t chi2 = 0;

    // Iterate over all bins:
    for(Int_t bin = 1; bin <= data.at(point)->GetNbinsX(); bin++) {
      Double_t chi = (data.at(point)->GetBinContent(bin) - mc.at(point)->GetBinContent(bin))/sqrt(data.at(point)->GetBinContent(bin));
      chi2 += chi*chi;
    }
    chisquare->SetPoint(point, masses.at(point), chi2);
  }

  chisquare->Fit("pol2","Q");
  if(storeHistograms) {
    chisquare->SetMarkerStyle(20);
    chisquare->GetXaxis()->SetTitle("m_{t} [GeV]");
    chisquare->GetYaxis()->SetTitle("#chi^{2}");
    chisquare->Draw("AP");
    DrawDecayChLabel(getChannelLabel(channel));
    DrawCMSLabels();
    chisquare->Write(name+channel);
    c->Write();
  }

  return chisquare->GetFunction("pol2");
}

Double_t extractor::getMinimum(TF1 * fit) {
  
  // For now, just return the function's minimum:
  Double_t chi2min = fit->GetMinimum(0,330);
  Double_t x_chi2min = fit->GetMinimumX(0,330);

  // Statictical error: vary chi2 by +-1 (going up the curve left and right by dChi2 = 1):
  Double_t x_left = x_chi2min - fit->GetX(chi2min+1,0,x_chi2min);
  Double_t x_right = fit->GetX(chi2min+1, x_chi2min,330) - x_chi2min;

  LOG(logDEBUG) << "Minimum Chi2 is " << chi2min << " at " << x_chi2min << " +" << x_right << "-" << x_left;
  
  // Store the averaged statistical error:
  statError = (x_right+x_left)/2;
  return x_chi2min;
}

Double_t extractor::getStatError() {

  // Just return the statistical error calculated from chi2:
  return statError;
}

Double_t extractor::getTopMass() {

  std::vector<TH1D*> data_hists;
  std::vector<TH1D*> mc_hists;
  std::vector<Double_t> masses;

  for(std::vector<TString>::iterator sample = samples.begin(); sample != samples.end(); ++sample) {
    // Input files:
    TString filename = "preunfolded/" + (*sample) + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    TFile * datafile = new TFile(filename);

    if(!datafile->IsOpen()) {
      LOG(logINFO) << "Failed to access data file " << filename;
      throw 1;
    }

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

    datafile->Close();
  }

  std::vector<std::vector<Double_t> > separated_data = splitBins(data_hists);
  std::vector<std::vector<Double_t> > separated_mc = splitBins(mc_hists);

  TFile output;
  if(storeHistograms) {
    output.Open("MassFitRates.root","update");
    gDirectory->pwd();
  }

  std::vector<std::vector<TF1*> > fits;
  for(UInt_t bin = 0; bin < separated_data.size(); ++bin) {
    fits.push_back(fitMassBins(channel,bin+1,masses,separated_data.at(bin),separated_mc.at(bin)));
  }

  TF1 * fit = getChiSquare(channel,masses,data_hists,mc_hists);
  extractedMass = getMinimum(fit);

  //getChiSquareFitted(channel,masses,data_hists,mc_hists);

  output.Close();
  return extractedMass;
}

void extractorOtherSamples::calcDifferenceToNominal(TString nominal, TString systematic) {

  if(systematic.Contains("HAD")) {
    if(systematic.Contains("UP")) { systematic = "MCATNLO"; }
    else if(systematic.Contains("DOWN")) { systematic = "POWHEG"; }
  }

  // Input files:
  TString nfilename = "preunfolded/" + nominal + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TString sfilename = "preunfolded/" + systematic + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";

  TFile * nominalfile = new TFile(nfilename);
  TFile * systematicfile = new TFile(sfilename);

  if(!nominalfile->IsOpen() || !systematicfile->IsOpen()) {
    LOG(logINFO) << "Failed to access file";
    throw 1;
  }

  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  TH1D * nominalReco = static_cast<TH1D*>(nominalfile->Get("aRecHist"));
  TH1D * varReco = static_cast<TH1D*>(systematicfile->Get("aRecHist"));

  TH1D * nominalBgr = static_cast<TH1D*>(nominalfile->Get("aBgrHist"));
  TH1D * varBgr = static_cast<TH1D*>(systematicfile->Get("aBgrHist"));

  TH1D * nominalTtbgr = static_cast<TH1D*>(nominalfile->Get("aTtBgrHist"));
  TH1D * varTtbgr = static_cast<TH1D*>(systematicfile->Get("aTtBgrHist"));

  for(Int_t bin = 1; bin <= nominalReco->GetNbinsX(); bin++) {
    Double_t rec = nominalReco->GetBinContent(bin) - varReco->GetBinContent(bin);
    deltaRec.push_back(rec);
    
    Double_t bgr = nominalBgr->GetBinContent(bin) - varBgr->GetBinContent(bin);
    deltaBgr.push_back(bgr);

    Double_t ttbgr = nominalTtbgr->GetBinContent(bin) - varTtbgr->GetBinContent(bin);
    deltaTtbgr.push_back(ttbgr);

    LOG(logDEBUG2) << "Diff bin #" << bin << " reco: " << nominalReco->GetBinContent(bin) << " - " << varReco->GetBinContent(bin) << " = " << rec << " dbgr=" << bgr << " dttbgr=" << ttbgr;
  }

  nominalfile->Close();
  systematicfile->Close();
}

Double_t extractor::getMassFromSample(TString sample) {

  Double_t topmass = nominalmass;
  // The mass samples:
  if(sample.Contains("MASS")) {
    if(sample.Contains("UP")) {
      if(sample.Contains("1GEV")) topmass += 1;
      else if(sample.Contains("3GEV")) topmass += 3;
      else if(sample.Contains("6GEV")) topmass += 6;
    }
    else if(sample.Contains("DOWN")) {
      if(sample.Contains("1")) topmass -= 1;
      else if(sample.Contains("3GEV")) topmass -= 3;
      else if(sample.Contains("6GEV")) topmass -= 6;
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

void extractor::setClosureSample(TString closure) {

  // Enable closure test:
  doClosure = true;

  // Input file:
  TString filename = "preunfolded/" + closure + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TFile * closureFile = new TFile(filename);

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

  closureFile->Close();
  delete closureFile;
}

extractor::extractor(TString ch, TString sample, bool storeHistos) : channel(ch), samples(), extractedMass(0), statError(0), storeHistograms(storeHistos), doClosure(false) {

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
  
  LOG(logDEBUG) << "Initialized.";
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

TString extractor::getChannelLabel(TString channel) {
 
  TString label = "";
  if(channel =="ee") { label = "ee"; }
  if(channel =="mumu"){ label = "#mu#mu"; }
  if(channel =="emu"){ label = "e#mu"; }
  if(channel =="combined"){ label = "Dilepton Combined"; }

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

void extract() {

  Log::ReportingLevel() = Log::FromString("INFO");

  bool closure = true;
  TString closure_sample = "Nominal";

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");

  // Do the same (or similar things) for the systematic uncertainties:
  std::vector<TString> syst_on_nominal;
  syst_on_nominal.push_back("MATCH_UP");
  syst_on_nominal.push_back("MATCH_DOWN");
  syst_on_nominal.push_back("SCALE_UP");
  syst_on_nominal.push_back("SCALE_DOWN");
  syst_on_nominal.push_back("HAD_UP");
  syst_on_nominal.push_back("HAD_DOWN");

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


    extractor * mass_samples = new extractor(*ch,"Nominal",true);
    if(closure) mass_samples->setClosureSample(closure_sample);

    Double_t topmass = mass_samples->getTopMass();
    Double_t total_stat = mass_samples->getStatError();
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +- " << total_stat;
    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;

    // Systematic Variations with own samples:
    for(std::vector<TString>::iterator syst = syst_on_nominal.begin(); syst != syst_on_nominal.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorOtherSamples * matchscale_samples = new extractorOtherSamples(*ch,"Nominal",false,(*syst));
      if(closure) matchscale_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = matchscale_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("UP")) SystOutputFile << (*syst) << " & $^{" << setprecision(3) << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
    }


    // Getting DrellYan/Background variations:
    for(std::vector<TString>::iterator syst = syst_bg.begin(); syst != syst_bg.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorBackground * bg_samples = new extractorBackground(*ch,"Nominal",false,(*syst));
      if(closure) bg_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = bg_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;
      
      if(syst->Contains("UP")) SystOutputFile << (*syst) << " & $^{" << setprecision(3) << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
    }

    // Systematic Variations produced by varying nominal samples:
    for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractor * variation_samples = new extractor(*ch,*syst,false);
      if(closure) variation_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = variation_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("UP")) SystOutputFile << (*syst) << " & $^{+" << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(3) <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
    }

    total_syst_pos = sqrt(total_syst_pos);
    total_syst_neg = sqrt(total_syst_neg);
    LOG(logRESULT) << "Channel " << *ch << ": m_t = " << topmass << " +-" << total_stat << " +" << total_syst_pos << " -" << total_syst_neg << " GeV";

    SystOutputFile << "Channel " << *ch << ": m_t = " << topmass << " +-" << total_stat << " +" << total_syst_pos << " -" << total_syst_neg << " GeV";
    SystOutputFile.close();
  }

  return;
}
