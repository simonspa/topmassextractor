#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include <TROOT.h>

#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "helpers.h"
#include "plotter.h"
#include "log.h"

using namespace std;
using namespace massextractor;
using namespace unilog;

TH1D * getNiceHistogram(std::vector<Double_t> binning, TFile * file, TString hist, int nbins) {

  TH1D *sample = (TH1D*)(file->Get(hist)->Clone());
  sample->Sumw2();
  TH1D *binned;
  if(!binning.empty()) {
    nbins = binning.size() - 1;
    LOG(logDEBUG2) << nbins << " bins. Rebinning now.";
    binned = dynamic_cast<TH1D*>(sample->Rebin(nbins,"input",&binning.front()));
    for (Int_t bin=0; bin < nbins; bin++) {
      // Divide rebinned histogram's bin content by bin width factor (new/old):
      binned->SetBinError(bin+1,sqrt(binned->GetBinContent(bin+1))/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
      binned->SetBinContent(bin+1,binned->GetBinContent(bin+1)/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
      LOG(logDEBUG3) << "Bin " << bin+1 << ": " << binned->GetBinContent(bin+1) << "+-" << binned->GetBinError(bin+1);
    }
  }
  else {
    LOG(logDEBUG2) << "Dividing binning by " << nbins << ". Rebinning now.";
    binned = dynamic_cast<TH1D*>(sample->Rebin(nbins,"input"));
  }

  LOG(logDEBUG2) << "Rebinned, returning histogram ptr " << binned;
  binned->Scale(1/binned->Integral("width"));
  binned->SetDirectory(0);
  return binned;
}

int main(int argc, char* argv[]) {

  Log::ReportingLevel() = Log::FromString("INFO");

  std::vector<std::string> files;
  std::vector<std::string> histograms;
  std::string outputpath = "default.root";
  Int_t output_bins = 0;
  std::string datahistogram = "";
  bool have_data = false;
  bool rebin_data = false;
  bool mass_comparison = false;
  bool combined = false;

  for (int i = 1; i < argc; i++) {
    // Read and tokeinze the files:
    if(!strcmp(argv[i],"-f")) { files = split(string(argv[++i]), ','); }
    // Set histogram name:
    else if (!strcmp(argv[i],"-h")) { histograms = split(string(argv[++i]), ','); }
    // Set output path option
    else if(!strcmp(argv[i],"-o")) { outputpath = string(argv[++i]); }
    // Set "data" flag to rebin from data histogram:
    else if(!strcmp(argv[i],"-d")) { have_data = true; datahistogram = string(argv[++i]); }
    // Set "rebin-only" flag to rebin from data histogram without adding it:
    else if(!strcmp(argv[i],"-ro")) { rebin_data = true; datahistogram = string(argv[++i]); }
    // Set "masses" flag to draw differently:
    else if(!strcmp(argv[i],"-m")) { mass_comparison = true; }
    // Set some number of bins to rebin to:
    else if(!strcmp(argv[i],"-b")) { output_bins = atoi(argv[++i]); }
    // Combine three histograms
    else if(!strcmp(argv[i],"-c")) { combined = true; }
     // Set the verbosity level:
    else if (!strcmp(argv[i],"-v")) { Log::ReportingLevel() = Log::FromString(argv[++i]); }
    else { LOG(logERROR) << "Unrecognized command line argument \"" << argv[i] << "\".";}
  }


  massextractor::setHHStyle(*gStyle);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();

  // Open the outout file:
  TFile output(outputpath.c_str(),"recreate");
  gDirectory->pwd();

  // Loop over all requested histograms:
  for(std::vector<std::string>::iterator hgrm = histograms.begin(); hgrm != histograms.end(); hgrm++) {

    TString hist = *hgrm;
    LOG(logINFO) << "Plotting " << (*hgrm) << "...";

    // Create a subdirectory in the output file
    output.cd();
    TDirectory *folder = output.mkdir(hist);

    // Get the reference sample:
    LOG(logDEBUG) << "Getting binning from " << files.front();
    TFile * binningFile = new TFile(files.front().c_str());
    TH1D * binningHist;
    Int_t rebinning;

    std::vector<Double_t> binning;
    // If we have data, get the binning from there:
    if(have_data || rebin_data) {
      LOG(logDEBUG2) << "have data, fetching " << datahistogram;
      binningHist = (TH1D*)(binningFile->Get(datahistogram.c_str())->Clone());
      LOG(logDEBUG2) << "Fetched successfully.";

      int nbins = binningHist->GetNbinsX();
      for (Int_t bin = 1; bin <= nbins; bin++) { 
	binning.push_back(binningHist->GetBinLowEdge(bin));
      }
      binning.push_back(binningHist->GetBinLowEdge(nbins) + binningHist->GetBinWidth(nbins));
      LOG(logINFO) << "Rebinning with " << binning.size()-1 << " bins.";

    }
    // Otherwise rebin to CLI option or to default:
    else {
      Int_t bins_now = ((TH1D*)binningFile->Get(hist))->GetNbinsX();
      if(output_bins == 0) output_bins = bins_now;
      LOG(logDEBUG2) << "Binning reference has " << bins_now << " bins, requested " << output_bins << ", dividing by " << (bins_now/output_bins);
      rebinning = bins_now/output_bins;
    }

    // Getting the reference histogram:
    bool have_reference = false;
    TString referencename;
    if(have_data) referencename = "Data";
    else referencename = "MadGraph";
    //else referencename = getSampleLabel(files.front().c_str());

    TH1D * reference;
    if (have_data) {
      TFile * referenceFile = new TFile(files.front().c_str());
      LOG(logDEBUG2) << "Have data, fetching " << datahistogram;
      reference = (TH1D*)(referenceFile->Get(datahistogram.c_str())->Clone());
      LOG(logDEBUG2) << "Fetched successfully.";
      reference->Scale(1/reference->Integral("width"));

      have_reference = true;
      files.erase(files.begin());
    }
    else if(rebin_data) {
      // First histogram is data, but we don't use it:
      files.erase(files.begin());
    }
    LOG(logDEBUG2) << "Skip reference file, " << files.size() << " files remaining.";

    LOG(logINFO) << "Fetching histograms " << *hgrm;
    // Prepare all input histograms:
    std::vector<TH1D*> myhistograms;
    std::vector<TString> myhistogramnames;
    TH1D * tempBinned, * tempBinned2, * tempBinned3;
    for(std::vector<std::string>::iterator file = files.begin(); file != files.end(); file++) {
      LOG(logDEBUG) << "Attempting to open " << *file;
      TFile * tempFile = new TFile((*file).c_str());
      LOG(logDEBUG2) << "Opened successfully, attempting to read " << *hgrm;
      tempBinned = getNiceHistogram(binning, tempFile, hist, rebinning);
      LOG(logDEBUG2) << "Storing histogram pointer: " << tempBinned;
      LOG(logDEBUG2) << "First bin contains: " << tempBinned->GetBinContent(1);

      // Sum three for combined channel:
      if(combined) {
	delete tempFile;
	file++;
	LOG(logDEBUG) << "Getting the others: " << *file;
	tempFile = new TFile((*file).c_str());
	tempBinned2 = getNiceHistogram(binning, tempFile, hist, rebinning);
	delete tempFile;
	file++;
	LOG(logDEBUG) << "Getting the others: " << *file;
	tempFile = new TFile((*file).c_str());
	tempBinned3 = getNiceHistogram(binning, tempFile, hist, rebinning);

	LOG(logDEBUG) << "Summing channels...";
	tempBinned->Add(tempBinned2);
	tempBinned->Add(tempBinned3);
	/*Int_t nbins = tempBinned->GetNbinsX();
	for(Int_t bin = 1; bin <= nbins; bin++) {
	  tempBinned->SetBinContent(bin,((tempBinned->GetBinContent(bin)/(tempBinned->GetBinError(bin)*tempBinned->GetBinError(bin))
					  + tempBinned2->GetBinContent(bin)/(tempBinned2->GetBinError(bin)*tempBinned2->GetBinError(bin))
					  + tempBinned3->GetBinContent(bin)/(tempBinned3->GetBinError(bin)*tempBinned3->GetBinError(bin)))/
					 ((1/(tempBinned->GetBinError(bin)*tempBinned->GetBinError(bin)))
					  + (1/(tempBinned2->GetBinError(bin)*tempBinned2->GetBinError(bin)))
					  + (1/(tempBinned3->GetBinError(bin)*tempBinned3->GetBinError(bin))))));
	  tempBinned->SetBinError(bin,(1/(TMath::Sqrt((1/(tempBinned->GetBinError(bin)*tempBinned->GetBinError(bin)))
						      + (1/(tempBinned2->GetBinError(bin)*tempBinned2->GetBinError(bin)))
						      + (1/(tempBinned3->GetBinError(bin)*tempBinned3->GetBinError(bin)))))));
						      }*/
	LOG(logDEBUG2) << "Summed, first bin is: " << tempBinned->GetBinContent(1);
	tempBinned->Scale(1/tempBinned->Integral("width"));
	LOG(logDEBUG2) << "First bin now contains: " << tempBinned->GetBinContent(1);
       }

      // Get the name:
      TString name = (*file);
      TString title;
      if(!mass_comparison) {
	if(name.Contains("POWHEG_13")) { title = "POWHEG ttJ"; }
	else if(name.Contains("alpha")) { title = "POWHEG ttJ, #alpha = 0.1273"; }
	else if(name.Contains("POWHEG_15")) { title = "POWHEG ttJ (15)"; }
	else if(name.Contains("POWHEG_16")) { title = "POWHEG ttJ (16)"; }
	else if(name.Contains("POWHEG_17")) { title = "POWHEG ttJ, #alpha = 0.1365"; }
	else if(name.Contains("SCALEUP")) { title = "ttJ PS Scale Up"; }
	else if(name.Contains("SCALEDOWN")) { title = "ttJ PS Scale Down"; }
	else if(name.Contains("SCALE_UP_nopt")) { title = "ttJ ME Scale Up (ptHard=0)"; }
	else if(name.Contains("SCALE_DOWN_nopt")) { title = "ttJ ME Scale Down (ptHard=0)"; }
	else if(name.Contains("SCALE_UP")) { title = "ttJ ME Scale Up"; }
	else if(name.Contains("SCALE_DOWN")) { title = "ttJ ME Scale Down"; }
	else if(name.Contains("nopthard")) { title = "POWHEG ttJ ptHard=0"; }
	else if(name.Contains("powhegbox")) { title = "POWHEG ttJ"; }
	else if(name.Contains("fxfx")) { title = "aMCatNLO fxfx"; }
	else if(name.Contains("powheg")) { title = "POWHEG"; }
	else if(name.Contains("Nominal")) { title = "MadGraph"; }
	else { title = getSampleLabel(name); }
      }
      else {
	//if(name.Contains("MASS")) {
	  Double_t mass = 0;
	  if(name.Contains("178")) mass = 178.5;
	  else if(name.Contains("175")) mass = 175.5;
	  else if(name.Contains("169")) mass = 169.5;
	  else if(name.Contains("166")) mass = 166.5;
	  else if(name.Contains("massup")) mass = 173.5;
	  else if(name.Contains("massdown")) mass = 172.5;
	  else mass = 172.5;
	  
	  //if(name.Contains("fxfx")) { title = Form("aMCatNLO fxfx, %3.1f GeV",mass); }
	  if(name.Contains("fxfx")) { title = "aMCatNLO fxfx"; }
	  else title = Form("m_{t}^{MC} = %3.1f GeV",mass); 
	  //}
	  //else { myhistogramnames.push_back(getSampleLabel(name)); }
      }
      LOG(logDEBUG) << "Histogram read was: " << title;

      if(!have_reference) {
	reference = tempBinned;
	referencename = title;
	have_reference = true;
      }
      else {
	myhistograms.push_back(tempBinned);
	myhistogramnames.push_back(title);
      }
      
      delete tempFile;
    }

    // Descend into the folder:
    folder->cd();
    TCanvas *c = new TCanvas(hist,"");
    gStyle->SetOptStat(0);

    TLegend *leg;
    if(hist.Contains("Rapidity")) { leg = new TLegend(0.15,0.6,0.4,.85); }
    else { leg = new TLegend(0.5,0.6,0.85,.85); }
    massextractor::setLegendStyle(leg);

    // Plot the reference:
    LOG(logDEBUG) << "Plotting reference...";
    massextractor::setStyle(reference);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    if(!have_data) reference->SetMarkerStyle(0);
    setStyle(reference,"reference");
    //if(mass_comparison) { reference->SetMarkerStyle(0); }
    reference->SetLineColor(kBlack);
    reference->SetLineWidth(2);
    reference->SetTitle("");
    //if(!mass_comparison) { 
    leg->AddEntry(reference,referencename,"lp"); 
    //}

    // Set axis label:
    std::string quantity = "";
    std::string unit = "";
    if(hist.Contains("Eta")) {             quantity = "#eta"; }
    else if(hist.Contains("pT")) {         quantity = "p_{T}"; unit = "GeV"; }
    else if(hist.Contains("1stJetMass")) { quantity = "#rho_{s}"; }
    else if(hist.Contains("Mass")) {       quantity = "m_{tt#bar}"; unit = "GeV"; }
    else if(hist.Contains("Rapidity")) {   quantity = "y"; }
    
    if(unit != "") {
      reference->GetXaxis()->SetTitle(Form("%s [%s]",quantity.c_str(),unit.c_str()));
      reference->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s} [%s^{-1}]",quantity.c_str(),unit.c_str()));
    }
    else {
      reference->GetXaxis()->SetTitle(Form("%s",quantity.c_str()));
      reference->GetYaxis()->SetTitle(Form("#frac{1}{#sigma} #frac{d#sigma}{d%s}",quantity.c_str()));
    }
    reference->GetYaxis()->SetTitleOffset(1.8);
    reference->GetYaxis()->SetLabelOffset(0.01);

    reference->Draw("e");
    if(have_data) DrawCMSLabels();
    else DrawFreeCMSLabels("Simulation");
    reference->Write(referencename);


    TString name = files.front().c_str();
    if(combined || name.Contains("combined")) DrawDecayChLabel("combined");
    else if(name.Contains("ee")) DrawDecayChLabel("ee");
    else if(name.Contains("emu")) DrawDecayChLabel("emu");
    else if(name.Contains("mumu")) DrawDecayChLabel("mumu");

    std::vector<TH1*> denominators;

    // Plot all others:
    LOG(logDEBUG) << "Plotting all other histograms...";
    for(size_t i = 0; i < myhistograms.size(); i++) {
      LOG(logDEBUG2) << "Attempting to plot histogram " << i;
      TH1D* plothist = myhistograms.at(i);
      LOG(logDEBUG2) << "Histogram pointer: " << plothist;
      if(!mass_comparison) {
	if(i == 0) {
	  plothist->SetLineColor(kGray+2);
	  plothist->SetLineStyle(7);
	}
	else if(i == 1) {	plothist->SetLineColor(4); }
	else if(i == 2) {	plothist->SetLineColor(kGreen+1); }
	else { plothist->SetLineColor(i+2); }
      }
      else {
	if(i == 2) {
	  plothist->SetLineColor(kGray+2);
	}

	if(i == 1 || i == 3) { plothist->SetLineStyle(7); }

	if(i == 0) { plothist->SetLineColor(kRed+1); }
	else if(i == 1) { plothist->SetLineColor(kOrange+8); }
	//else if(i == 2) { plothist->SetLineColor(kOrange-4); }
	//else if(i == 3) { plothist->SetLineColor(kAzure-3); }
	else if(i == 3) { plothist->SetLineColor(kBlue-4); }
	else if(i == 4) { plothist->SetLineColor(kBlue+2); }
	else if(i == 5) { plothist->SetLineColor(kGreen+1); }
	//else { plothist->SetLineColor(i+2); }
      }

      plothist->SetMarkerStyle(0);
      plothist->SetLineWidth(2);
      leg->AddEntry(plothist,myhistogramnames.at(i),"l");
      plothist->Draw("same e");
      plothist->Write(myhistogramnames.at(i));
      if(i < 7) denominators.push_back(plothist);

      //if(mass_comparison && i == (myhistograms.size()-1)/2) { leg->AddEntry(reference,referencename,"l"); }
    }
    leg->Draw();

    while(denominators.size() < 7) denominators.push_back(NULL);

    // Prepare the uncertainty band:
    TH1* uncBand = NULL;
    uncBand = dynamic_cast<TH1*>(reference->Clone("uncBand"));
    uncBand->SetFillStyle(3354);
    uncBand->SetFillColor(kBlack);
    uncBand->SetLineColor(kBlack);
    gStyle->SetHatchesLineWidth(1);
    gStyle->SetHatchesSpacing(0.8);
    uncBand->SetMarkerStyle(0);

    TGraphAsymmErrors * ratio_stat;
    TGraphAsymmErrors * ratio_syst;
    
    const Int_t n = 5;
    Double_t x[n]   = {0.1,0.25,0.38,0.53,0.8};
    Double_t y[n]   = {1.0,1.0,1.0,1.0,1.0};
    Double_t exl[n] = {0.1,.05,.08,.08,.2};
    Double_t exh[n] = {0.1,.05,.08,.08,.2};

    // STAT:
    Double_t eyl[n] = {0.13,0.052,0.03,0.025,0.051};
    // TOT
    Double_t eyh[n] = {0.199,0.117,0.061,0.048,0.097};

    ratio_stat = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyl);
    //ratio_stat->SetFillStyle(3005);
    ratio_stat->SetFillColor(kOrange-4);
    ratio_stat->SetLineColor(0);

    ratio_syst = new TGraphAsymmErrors(n,x,y,exl,exh,eyh,eyh);
    //ratio_syst->SetFillStyle(3005);
    ratio_syst->SetFillColor(kGray+2);
    ratio_syst->SetLineColor(0);

    // Draw the ratio pad:
    massextractor::drawRatio(reference, denominators.at(0), 
			     uncBand,
			     //NULL, 
			     //ratio_stat, ratio_syst, 
			     NULL, NULL,
			     denominators.at(1), denominators.at(2), denominators.at(3), denominators.at(4), denominators.at(5), denominators.at(6), 
			     0.7,1.5,have_data);


    c->Update();
    c->Modified();
    c->Write();
  }

  return 0;
}
