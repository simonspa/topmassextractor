#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
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
  Int_t output_bins;
  std::string datahistogram = "";
  bool have_data = false;

  for (int i = 1; i < argc; i++) {
    // Read and tokeinze the files:
    if(!strcmp(argv[i],"-f")) { files = split(string(argv[++i]), ','); }
    // Set histogram name:
    else if (!strcmp(argv[i],"-h")) { histograms = split(string(argv[++i]), ','); }
    // Set output path option
    else if(!strcmp(argv[i],"-o")) { outputpath = string(argv[++i]); }
    // Set "data" flag to rebin from data histogram:
    else if(!strcmp(argv[i],"-d")) { have_data = true; datahistogram = string(argv[++i]); }
   // Set some number of bins to rebin to:
    else if(!strcmp(argv[i],"-b")) { output_bins = atoi(argv[++i]); }
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
    LOG(logDEBUG) << "Getting reference from " << files.front();
    TFile * referenceFile = new TFile(files.front().c_str());
    TH1D * reference;
    Int_t rebinning;
    TString referencename;
    if(have_data) referencename = "CMS Data";
    else referencename = "MadGraph";
    //else referencename = getSampleLabel(files.front().c_str());

    std::vector<Double_t> binning;
    // If we have data, get the binning from there:
    if(have_data) {
      LOG(logDEBUG2) << "have data, fetching " << datahistogram;
      reference = (TH1D*)(referenceFile->Get(datahistogram.c_str())->Clone());
      LOG(logDEBUG2) << "Fetched successfully.";
      reference->Scale(1/reference->Integral("width"));

      int nbins = reference->GetNbinsX();
      for (Int_t bin = 1; bin <= nbins; bin++) { 
	binning.push_back(reference->GetBinLowEdge(bin));
      }
      binning.push_back(reference->GetBinLowEdge(nbins) + reference->GetBinWidth(nbins));
      LOG(logINFO) << "Rebinning with " << binning.size()-1 << " bins.";

    }
    // Otherwise rebin to CLI option of to default:
    else {
      Int_t bins_now = ((TH1D*)referenceFile->Get(hist))->GetNbinsX();
      rebinning = bins_now/output_bins;
      LOG(logDEBUG2) << "Reference has " << bins_now << " bins, dividing by " << rebinning;
      LOG(logDEBUG2) << "Fetching " << *hgrm;
      reference = getNiceHistogram(binning, referenceFile, hist, rebinning);
      LOG(logDEBUG2) << "Fetched successfully.";
    }
    LOG(logDEBUG2) << "Skip reference file, " << files.size() << " files remaining.";


    LOG(logINFO) << "Fetching histograms " << *hgrm;
    // Prepare all input histograms:
    std::vector<TH1D*> myhistograms;
    std::vector<TString> myhistogramnames;
    TH1D * tempBinned;
    for(std::vector<std::string>::iterator file = files.begin() + 1; file != files.end(); file++) {
      LOG(logDEBUG) << "Attempting to open " << *file;
      TFile * tempFile = new TFile((*file).c_str());
      LOG(logDEBUG2) << "Opened successfully, attempting to read " << *hgrm;
      tempBinned = getNiceHistogram(binning, tempFile, hist, rebinning);
      LOG(logDEBUG2) << "Storing histogram pointer: " << tempBinned;
      LOG(logDEBUG2) << "First bin contains: " << tempBinned->GetBinContent(1);
      myhistograms.push_back(tempBinned);

      // Get the name:
      TString name = (*file);
      if(name.Contains("POWHEG_15")) { myhistogramnames.push_back("POWHEG ttJ w/o hdamp"); }
      else if(name.Contains("powhegbox")) { myhistogramnames.push_back("POWHEG ttJ"); }
      else if(name.Contains("powheg")) { myhistogramnames.push_back("POWHEG"); }
      else if(name.Contains("Nominal")) { myhistogramnames.push_back("MadGraph"); }
      else { myhistogramnames.push_back(getSampleLabel(name)); }
      LOG(logDEBUG) << "Histogram read was: " << myhistogramnames.back();

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
    reference->SetLineColor(kBlack);
    reference->SetLineWidth(2);
    reference->SetTitle("");
    leg->AddEntry(reference,referencename,"lp");

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

    std::vector<TH1*> denominators;

    // Plot all others:
    LOG(logDEBUG) << "Plotting all other histograms...";
    for(size_t i = 0; i < myhistograms.size(); i++) {
      LOG(logDEBUG2) << "Attempting to plot histogram " << i;
      TH1D* plothist = myhistograms.at(i);
      LOG(logDEBUG2) << "Histogram pointer: " << plothist;
      if(i == 0) {
	plothist->SetLineColor(kGray+2);
	plothist->SetLineStyle(7);
      }
      else if(i == 1) {	plothist->SetLineColor(4); }
      else if(i == 2) {	plothist->SetLineColor(kRed+1); }
      else { plothist->SetLineColor(i+2); }
      plothist->SetMarkerStyle(0);
      plothist->SetLineWidth(2);
      leg->AddEntry(plothist,myhistogramnames.at(i),"l");
      plothist->Draw("same e");
      plothist->Write(myhistogramnames.at(i));
      if(i < 7) denominators.push_back(plothist);
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

    // Draw the ratio pad:
    massextractor::drawRatio(reference, denominators.at(0), uncBand,
			     NULL, NULL, 
			     denominators.at(1), denominators.at(2), denominators.at(3), denominators.at(4), denominators.at(5), denominators.at(6), 
			     0.5,1.5,have_data);
    c->Update();
    c->Modified();
    c->Write();
  }

  return 0;
}
