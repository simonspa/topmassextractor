#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
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
  TH1D *binned;
  if(!binning.empty()) {
    nbins = binning.size() - 1;
    LOG(logDEBUG2) << nbins << " bins. Rebinning now.";
    binned = dynamic_cast<TH1D*>(sample->Rebin(nbins,"madgraph",&binning.front()));
    for (Int_t bin=0; bin < nbins; bin++) {
      // Divide rebinned histogram's bin content by bin width factor (new/old):
      binned->SetBinError(bin+1,sqrt(binned->GetBinContent(bin+1))/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
      binned->SetBinContent(bin+1,binned->GetBinContent(bin+1)/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
    }
  }
  else {
    LOG(logDEBUG2) << "Dividing binning by " << nbins << ". Rebinning now.";
    binned = dynamic_cast<TH1D*>(sample->Rebin(nbins,"madgraph"));
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
  int rebinning = 50;
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
    else if(!strcmp(argv[i],"-b")) { rebinning = atoi(argv[++i]); }
     // Set the verbosity level:
    else if (!strcmp(argv[i],"-v")) { Log::ReportingLevel() = Log::FromString(argv[++i]); }
    else { LOG(logERROR) << "Unrecognized command line argument \"" << argv[i] << "\".";}
  }


  massextractor::setHHStyle(*gStyle);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Loop over all requested histograms:
  for(std::vector<std::string>::iterator hgrm = histograms.begin(); hgrm != histograms.end(); hgrm++) {

    TString hist = *hgrm;

    LOG(logINFO) << "Plotting " << (*hgrm) << "...";


    // Get the reference sample:
    LOG(logDEBUG) << "Getting reference from " << files.front();
    TFile * referenceFile = new TFile(files.front().c_str());
    TH1D * reference;
    TString referencename = getSampleLabel(files.front().c_str());

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
      rebinning = bins_now/rebinning;
      LOG(logDEBUG2) << "Reference has " << bins_now << " bins, dividing by " << rebinning;

      LOG(logDEBUG2) << "Fetching " << *hgrm;
      reference = getNiceHistogram(binning, referenceFile, hist, rebinning);
      LOG(logDEBUG2) << "Fetched successfully.";

      // FIXME - no rebinning done!
    }
    // Remove data file from others:
    files.erase(files.begin());
    LOG(logDEBUG2) << "Removed reference file from list, " << files.size() << " remaining.";


    LOG(logINFO) << "Fetching histograms " << *hgrm;
    // Prepare all input histograms:
    std::vector<TH1D*> myhistograms;
    std::vector<TString> myhistogramnames;
    TH1D * tempBinned;
    for(std::vector<std::string>::iterator file = files.begin(); file != files.end(); file++) {
      LOG(logDEBUG) << "Attempting to open " << *file;
      TFile * tempFile = new TFile((*file).c_str());
      LOG(logDEBUG2) << "Opened successfully, attempting to read " << *hgrm;
      tempBinned = getNiceHistogram(binning, tempFile, hist, rebinning);
      LOG(logDEBUG2) << "Storing histogram pointer: " << tempBinned;
      LOG(logDEBUG2) << "First bin contains: " << tempBinned->GetBinContent(1);
      myhistograms.push_back(tempBinned);

      // Get the name:
      TString name = (*file);
      myhistogramnames.push_back(getSampleLabel(name));
      LOG(logDEBUG) << "Histogram read was: " << myhistogramnames.back();

      delete tempFile;
    }

    TFile output(outputpath.c_str(),"recreate");
    gDirectory->pwd();
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas(hist, "");

    TLegend *leg;
    if(hist.Contains("Rapidity")) { leg = new TLegend(0.15,0.6,0.4,.85); }
    else { leg = new TLegend(0.5,0.6,0.85,.85); }
    massextractor::setLegendStyle(leg);

    // Plot the reference:
    LOG(logDEBUG) << "Plotting reference...";
    massextractor::setStyle(reference);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    reference->SetLineColor(kGray+2);
    reference->SetLineWidth(2);
    leg->AddEntry(reference,referencename,"lp");

    // Set axis label:
    if(hist.Contains("Eta")) reference->GetXaxis()->SetTitle("#eta");
    else if(hist.Contains("pT")) reference->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if(hist.Contains("1stJetMass")) reference->GetXaxis()->SetTitle("#rho_{s}");
    else if(hist.Contains("Mass")) reference->GetXaxis()->SetTitle("m_{tt#bar} [GeV]");
    else if(hist.Contains("Rapidity")) reference->GetXaxis()->SetTitle("Rapidity");

    reference->Draw();
    reference->Write();

    std::vector<TH1*> denominators;

    // Plot all others:
    LOG(logDEBUG) << "Plotting all other histograms...";
    for(size_t i = 0; i < myhistograms.size(); i++) {
      LOG(logDEBUG2) << "Attempting to plot histogram " << i;
      TH1D* plothist = myhistograms.at(i);
      LOG(logDEBUG2) << "Histogram pointer: " << plothist;
      plothist->SetMarkerColor(i+2);
      plothist->SetLineColor(i+2);
      plothist->SetLineWidth(2);
      leg->AddEntry(plothist,myhistogramnames.at(i),"l");
      plothist->Draw("samel");
      plothist->Write();
      if(i < 7) denominators.push_back(plothist);
    }
    leg->Draw();

    while(denominators.size() < 7) denominators.push_back(NULL);

    // Draw the ratio pad:
    massextractor::drawRatio(reference, denominators.at(0), 
			     NULL, NULL, 
			     denominators.at(1), denominators.at(2), denominators.at(3), denominators.at(4), denominators.at(5), denominators.at(6), 
			     0.5,1.5);
    c->Update();
    c->Modified();
    c->Write();
  }

  return 0;
}
