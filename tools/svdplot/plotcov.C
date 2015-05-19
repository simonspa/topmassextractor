#include <iostream>
#include <vector>
#include <string>
#include <sstream>
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

Double_t nominalmass = 172.5;
Double_t lumi = 19712;

TH2D * readCov(TString ch) {

  TString filename = "SVD/Nominal/Unfolding_" + ch + "_TtBar_Mass_HypTTBar1stJetMass.root";
  TString histogramname = "SVD_" + ch + "_TtBar_Mass_HypTTBar1stJetMass_Nominal_STATCOV";
  
  TFile * file = TFile::Open(filename,"read");
  if(!file->IsOpen()) {
    std::cout << "Failed to open file " << filename << std::endl;
    throw 1;
  }
  else std::cout << "Opened file " << filename << std::endl;
  
  // Histogram containing differential cross section from data:
  TH2D * statCovNorm = static_cast<TH2D*>(file->Get(histogramname));
  if(!statCovNorm) { std::cout << "Failed to get histogram " << histogramname << std::endl; }
  else { std::cout << "Fetched " << histogramname << std::endl; }
  
  std::vector<Double_t> binning;
  for (Int_t bin = 2; bin < statCovNorm->GetNbinsX(); bin++) { 
    binning.push_back(statCovNorm->GetBinLowEdge(bin));
      std::cout << " " << binning.back();
  }
  binning.push_back(statCovNorm->GetBinLowEdge(statCovNorm->GetNbinsX()-1) + statCovNorm->GetBinWidth(statCovNorm->GetNbinsX()-1));
  std::cout << " " << binning.back() << std::endl;
  
  // Cut away over and underflow bins:
  TH2D * cov = new TH2D("cov","cov",statCovNorm->GetNbinsX()-2,&binning.front(),statCovNorm->GetNbinsY()-2,&binning.front());
  for(Int_t y = 0; y <= cov->GetNbinsX(); y++) {
    for(Int_t x = 0; x <= cov->GetNbinsY(); x++) {
      cov->SetBinContent(x,y,statCovNorm->GetBinContent(x+1,y+1));
    }
  }

  return cov;
}

void plotcov() {

  // Compile the plotter utils:
  gROOT->ProcessLine(".L ../../utils/plotter.cc+");

  bool drawSys = false;

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");
  
  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    massextractor::setHHStyle(*gStyle);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.16);

    TH2D *cov;
    if(*ch != "combined") cov = readCov(*ch);
    // Combined channel - add up the three contributions:
    else {
      TH2D * cov_ee = readCov("ee");
      TH2D * cov_emu = readCov("emu");
      TH2D * cov_mumu = readCov("mumu");

      std::vector<Double_t> binning;
      for (Int_t bin = 1; bin <= cov_ee->GetNbinsX(); bin++) { 
	binning.push_back(cov_ee->GetBinLowEdge(bin));
	std::cout << " " << binning.back();
      }
      binning.push_back(cov_ee->GetBinLowEdge(cov_ee->GetNbinsX()) + cov_ee->GetBinWidth(cov_ee->GetNbinsX()));
      std::cout << " " << binning.back() << std::endl;

      // New histo containing the sum:
      cov = new TH2D("cov","cov",cov_ee->GetNbinsX(),&binning.front(),cov_ee->GetNbinsY(),&binning.front());
      for(Int_t y = 0; y <= cov->GetNbinsX(); y++) {
	for(Int_t x = 0; x <= cov->GetNbinsY(); x++) {
	  cov->SetBinContent(x,y,cov_ee->GetBinContent(x,y) + cov_emu->GetBinContent(x,y) + cov_mumu->GetBinContent(x,y));
	}
      }
    }

    massextractor::setStyle(cov,"cov");

    TCanvas * c = new TCanvas("statcov_" + *ch,"statcov_" + *ch);
    c->cd();
    cov->Draw("colz text");
    TGaxis::SetMaxDigits(3);

    //massextractor::DrawCMSLabels();
    massextractor::DrawDecayChLabel(*ch);
    c->Print("Unfolding_" + *ch + "_TtBar_Mass_HypTTBar1stJetMass_ResultCovStat.eps");
  }

  return;
}
