#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <Riostream.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TString.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

Double_t nominalmass = 172.5;
Double_t lumi = 19712;

TH2D * readCorr(TString ch) {

  TString filename = "SVD/Nominal/Unfolding_" + ch + "_TtBar_Mass_HypTTBar1stJetMass.root";
  TString histogramname = "SVD_" + ch + "_TtBar_Mass_HypTTBar1stJetMass_Nominal_STATCORRORM";
  
  TFile * file = TFile::Open(filename,"read");
  if(!file->IsOpen()) {
    std::cout << "Failed to open file " << filename << std::endl;
    throw 1;
  }
  else std::cout << "Opened file " << filename << std::endl;
  
  // Histogram containing differential cross section from data:
  TH2D * statCorrNorm = static_cast<TH2D*>(file->Get(histogramname));
  if(!statCorrNorm) { std::cout << "Failed to get histogram " << histogramname << std::endl; }
  else { std::cout << "Fetched " << histogramname << std::endl; }
  
  std::vector<Double_t> binning;
  for (Int_t bin = 2; bin < statCorrNorm->GetNbinsX(); bin++) { 
    binning.push_back(statCorrNorm->GetBinLowEdge(bin));
      std::cout << " " << binning.back();
  }
  binning.push_back(statCorrNorm->GetBinLowEdge(statCorrNorm->GetNbinsX()-1) + statCorrNorm->GetBinWidth(statCorrNorm->GetNbinsX()-1));
  std::cout << " " << binning.back() << std::endl;
  
  // Cut away over and underflow bins:
  TH2D * corr = new TH2D("corr","corr",statCorrNorm->GetNbinsX()-2,&binning.front(),statCorrNorm->GetNbinsY()-2,&binning.front());
  for(Int_t y = 0; y <= corr->GetNbinsX(); y++) {
    for(Int_t x = 0; x <= corr->GetNbinsY(); x++) {
      corr->SetBinContent(x,y,statCorrNorm->GetBinContent(x+1,y+1));
    }
  }

  return corr;
}

void plotcorr() {

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

    TH2D *corr;
    if(*ch != "combined") corr = readCorr(*ch);
    // Combined channel - add up the three contributions:
    else {
      TH2D * corr_ee = readCorr("ee");
      TH2D * corr_emu = readCorr("emu");
      TH2D * corr_mumu = readCorr("mumu");

      std::vector<Double_t> binning;
      for (Int_t bin = 1; bin <= corr_ee->GetNbinsX(); bin++) { 
	binning.push_back(corr_ee->GetBinLowEdge(bin));
	std::cout << " " << binning.back();
      }
      binning.push_back(corr_ee->GetBinLowEdge(corr_ee->GetNbinsX()) + corr_ee->GetBinWidth(corr_ee->GetNbinsX()));
      std::cout << " " << binning.back() << std::endl;

      // New histo containing the sum:
      corr = new TH2D("corr","corr",corr_ee->GetNbinsX(),&binning.front(),corr_ee->GetNbinsY(),&binning.front());
      for(Int_t y = 0; y <= corr->GetNbinsX(); y++) {
	for(Int_t x = 0; x <= corr->GetNbinsY(); x++) {
	  corr->SetBinContent(x,y,(corr_ee->GetBinContent(x,y) + corr_emu->GetBinContent(x,y) + corr_mumu->GetBinContent(x,y))/3);
	}
      }
    }

    massextractor::setStyle(corr,"corr");
    gStyle->SetPaintTextFormat("2.1f");

    TCanvas * c = new TCanvas("statcorr_" + *ch,"statcorr_" + *ch);
    c->cd();
    //corr->GetZaxis()->SetTitle("stat. corr.");
    //corr->Draw("colz text");
    corr->Draw("text");
    TGaxis::SetMaxDigits(3);

    //superimpose lines at the xbins positions
    std::vector<Double_t> binning;
    for (Int_t bin = 2; bin < corr->GetNbinsX(); bin++) { 
      binning.push_back(corr->GetBinLowEdge(bin));
      std::cout << " " << binning.back();
    }
    binning.push_back(corr->GetBinLowEdge(corr->GetNbinsX()-1) + corr->GetBinWidth(corr->GetNbinsX()-1));
    std::cout << " " << binning.back() << std::endl;

    TLine l;
    c->Update();
    Double_t ymin = c->GetUymin();
    Double_t ymax = c->GetUymax();
    Double_t xmin = c->GetUxmin();
    Double_t xmax = c->GetUxmax();
    l.SetLineStyle(3);
    for (Int_t bin = 0; bin < binning.size();bin++) {
      l.DrawLine(binning.at(bin),ymin,binning.at(bin),ymax);
      l.DrawLine(xmin,binning.at(bin),xmax,binning.at(bin));
    }

    //massextractor::DrawFreeCMSLabels(Form("%s, %2.1f fb^{-1} (8 TeV)",*ch,lumi/1000),0.045);
    if(*ch == "ee") massextractor::DrawFreeCMSLabels(Form("ee, %2.1f fb^{-1} (8 TeV)",lumi/1000),0.045);
    else if(*ch == "emu") massextractor::DrawFreeCMSLabels(Form("e#mu, %2.1f fb^{-1} (8 TeV)",lumi/1000),0.045);
    else if(*ch == "mumu") massextractor::DrawFreeCMSLabels(Form("#mu#mu, %2.1f fb^{-1} (8 TeV)",lumi/1000),0.045);
    else massextractor::DrawFreeCMSLabels(Form("Dilepton combined, %2.1f fb^{-1} (8 TeV)",lumi/1000),0.045);
    //massextractor::DrawDecayChLabel(*ch);
    c->Print("Unfolding_" + *ch + "_HypTTBar1stJetMass_CorrStat.pdf");
  }

  return;
}
