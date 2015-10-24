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

#include "plotter.h"

TH1D * getNiceHistogram(int nbins, std::vector<Double_t> binning, TFile * file, TString hist) {
  TH1D *sample = (TH1D*)(file->Get(hist)->Clone());
  TH1D *binned = dynamic_cast<TH1D*>(sample->Rebin(nbins,"madgraph",&binning.front()));
  for (Int_t bin=0; bin < nbins; bin++) {
    // Divide rebinned histogram's bin content by bin width factor (new/old):
    binned->SetBinError(bin+1,sqrt(binned->GetBinContent(bin+1))/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
    binned->SetBinContent(bin+1,binned->GetBinContent(bin+1)/((binning.at(bin+1)-binning.at(bin))/sample->GetBinWidth(1)));
  }
  binned->Scale(1/binned->Integral("width"));
  return binned;
}

int main() {

  massextractor::setHHStyle(*gStyle);

  // Compare MADGRAPH with DATA, POWHEG_CMS, POWHEGv2_CMS, POWHEG_BOX, POWHEGv2_BOX
  TFile * data = new TFile("/nfs/dust/cms/user/spanns/N010_Carmen/UnfoldingResults/Nominal/emu/HypTTBar1stJetMassResults.root");

  TFile * nominal = new TFile("/nfs/dust/cms/user/spanns/N010_Carmen/selectionRoot/Nominal/emu/emu_ttbarsignalplustau.root");
  TFile * pwhg = new TFile("/nfs/dust/cms/user/spanns/N010_Carmen/selectionRoot/POWHEG/emu/emu_ttbarsignalplustau_powheg.root");
  TFile * pwhg2 = new TFile("/nfs/dust/cms/user/spanns/N010_Carmen_experiments/selectionRoot/POWHEGV2/emu/emu_ttbarsignalplustau_powhegv2.root");

  TFile * pwhgbx = new TFile("/nfs/dust/cms/user/spanns/TopMass_SelectionRoot_outdated/ResultsPseudo172/selectionRoot/POWHEG/emu/emu_ttbarsignalplustau_pwhgbox_m172p5.root");
  TFile * pwhgbx2 = new TFile("/nfs/dust/cms/user/spanns/N010_PowhegBox/selectionRoot/POWHEG/emu/emu_ttbarsignalplustau_powhegbox_r9_m172p5.root");

  TFile * pwhgbx2r11 = new TFile("/nfs/dust/cms/user/spanns/N010_PowhegBox/selectionRoot/POWHEG/emu/emu_ttbarsignalplustau_powhegbox_m172p5.root");


  TFile output("comparison_pwhg_data.root","recreate");
  gDirectory->pwd();

  std::vector<TString> histograms;

  histograms.push_back("VisGenTTBarpT"); // TTBar system pT
  histograms.push_back("VisGenTTBarMass"); // TTBar system Mass
  histograms.push_back("VisGenTTBarRapidity"); // TTBar Rapidity
  histograms.push_back("VisGenLeptonpT"); // Lepton pT
  histograms.push_back("VisGenLeptonEta"); // Lepton Eta
  histograms.push_back("VisGenBJetpT"); // BJet pT
  histograms.push_back("VisGenBJetRapidity"); // BJet Rapidity
  histograms.push_back("VisGenExtraJetpT"); // Additional Jet pT
  histograms.push_back("VisGenTTBar1stJetMass"); // Invariant Mass TTBar+1Jet  


  for(std::vector<TString>::iterator hist = histograms.begin(); hist != histograms.end(); hist++) {

    std::cout << "Plotting " << (*hist) << "..." << std::endl;

    TH1D *sample_data = (TH1D*)(nominal->Get("VisGenTTBarpT")->Clone());
    //    sample_data->Scale(1/sample_data->Integral("width"));    
    // Get binning from dat:
    std::vector<Double_t> binning;

    int nbins = sample_data->GetNbinsX();
    for (Int_t bin = 1; bin <= nbins; bin++) { 
      binning.push_back(sample_data->GetBinLowEdge(bin));
    }
    binning.push_back(sample_data->GetBinLowEdge(nbins) + sample_data->GetBinWidth(nbins));
    std::cout << "Rebinning with " << binning.size()-1 << " bins." << std::endl;

    TH1D *nominalBinned = getNiceHistogram(nbins, binning, nominal, *hist);

    TH1D *pwhgBinned = getNiceHistogram(nbins, binning, pwhg, *hist);
    TH1D *pwhg2Binned = getNiceHistogram(nbins, binning, pwhg2, *hist);

    TH1D *pwhgbxBinned = getNiceHistogram(nbins, binning, pwhgbx, *hist);
    TH1D *pwhgbx2Binned = getNiceHistogram(nbins, binning, pwhgbx2, *hist);

    TH1D *pwhgbx2r11Binned = getNiceHistogram(nbins, binning, pwhgbx2r11, *hist);

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas(*hist, "");

    TLegend *leg;
    if(hist->Contains("Rapidity")) { leg = new TLegend(0.15,0.6,0.4,.85); }
    else { leg = new TLegend(0.5,0.6,0.85,.85); }

    massextractor::setStyle(sample_data);
    sample_data->SetLineColor(kBlack);
    sample_data->SetLineWidth(3);
    sample_data->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#rho_{s}}");
    //sample_data->Draw();
    sample_data->Write();
    leg->AddEntry(sample_data,Form("Data",sample_data->GetMean()),"l");

    massextractor::setTheoryStyleAndFillLegend(nominalBinned, "madgraph", leg);
    nominalBinned->Draw();
    nominalBinned->Write();
    //leg->AddEntry(nominalBinned,Form("MadGraph Nominal (Mean %1.2f)",nominalBinned->GetMean()),"l");

    if(hist->Contains("Eta")) nominalBinned->GetXaxis()->SetTitle("#eta");
    if(hist->Contains("pT")) nominalBinned->GetXaxis()->SetTitle("p_{T} [GeV]");
    if(hist->Contains("Mass")) nominalBinned->GetXaxis()->SetTitle("m_{tt#bar} [GeV]");
    if(hist->Contains("1stJetMass")) nominalBinned->GetXaxis()->SetTitle("m_{t#bar{t}+Jet} [GeV]");
    if(hist->Contains("Rapidity")) nominalBinned->GetXaxis()->SetTitle("Rapidity");

    massextractor::setTheoryStyleAndFillLegend(pwhgBinned, "powhegpythia", leg);
    pwhgBinned->Draw("same");
    pwhgBinned->Write();

    massextractor::setTheoryStyleAndFillLegend(pwhgbxBinned, "powhegbox", leg);
    pwhgbxBinned->Draw("same");
    pwhgbxBinned->Write();

    massextractor::setTheoryStyleAndFillLegend(pwhg2Binned, "powheg2", leg);
    pwhg2Binned->Draw("same");
    pwhg2Binned->Write();

    massextractor::setTheoryStyleAndFillLegend(pwhgbx2Binned, "powheg2box", leg);
    pwhgbx2Binned->Draw("same");
    pwhgbx2Binned->Write();

    massextractor::setTheoryStyleAndFillLegend(pwhgbx2r11Binned, "powheg2r11box", leg);
    pwhgbx2r11Binned->Draw("same");
    pwhgbx2r11Binned->Write();
    
    massextractor::setLegendStyle(leg);
    leg->Draw();

    massextractor::drawRatio(sample_data, nominalBinned, NULL, NULL, 
			     pwhgBinned, pwhgbxBinned, 
			     pwhg2Binned, pwhgbx2Binned,
			     pwhgbx2r11Binned, NULL,
			     0.3,3);

    c->Update();
    c->Modified();

    c->Write();
  }

  return 0;
}
