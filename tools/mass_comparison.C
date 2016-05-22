#include "TGraphAsymmErrors.h"
#include "plotter.C"

#define LEFT 155.0
#define RIGHT 186.0

TGraphAsymmErrors * getGraph(Int_t i, Double_t m, Double_t errp, Double_t errm) {

  Double_t x[1]   = {m};
  Double_t y[1]   = {static_cast<Double_t>(i)};
  Double_t exl[1] = {errm};
  Double_t eyl[1] = {0};
  Double_t exh[1] = {errp};
  Double_t eyh[1] = {0};
  TGraphAsymmErrors * gr = new TGraphAsymmErrors(1,x,y,exl,exh,eyl,eyh);
  gr->SetLineWidth(2);
  gr->GetXaxis()->SetTitle("m_{t}^{pole} [GeV]");
  gr->GetYaxis()->SetTickLength(0);
  gr->GetYaxis()->SetLabelOffset(999);
  return gr;
}

TPaveText * getName(Int_t i, TString text) {
  TPaveText *label = new TPaveText(LEFT+1,(i+0.15),LEFT+10,(i+0.15));
  label->SetTextFont(42);
  label->SetTextSize(0.035);
  label->AddText(text);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextAlign(12);
  return label;
}

TPaveText * getJournal(Int_t i, TString text) {
  TPaveText *label = new TPaveText(LEFT+1,(i-0.15),LEFT+10,(i-0.15));
  label->SetTextFont(42);
  label->SetTextSize(0.023);
  label->AddText(text);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextAlign(12);
  return label;
}

TPaveText * getLabel(Int_t i, Double_t mass, Double_t errp, Double_t errm) {
  TPaveText *label = new TPaveText(RIGHT-10,static_cast<Double_t>(i),RIGHT-1,static_cast<Double_t>(i));
  label->SetTextFont(42);
  label->SetTextSize(0.03);
  label->AddText(Form("%3.1f ^{+%1.1f}_{-%1.1f} GeV",mass,errp,errm));
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextAlign(32);
  return label;
}

void drawMeasurement(Int_t i, TString name, TString journal, Double_t mass, Double_t totp, Double_t totm, Double_t statp = 0, Double_t statm = 0, bool mark = false, bool first = false) {

  TGraphAsymmErrors *gr = getGraph(static_cast<Double_t>(i),mass,totp,totm);

  if(mark) {
    gr->SetLineColor(kRed+1);
    gr->SetMarkerColor(kRed+1);
  }
  
  if(first) {
    cout << "first";
    gr->GetXaxis()->SetLimits(LEFT-10,RIGHT+10);
    gr->Draw("ALP");
    gr->GetXaxis()->SetRangeUser(LEFT,RIGHT);
    gr->GetYaxis()->SetRangeUser(0,8);
  }
  else { gr->Draw("LP"); }

  if(statp != 0 && statm != 0) {
    gr = getGraph(static_cast<Double_t>(i),mass,statp,statm);
    gr->Draw("LP");
  }

  TPaveText *tmass = getLabel(static_cast<Double_t>(i),mass,totp,totm);
  tmass->Draw("same");
  TPaveText *mname = getName(static_cast<Double_t>(i),name);
  mname->Draw("same");

  if(journal != "") {
    TPaveText *label = getJournal(static_cast<Double_t>(i),journal);
    label->Draw("same");
  }
  else {
    mname->SetY1(i);
    mname->SetY2(i);
  }

  if(mark) {
    gr->SetLineColor(kRed+1);
    gr->SetMarkerColor(kRed+1);
    tmass->SetTextColor(kRed+1);
    mname->SetTextColor(kRed+1);
  }
}

void mass_comparison() {

  // Set the histogram styles:
  setHHStyle(*gStyle);

  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  TCanvas *c = new TCanvas("c", "c",2277,231,1000,700);

  drawMeasurement(7, "D0 5.3 fb^{-1}", "PLB 703 (2011) 422", 167.5, 5.4, 4.9, 0, 0, false, true);
  drawMeasurement(6, "CMS 7 TeV", "PLB 738 (2014) 526", 176.7, 3.0, 2.8);
  drawMeasurement(5, "ATLAS 7+8 TeV", "EPJ C74 (2014) 3109", 172.9, 2.5, 2.6);
  drawMeasurement(4, "ATLAS tt+jet", "JHEP 10 (2015) 121", 173.7, 2.3, 2.1, 1.5, 1.5);
  drawMeasurement(3, "D0 9.7 fb^{-1} Preliminary", "D0 6453-CONF", 169.5, 3.3, 3.4);
  drawMeasurement(2, "CMS 7+8 TeV Preliminary", "PAS-TOP-13-004", 173.6, 1.7, 1.8);
  drawMeasurement(1, "This Measurement", "", 168.2, 4.7, 2.1, 1.1, 1.1, true);

}
