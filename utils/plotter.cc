#include <TPaveText.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TF1.h>

#include <iostream>

#include "plotter.h"

using namespace massextractor;

void massextractor::setHHStyle(TStyle& HHStyle)
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
    HHStyle.SetHistLineWidth(2);
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
void massextractor::DrawDecayChLabel(TString decaychannel, bool drawchannel, Int_t bin, std::vector<Double_t> boundaries, int cmsprelim, double textSize) {
  /*
    TPaveText *cms = new TPaveText();
    cms->AddText("CMS");

    cms->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    cms->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    cms->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    cms->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    cms->SetFillStyle(0);
    cms->SetBorderSize(0);
    if (textSize!=0) cms->SetTextSize(textSize*1.1);
    cms->SetTextAlign(12);
    cms->SetTextFont(61);
    cms->Draw("same");

    if(cmsprelim > 0) {
      TPaveText *extra = new TPaveText();
      if(cmsprelim == 2) { extra->AddText("Private Work"); }
      else { extra->AddText("Preliminary"); }

      extra->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
      extra->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10 );
      extra->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
      extra->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );

      extra->SetFillStyle(0);
      extra->SetBorderSize(0);
      if (textSize!=0) extra->SetTextSize(textSize);
      extra->SetTextAlign(12);
      extra->SetTextFont(52);
      extra->Draw("same");
    }
  */
    if(bin > 0) {
      TPaveText *bintxt = new TPaveText();
      if(!boundaries.empty()) {
	// Lower and upper edge of the bin:
	Double_t bin_low = boundaries.at(bin-1);
	Double_t bin_high = boundaries.at(bin);
	bintxt->AddText(Form("%1.2f < #rho_{S} < %1.2f",bin_low,bin_high));
      }
      else {
	bintxt->AddText(Form("bin %i",bin));
      }
      bintxt->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
      bintxt->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10 );
      bintxt->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength()       );
      bintxt->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );

      bintxt->SetFillStyle(0);
      bintxt->SetBorderSize(0);
      if (textSize!=0) bintxt->SetTextSize(textSize);
      bintxt->SetTextAlign(32);
      bintxt->Draw("same");
    }

    if(drawchannel) {
      TPaveText *decch = new TPaveText();
      decch->AddText(getChannelLabel(decaychannel));

      decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
      decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
      decch->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength()       );
      decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

      decch->SetFillStyle(0);
      decch->SetBorderSize(0);
      if (textSize!=0) decch->SetTextSize(textSize);
      decch->SetTextAlign(32);
      decch->Draw("same");
    }
}

void massextractor::DrawFreeCMSLabels(char* text, double textSize) {

  //const char *text = "%2.1f #times 10^{6} clusters (fiducial) (%1.1f GeV)";
  //char *text = "%2.1f #times 10^{6} clusters (fiducial)";
    
  TPaveText *label = new TPaveText();
  label->SetX1NDC(gStyle->GetPadLeftMargin());
  label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label->SetY2NDC(1.0);
  label->SetTextFont(42);
  label->AddText(text);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  if (textSize!=0) label->SetTextSize(textSize);
  label->SetTextAlign(32);
  label->Draw("same");
}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void massextractor::DrawCMSLabels(double lumi, double energy, double textSize) {

    const char *text = "%2.1f fb^{-1} (%2.f TeV)";
    
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

void massextractor::setStyle(TGraphErrors *hist, TString name)
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

  if(name == "data") {
    hist->SetLineWidth(2);
    hist->SetFillStyle(3005);
    hist->SetFillColor(kBlack);
  }
  else if(name == "madgraph" || name == "powheg") {
    hist->SetMarkerColor(kRed+1);
    hist->SetLineColor(kRed+1);
    hist->SetFillStyle(3004);
    hist->SetFillColor(kRed+1);
  }
}

void massextractor::setTheoryStyleAndFillLegend(TH1* histo, TString theoryName, TLegend *leg) {

    histo->GetXaxis()->SetTitleOffset(1.08);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelOffset(0.007);
    histo->GetXaxis()->SetLabelSize(0.04);

    histo->GetYaxis()->SetTitleOffset(1.7);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelOffset(0.007);
    histo->GetYaxis()->SetLabelSize(0.04);

    histo->SetLineWidth(2);
    if(theoryName != "data"){
        histo->SetMarkerSize(0);
        histo->SetMarkerStyle(1);
    }

    if(theoryName == "madgraph"){
        histo->SetLineColor(kRed+1);
        histo->SetLineStyle(1);
        if(leg) leg->AddEntry(histo, "MadGraph+Pythia",  "l");
    }
    if(theoryName == "powheg") {
        histo->SetLineColor(kRed+1);
        histo->SetLineStyle(1);
        if(leg) leg->AddEntry(histo, "Powheg+Pythia8",  "l");
    }
    if(theoryName == "powhegpythia"){
      histo->SetLineColor(kGreen+1);
      histo->SetLineStyle(7);
      if(leg) leg->AddEntry(histo, "Powheg+Pythia",  "l");
    }
    if(theoryName == "powhegbox") {
        histo->SetLineColor(kBlue);
        histo->SetLineStyle(5);
        if(leg) leg->AddEntry(histo, "PowhegBox",  "l");
    }
    if(theoryName == "powheg2") {
        histo->SetLineColor(kRed-7);
        histo->SetLineStyle(5);
        if(leg) leg->AddEntry(histo, "Powheg (hdamp)",  "l");
    }
    if(theoryName == "powheg2box") {
        histo->SetLineColor(kCyan+2);
        histo->SetLineStyle(1);
        if(leg) leg->AddEntry(histo, "PowhegBox (hdamp)",  "l");
    }
    if(theoryName == "powheg2r11box") {
        histo->SetLineColor(kRed+2);
        histo->SetLineStyle(1);
        if(leg) leg->AddEntry(histo, "PowhegBox (hdamp) r11",  "l");
    }
}

void massextractor::setStyle(TH1D *hist, TString name)
{
  hist->SetLineWidth(1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.08);
  hist->GetZaxis()->SetTitleOffset(1.08);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelOffset(0.007);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetTitle("#rho_{s}");
  hist->GetYaxis()->SetTitle("#rho_{s}");
}

void massextractor::setStyle(TH2D *hist, TString name)
{
  hist->SetLineWidth(1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.08);
  hist->GetZaxis()->SetTitleOffset(1.08);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelOffset(0.007);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.2);
  hist->SetLineWidth(2);
  
  hist->GetXaxis()->SetTitle("#rho_{s}");
  hist->GetYaxis()->SetTitle("#rho_{s}");

  if(name == "cov") {
    hist->SetLineWidth(2);
    hist->SetFillStyle(3005);
    hist->SetFillColor(kBlack);
    hist->GetZaxis()->SetTitle("Statisical Covariance (Absolute)");
  }
}

void massextractor::setStyleAndFillLegend(TGraphErrors* hist, TString name, TLegend *leg, bool closure) {

  setStyle(hist,name);

  if(name == "data") {
    if(leg && closure) leg->AddEntry(hist, "Pseudo Data",  "l");
    else if(leg) leg->AddEntry(hist, "CMS Data",  "l");
  }
  else if(name == "madgraph") {
    if(leg) leg->AddEntry(hist, "MadGraph+Pythia",  "p");
  }
  else if(name == "powheg") {
    if(leg) leg->AddEntry(hist, "Powheg+Pythia8",  "p");
  }
}

void massextractor::setLegendStyle(TLegend *leg) {

  Double_t x1 = 0.560;
  Double_t y1 = 1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.2;
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

void massextractor::rescaleGraph(TGraphErrors * g, Double_t up, Double_t down) {
  
  Double_t low = TMath::MinElement(g->GetN(),g->GetY())*down;
  Double_t high = TMath::MaxElement(g->GetN(),g->GetY())*up;
  g->GetYaxis()->SetRangeUser(low,high);
}

TString massextractor::getChannelLabel(TString channel) {
 
  TString label = channel;
  if(channel =="ee") { label = "ee"; }
  if(channel =="mumu"){ label = "#mu#mu"; }
  if(channel =="emu"){ label = "e#mu"; }
  if(channel =="combined"){ label = "Dilepton Combined"; }

  return label;
}

void massextractor::drawRatio(const TH1* histNumerator, const TH1* histDenominator1, const TH1* uncband, 
			      TGraphAsymmErrors *ratio_stat, TGraphAsymmErrors *ratio_total, 
			      const TH1* histDenominator2, const TH1* histDenominator3, 
			      const TH1* histDenominator4, const TH1* histDenominator5, 
			      const TH1* histDenominator6, const TH1* histDenominator7, 
			      const Double_t& ratioMin, const Double_t& ratioMax, const bool data,
			      TStyle myStyle) {
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator_i' (_i = 1-5)
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // the histogram style is transferred from 'histDenominator_i' to the 'ratio_i'
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    /// check that histos exist and have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator1->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        return;
    }

    /// create ratio
    TH1F *ratio1 = 0, *ratio2 = 0, *ratio3 = 0, *ratio4 = 0, *ratio5 = 0, *ratio6 = 0, *ratio7 = 0; 

    ratio1 = (TH1F*)histDenominator1->Clone();
    ratio1->SetLineColor(histDenominator1->GetLineColor());
    ratio1->SetLineStyle(histDenominator1->GetLineStyle());
    ratio1->SetLineWidth(histDenominator1->GetLineWidth());
    ratio1->Divide(histNumerator);

    if (histDenominator2){
        ratio2 = (TH1F*)histDenominator2->Clone();
        ratio2->SetLineColor(histDenominator2->GetLineColor());
        ratio2->SetLineStyle(histDenominator2->GetLineStyle());
        ratio2->SetLineWidth(histDenominator2->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator2->GetNbinsX()){ratio2 = 0;}
        else {ratio2->Divide(histNumerator);}
    };
    if (histDenominator3){
        ratio3 = (TH1F*)histDenominator3->Clone();
        ratio3->SetLineColor(histDenominator3->GetLineColor());
        ratio3->SetLineStyle(histDenominator3->GetLineStyle());
        ratio3->SetLineWidth(histDenominator3->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator3->GetNbinsX()){ratio3 = 0;}
        else {ratio3->Divide(histNumerator);}
    };
    if (histDenominator4){
        ratio4 = (TH1F*)histDenominator4->Clone();
        ratio4->SetLineColor(histDenominator4->GetLineColor());
        ratio4->SetLineStyle(histDenominator4->GetLineStyle());
        ratio4->SetLineWidth(histDenominator4->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator4->GetNbinsX()){ratio4 = 0;}
        else {ratio4->Divide(histNumerator);}
    };
    if (histDenominator5){
        ratio5 = (TH1F*)histDenominator5->Clone();
        ratio5->SetLineColor(histDenominator5->GetLineColor());
        ratio5->SetLineStyle(histDenominator5->GetLineStyle());
        ratio5->SetLineWidth(histDenominator5->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator5->GetNbinsX()){ratio5 = 0;}
        else {ratio5->Divide(histNumerator);}
    };
    if (histDenominator6){
        ratio6 = (TH1F*)histDenominator6->Clone();
        ratio6->SetLineColor(histDenominator6->GetLineColor());
        ratio6->SetLineStyle(histDenominator6->GetLineStyle());
        ratio6->SetLineWidth(histDenominator6->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator6->GetNbinsX()){ratio6 = 0;}
        else {ratio6->Divide(histNumerator);}
    };
    if (histDenominator7){
        ratio7 = (TH1F*)histDenominator7->Clone();
        ratio7->SetLineColor(histDenominator7->GetLineColor());
        ratio7->SetLineStyle(histDenominator7->GetLineStyle());
        ratio7->SetLineWidth(histDenominator7->GetLineWidth());
        if(histNumerator->GetNbinsX()!=histDenominator7->GetNbinsX()){ratio7 = 0;}
        else {ratio7->Divide(histNumerator);}
    };

    /// calculate error for ratio only gaussian error of histNumerator
    for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
        if (ratio1) ratio1->SetBinError(bin, sqrt(histDenominator1->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio2) ratio2->SetBinError(bin, sqrt(histDenominator2->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio3) ratio3->SetBinError(bin, sqrt(histDenominator3->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio4) ratio4->SetBinError(bin, sqrt(histDenominator4->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio5) ratio5->SetBinError(bin, sqrt(histDenominator5->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio6) ratio6->SetBinError(bin, sqrt(histDenominator6->GetBinContent(bin))/histNumerator->GetBinContent(bin));
        if (ratio7) ratio7->SetBinError(bin, sqrt(histDenominator7->GetBinContent(bin))/histNumerator->GetBinContent(bin));
    }

    Int_t    logx  = myStyle.GetOptLogx();
    Double_t left  = myStyle.GetPadLeftMargin();
    Double_t right = myStyle.GetPadRightMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(0);

    /// create new pad for ratio plot
    TPad *rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->SetFillColor(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    
    /// configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio1->SetStats(kFALSE);
    ratio1->SetTitle("");
    ratio1->SetMaximum(ratioMax);
    ratio1->SetMinimum(ratioMin);
    
    /// configure axis of ratio plot
    ratio1->GetXaxis()->SetTitleSize(histNumerator->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
    ratio1->GetXaxis()->SetTitleOffset(histNumerator->GetXaxis()->GetTitleOffset()*0.9);
    ratio1->GetXaxis()->SetLabelSize(histNumerator->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
    ratio1->GetXaxis()->SetTitle(histNumerator->GetXaxis()->GetTitle());


    ratio1->GetYaxis()->CenterTitle();
    if(data) ratio1->GetYaxis()->SetTitle("#frac{Theory}{Data}");
    else ratio1->GetYaxis()->SetTitle("#frac{Theory}{MadGraph}");
    ratio1->GetYaxis()->SetTitleSize(histNumerator->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio1->GetYaxis()->SetTitleOffset(histNumerator->GetYaxis()->GetTitleOffset()/scaleFactor);
    ratio1->GetYaxis()->SetLabelSize(histNumerator->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio1->GetYaxis()->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset());
    ratio1->GetYaxis()->SetTickLength(0.03);
    ratio1->GetYaxis()->SetNdivisions(505);
    ratio1->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    ratio1->GetXaxis()->SetNoExponent(kTRUE);
    
    /// delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    histNumerator->GetXaxis()->SetNoExponent(kFALSE);

    /// draw ratio plot
    ratio1->DrawClone("Histo");
    rPad->Update();
    rPad->Modified();
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(0.15*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    
    /// draw grid
    rPad->SetGrid(0,1);

    // draw a horizontal lines on a given histogram
    // a) at 1
    Double_t xmin = ratio1->GetXaxis()->GetXmin();
    Double_t xmax = ratio1->GetXaxis()->GetXmax();
    TString height = ""; height += 1;
    TF1 *f = new TF1("f", height.Data(), xmin, xmax);
    f->SetLineStyle(1);//this is frustrating and stupid but apparently necessary...
    f->SetLineWidth(1);
    f->SetLineColor(kBlack);
    f->Draw("L same");
    // b) at upper end of ratio pad
    TString height2 = ""; height2 += ratioMax;
    TF1 *f2 = new TF1("f2", height2.Data(), xmin, xmax);
    f2->SetLineStyle(1);
    f2->SetLineWidth(1);
    f2->SetLineColor(kBlack);
    f2->Draw("L same");

    // create ratio of uncertainty band
    TH1 *band = nullptr;
    if (uncband) {
        band = (TH1*)uncband->Clone("band");
        for(int i=0; i<= 1+uncband->GetNbinsX(); i++) {
            double error = 0;
            double content = 1;
            if(band->GetBinContent(i)) {
                error = band->GetBinError(i) / band->GetBinContent(i);
                content = band->GetBinContent(i) /band->GetBinContent(i);
            }
            band->SetBinError(i, error/content);
            band->SetBinContent(i, content/content);
        }
    }
    if(band) band->Draw("same,e2");

    if(ratio_stat) {
        TLegend *leg_band = new TLegend();
        if(ratio_total) leg_band->AddEntry(ratio_total, "Stat. #oplus Syst.", "f");
        leg_band->AddEntry(ratio_stat, "Stat.", "f");
        leg_band->SetX1NDC(0.22);
        leg_band->SetY1NDC(0.97);
        leg_band->SetX2NDC(0.46);
        leg_band->SetY2NDC(0.77);
        leg_band->SetFillStyle(1001);
        leg_band->SetFillColor(10);
        leg_band->SetBorderSize(0);
        leg_band->SetTextSize(0.1);
        leg_band->SetTextAlign(12);
        leg_band->Draw("same");
        if(ratio_total)ratio_total->Draw("same,e2");
        ratio_stat->Draw("same,e2");
    }
    f->Draw("l,same");
    f2->Draw("l,same");
    gPad->RedrawAxis();
    gPad->Update();
    gPad->Modified();
    ratio1->Draw("histo,same");
    if (ratio2) ratio2->Draw("Histo,same");
    if (ratio3) ratio3->Draw("Histo,same");
    if (ratio4) ratio4->Draw("Histo,same");
    if (ratio5) ratio5->Draw("Histo,same");
    if (ratio6) ratio6->Draw("Histo,same");
    if (ratio7) ratio7->Draw("Histo,same");
    gPad->RedrawAxis();
    gPad->Update();
    gPad->Modified();
    rPad->RedrawAxis();
}
