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

void setHHStyle(TStyle& HHStyle)
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
void DrawDecayChLabel(TString decaychannel="", Int_t bin=0, int cmsprelim=1, double textSize=0.04) {

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

    if(bin > 0) {
      TPaveText *bintxt = new TPaveText();
      if(!bin_boundaries.empty()) {
	// Lower and upper edge of the bin:
	Double_t bin_low = bin_boundaries.at(bin-1);
	Double_t bin_high = bin_boundaries.at(bin);
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

    TPaveText *decch = new TPaveText();
    decch->AddText(decaychannel);

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

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void DrawCMSLabels(double energy=8, double textSize=0.045) {

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

TString getChannelLabel(TString channel) {
 
  TString label = "";
  if(channel =="ee") { label = "ee"; }
  if(channel =="mumu"){ label = "#mu#mu"; }
  if(channel =="emu"){ label = "e#mu"; }
  if(channel =="combined"){ label = "Dilepton Combined"; }

  return label;
}

TString getSampleLabel(TString systematic) {

  TString label = "";
  if(systematic.Contains("Nominal") || systematic.Contains("MASS")) { label = "Nominal"; }
  else if(systematic.Contains("BTAG")) { label = "B-Tagging"; }
  else if(systematic.Contains("JER")) { label = "Jet Energy Resolution"; }
  else if(systematic.Contains("JES")) { label = "Jet Energy Scale"; }
  else if(systematic.Contains("PU")) { label = "Pile-Up"; }
  else if(systematic.Contains("TRIG")) { label = "Trigger"; }
  else if(systematic.Contains("LEPT")) { label = "Lepton"; }
  else if(systematic.Contains("BG")) { label = "Background"; }
  else if(systematic.Contains("DY")) { label = "Drell-Yan"; }
  else if(systematic.Contains("KIN")) { label = "Kinematic Reconstruction"; }
  else if(systematic.Contains("MATCH")) { label = "Matching"; }
  else if(systematic.Contains("SCALE")) { label = "Scale"; }
  else if(systematic.Contains("HAD")) { label = "Hadronization"; }

  return label;
}

void setStyle(TGraphAsymmErrors *hist)
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


void setStyleAndFillLegend(TGraphErrors* hist, TString name, TLegend *leg) {


  if(name == "data"){
    hist->SetLineWidth(0);
    if(leg) leg->AddEntry(hist, "Pseudo Data",  "p");
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

void setLegendStyle(TLegend *leg)
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

std::string getLastLine(std::ifstream& in)
{
    std::string line;
    while (in >> std::ws && std::getline(in, line)) // skip empty lines
        ;

    return line;
}

void combine_closure() {

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");

  std::vector<TString> masses;
  masses.push_back("166");
  masses.push_back("169");
  masses.push_back("171");
  masses.push_back("172");
  masses.push_back("173");
  masses.push_back("175");
  masses.push_back("178");
  
  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    setHHStyle(*gStyle);

    TGraphAsymmErrors * graph = new TGraphAsymmErrors();
    TGraphAsymmErrors * graphtotal = new TGraphAsymmErrors();

    setStyle(graph);
    setStyle(graphtotal);
    
    for(std::vector<TString>::iterator m = masses.begin(); m != masses.end(); ++m) {
      TString filename = *m + "/MassFitDiffXSecSystematics_" + *ch + ".txt";
      std::ifstream infile(filename);

      if (infile) {
	std::stringstream stream;
	std::string dummy;
  	std::string line = getLastLine(infile);
	std::cout << line << '\n';
	stream.clear();
	stream << line;
	for(size_t i = 0; i < 4; ++i ) { stream >> dummy; }
	double outmass, statp, statn, systp, systn;
	stream >> outmass >> statp >> statn >> dummy >> systp >> systn;
	//stat = atof(dummy.erase(0,2).c_str());
	double inmass = (atof(*m)+0.5);

	double tot_pos = sqrt(statp*statp + systp*systp);
	double tot_neg = sqrt(statn*statn + systn*systn);

	std::cout << outmass << "@" << inmass << " | " 
		  << statp << " " << statn << " | " 
		  << systp << " " << systn << " | "
		  << tot_pos << " " << tot_neg << std::endl;

	// Add point:
	graph->SetPoint(m-masses.begin(),inmass,outmass);
	graphtotal->SetPoint(m-masses.begin(),inmass,outmass);

	//graphtotal->SetPointError(m-masses.begin(),0,0,tot_neg,tot_pos);
	graphtotal->SetPointError(m-masses.begin(),0,0,TMath::Abs(statn),TMath::Abs(statp));
	graph->SetPointError(m-masses.begin(),0,0,TMath::Abs(statn),TMath::Abs(statp));
      }
      else
        std::cout << "Unable to open file " << filename << ".\n";

      infile.close();
    }

    graphtotal->GetXaxis()->SetTitle("input m_{t} #left[GeV#right]");
    graphtotal->GetYaxis()->SetTitle("output m_{t} #left[GeV#right]");

    Double_t m_min = 159, m_max = 187, m_nominal = 172.5;
    graphtotal->GetXaxis()->SetLimits(m_min,m_max);
    graphtotal->GetYaxis()->SetRangeUser(m_min,m_max);

    graph->SetLineColor(kGray+1);
    graph->SetMarkerStyle(20);
    graphtotal->SetLineColor(kOrange-4);

    TLine *l = new TLine(m_min,m_min,m_max,m_max);
    l->SetLineColor(kGray+1);
    l->SetLineStyle(7);
    l->SetLineWidth(2);

    TLine *l2 = new TLine(m_min,m_nominal,m_max,m_nominal);
    l2->SetLineColor(kGray+1);
    l2->SetLineStyle(7);
    l2->SetLineWidth(2);

    TCanvas * c = new TCanvas("diffxs_closure_unfold_" + *ch,"diffxs_closure_unfold_" + *ch);
    c->cd();
    graphtotal->Draw("a p e1");
    graph->Draw("same e1");
    l->Draw();
    l2->Draw();
    DrawCMSLabels();
    DrawDecayChLabel(*ch);
    c->Print("diffxs_closure_unfold_" + *ch + ".pdf");
  }

  return;
}
