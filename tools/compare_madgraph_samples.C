void DrawDecayChLabel(TString decaychannel, double textSize=0.04) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength()  - 0.325 );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(1.0 - gStyle->GetPadRightMargin() + gStyle->GetTickLength());
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void DrawCMSLabels(int cmsprelim=1, double energy=8, double textSize=0.04) {

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
    label->AddText(Form(text, ((Double_t)19712/1000), energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}

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

void compare_madgraph_samples() {
  TFile mass166("preunfolded/MASS_DOWN_6GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass169("preunfolded/MASS_DOWN_3GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass171("preunfolded/MASS_DOWN_1GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass172("preunfolded/Nominal/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass173("preunfolded/MASS_UP_1GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass175("preunfolded/MASS_UP_3GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");
  TFile mass178("preunfolded/MASS_UP_6GEV/combined/HypTTBar1stJetMass_UnfoldingHistos.root");

  TFile output("comparison_madgraph_masses.root","recreate");

  gDirectory->pwd();

  std::vector<TString> histograms;
  histograms.push_back("aRecHist"); // Invariant Mass TTBar+1Jet

  setHHStyle(*gStyle);

  for(std::vector<TString>::iterator hist = histograms.begin(); hist != histograms.end(); hist++) {

    std::cout << "Plotting " << (*hist) << "..." << std::endl;

    TCanvas *c = new TCanvas(*hist, *hist);
    TH1D *sample_m166 = (TH1D*)mass166.Get(*hist);
    TH1D *sample_m169 = (TH1D*)mass169.Get(*hist);
    TH1D *sample_m171 = (TH1D*)mass171.Get(*hist);
    TH1D *sample_m172 = (TH1D*)mass172.Get(*hist);
    TH1D *sample_m173 = (TH1D*)mass173.Get(*hist);
    TH1D *sample_m175 = (TH1D*)mass175.Get(*hist);
    TH1D *sample_m178 = (TH1D*)mass178.Get(*hist);

    gStyle->SetOptStat(0);
    
    sample_m172->SetLineWidth(3);
    sample_m178->SetFillColor(0);
    sample_m175->SetFillColor(0);
    sample_m173->SetFillColor(0);
    sample_m172->SetFillColor(0);
    sample_m171->SetFillColor(0);
    sample_m169->SetFillColor(0);
    sample_m166->SetFillColor(0);

    sample_m171->SetLineColor(kRed+3);
    sample_m169->SetLineColor(kRed);
    sample_m166->SetLineColor(kOrange);
    sample_m166->SetLineWidth(2);

    sample_m173->SetLineColor(kGreen+3);
    sample_m175->SetLineColor(kGreen);
    sample_m178->SetLineColor(kCyan);
    sample_m178->SetLineWidth(2);
    
    sample_m178->GetYaxis()->SetRangeUser(0,9200);

    sample_m178->Draw();
    sample_m171->Draw("same");
    sample_m169->Draw("same");
    sample_m166->Draw("same");
    sample_m173->Draw("same");
    sample_m175->Draw("same");
    sample_m172->Draw("same");

    if(hist->Contains("aRec")) { 
      sample_m178->GetXaxis()->SetTitle("#rho_{s}");
      sample_m178->GetYaxis()->SetTitle("Events");
    }

    DrawDecayChLabel("Dilepton combined");
    DrawCMSLabels();

    TLegend *leg;
    leg = new TLegend();
    double x1 = 0.215, y1 = 0.6;
    double height = 0.27, width = 0.275;

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    leg->AddEntry(sample_m166,"M = 166.5 GeV","l");
    leg->AddEntry(sample_m169,"M = 169.5 GeV","l");
    leg->AddEntry(sample_m171,"M = 171.5 GeV","l");
    leg->AddEntry(sample_m172,"M = 172.5 GeV","l");
    leg->AddEntry(sample_m173,"M = 173.5 GeV","l");
    leg->AddEntry(sample_m175,"M = 175.5 GeV","l");
    leg->AddEntry(sample_m178,"M = 178.5 GeV","l");
    leg->Draw();

    c->Draw();
    c->Write();
    c->Print((*hist) + ".eps");
    c->Print((*hist) + ".pdf");
  }
}
