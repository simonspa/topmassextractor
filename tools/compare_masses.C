void compare_masses() {
  TFile mass166("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m166p5.root");
  TFile mass169("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m169p5.root");
  TFile mass171("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m171p5.root");
  TFile mass172("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m172p5.root");
  TFile mass173("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m173p5.root");
  TFile mass175("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m175p5.root");
  TFile mass178("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m178p5.root");

  TFile output("comparison_pwhg_masses.root","recreate");

  gDirectory->pwd();

  std::vector<TString> histograms;

  histograms.push_back("VisGenTTBarpT"); // TTBar system pT
  histograms.push_back("VisGenTTBarMass"); // TTBar system Mass
  histograms.push_back("VisGenTTBarRapidity"); // TTBar Rapidity
  histograms.push_back("VisGenLeptonpT"); // Lepton pT
  histograms.push_back("VisGenLeptonEta"); // Lepton Eta
  histograms.push_back("VisGenBJetpT"); // BJet pT
  histograms.push_back("VisGenBJetRapidity"); // BJet Rapidity
  //histograms.push_back("VisGenJetMulti"); // Jet Multiplicity
  histograms.push_back("VisGenExtraJetpT"); // Additional Jet pT
  histograms.push_back("VisGenTTBar1stJetMass"); // Invariant Mass TTBar+1Jet  

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

    sample_m166->Scale(1/sample_m166->Integral(),"width");
    sample_m169->Scale(1/sample_m169->Integral(),"width");
    sample_m171->Scale(1/sample_m171->Integral(),"width");
    sample_m172->Scale(1/sample_m172->Integral(),"width");
    sample_m173->Scale(1/sample_m173->Integral(),"width");
    sample_m175->Scale(1/sample_m175->Integral(),"width");
    sample_m178->Scale(1/sample_m178->Integral(),"width");


    bool rebin = true;
    
    if(rebin) {
      Int_t rebinning = 8;
      if((*hist) == "VisGenTTBar1stJetMass") { rebinning *= 4; }
      if((*hist) == "VisGenLeptonEta") { rebinning *= 2; }
      //if((*hist) == "VisGenTTBarRapidity") { rebinning /= 2; }
      if((*hist) == "VisGenTTBarpT") { rebinning *= 10; }
      if((*hist) == "VisGenTTBarMass") { rebinning *= 10; }

      sample_m166->Rebin(rebinning);
      sample_m169->Rebin(rebinning);
      sample_m171->Rebin(rebinning);
      sample_m172->Rebin(rebinning);
      sample_m173->Rebin(rebinning);
      sample_m175->Rebin(rebinning);
      sample_m178->Rebin(rebinning);
    }

    gStyle->SetOptStat(0);
    
    sample_m172->SetLineWidth(2);
    sample_m172->Draw();

    sample_m171->SetLineColor(kRed+3);
    sample_m169->SetLineColor(kRed);
    sample_m166->SetLineColor(kOrange);
    sample_m166->SetLineWidth(2);

    sample_m173->SetLineColor(kGreen+3);
    sample_m175->SetLineColor(kGreen);
    sample_m178->SetLineColor(kCyan);
    sample_m178->SetLineWidth(2);

    sample_m171->Draw("same");
    sample_m169->Draw("same");
    sample_m166->Draw("same");
    sample_m173->Draw("same");
    sample_m175->Draw("same");
    sample_m178->Draw("same");

    if(hist->Contains("Eta")) sample_m172->GetXaxis()->SetTitle("#eta");
    if(hist->Contains("pT")) sample_m172->GetXaxis()->SetTitle("p_{T} [GeV]");
    if(hist->Contains("Mass")) sample_m172->GetXaxis()->SetTitle("m_{tt#bar} [GeV]");
    if(hist->Contains("1stJetMass")) sample_m172->GetXaxis()->SetTitle("m_{t#bar{t}+Jet} [GeV]");
    if(hist->Contains("Rapidity")) sample_m172->GetXaxis()->SetTitle("Rapidity");

    TLegend *leg;
    if(hist->Contains("Rapidity")) { leg = new TLegend(0.15,0.5,0.3,.85); }
    else { leg = new TLegend(0.6,0.5,0.85,.85); }
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
