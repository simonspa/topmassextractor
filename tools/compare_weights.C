void compare_weights() {
  TFile pwhg("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m172p5.root");
  TFile pwhg_pos("POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m172p5_posweights.root");
  TFile output("comparison_pwhg_weights.root","recreate");

  gDirectory->pwd();

  std::vector<std::string> histograms;

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

  for(std::vector<std::string>::iterator hist = histograms.begin(); hist != histograms.end(); hist++) {
    
    std::cout << "Plotting " << (*hist) << "..." << std::endl;
    
    TCanvas *c = new TCanvas(hist->c_str(), hist->c_str());
    TH1D *sample_powhegbox = (TH1D*)pwhg.Get(hist->c_str());
    TH1D *sample_pwhg_pos = (TH1D*)pwhg_pos.Get(hist->c_str());
    sample_pwhg_pos->Scale(1/sample_pwhg_pos->Integral(),"width");
    sample_powhegbox->Scale(1/sample_powhegbox->Integral(),"width");

    bool rebin = true;
    
    if(rebin) {
      Int_t rebinning = 4;
      if((*hist) == "VisGenTTBar1stJetMass") { rebinning *= 4; }
      if((*hist) == "VisGenTTBarRapidity") { rebinning /= 2; }
      if((*hist) == "VisGenTTBarpT") { rebinning *= 10; }
      if((*hist) == "VisGenTTBarMass") { rebinning = 20; }

      sample_powhegbox->Rebin(rebinning);
      sample_pwhg_pos->Rebin(rebinning);
    }

    sample_powhegbox->Draw();
    sample_pwhg_pos->SetLineColor(kRed);
    sample_pwhg_pos->Draw("same");

    TLegend *leg = new TLegend(0.6,0.6,0.9,.9);
    leg->AddEntry(sample_powhegbox,"withnegweights 1","l");
    leg->AddEntry(sample_pwhg_pos,"withnegweights 0","l");
    leg->Draw();

    c->Draw();
    c->Write();
  }
}
