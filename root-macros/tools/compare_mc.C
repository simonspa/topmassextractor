void compare_mc() {

  TString base = "selectionRoot/";

  TString filenominal = "Nominal/ee/ee_ttbarsignalplustau.root";
  TString filenominal2 = "POWHEG/ee/ee_ttbarsignalplustau_powheg.root";

  TString filepowhegbox = "POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m172p5.root";
  TString filepowhegbox2 = "POWHEG/ee/ee_ttbarsignalplustau_pwhgbox_m172p5-secondrun.root";

  TFile pwhg(base+filepowhegbox);
  TFile pwhg_pos(base+filepowhegbox2);
  TFile nominal(base+filenominal);
  TFile nominal_pwhg(base+filenominal2);

  TFile output("comparison_pwhg_madg.root","recreate");

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

    TCanvas *c = new TCanvas(*hist, *hist);
    TH1D *sample_powhegbox = (TH1D*)pwhg.Get(*hist);
    TH1D *sample_powhegbox_pos = (TH1D*)pwhg_pos.Get(*hist);
    TH1D *sample_nominal = (TH1D*)nominal.Get(*hist);
    TH1D *sample_nominal_pwhg = (TH1D*)nominal_pwhg.Get(*hist);

    sample_nominal->Scale(1/sample_nominal->Integral(),"width");
    sample_nominal_pwhg->Scale(1/sample_nominal_pwhg->Integral(),"width");
    sample_powhegbox->Scale(1/sample_powhegbox->Integral(),"width");
    sample_powhegbox_pos->Scale(1/sample_powhegbox_pos->Integral(),"width");

    bool rebin = true;
  

    if(rebin) {
      Int_t rebinning = 8;
      if((*hist) == "VisGenTTBar1stJetMass") { rebinning *= 4; }
      if((*hist) == "VisGenLeptonEta") { rebinning *= 2; }
      //if((*hist) == "VisGenTTBarRapidity") { rebinning /= 2; }
      if((*hist) == "VisGenTTBarpT") { rebinning *= 10; }
      if((*hist) == "VisGenTTBarMass") { rebinning *= 10; }

      sample_powhegbox->Rebin(rebinning);
      sample_powhegbox_pos->Rebin(rebinning);
      sample_nominal->Rebin(rebinning);
      sample_nominal_pwhg->Rebin(rebinning);
    }

    gStyle->SetOptStat(0);

    sample_nominal->SetLineWidth(2);
    sample_nominal->Draw();

    sample_nominal_pwhg->SetLineColor(kGreen+2);
    sample_nominal_pwhg->SetLineWidth(2);
    sample_nominal_pwhg->Draw("same");

    sample_powhegbox->SetLineColor(kRed);
    sample_powhegbox->SetLineWidth(2);
    sample_powhegbox->Draw("same");

    sample_powhegbox_pos->SetLineColor(kBlack);
    sample_powhegbox_pos->SetLineWidth(2);
    sample_powhegbox_pos->Draw("same");

    if(hist->Contains("Eta")) sample_nominal->GetXaxis()->SetTitle("#eta");
    if(hist->Contains("pT")) sample_nominal->GetXaxis()->SetTitle("p_{T} [GeV]");
    if(hist->Contains("Mass")) sample_nominal->GetXaxis()->SetTitle("m_{tt#bar} [GeV]");
    if(hist->Contains("1stJetMass")) sample_nominal->GetXaxis()->SetTitle("m_{t#bar{t}+Jet} [GeV]");
    if(hist->Contains("Rapidity")) sample_nominal->GetXaxis()->SetTitle("Rapidity");

    TLegend *leg;
    if(hist->Contains("Rapidity")) { leg = new TLegend(0.15,0.6,0.4,.85); }
    else { leg = new TLegend(0.5,0.6,0.85,.85); }
    leg->AddEntry(sample_powhegbox,"POWHEG ttbar+1jet NLO CT10","l");
    leg->AddEntry(sample_powhegbox_pos,"POWHEG ttbar+1jet NLO CTEQ66","l");

    leg->AddEntry(sample_nominal,"Nominal: MadGraph","l");
    leg->AddEntry(sample_nominal_pwhg,"Nominal: POWHEG","l");

    leg->Draw();

    c->Draw();
    c->Write();
    //c->Print((*hist) + ".eps");
    //c->Print((*hist) + ".pdf");
  }
}
