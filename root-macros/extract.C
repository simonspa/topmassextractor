#include <iostream>
#include <vector>
#include <string>
#include "Riostream.h"
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "log.h"
#include "extract.h"

Double_t nominalmass = 172.5;

using namespace unilog;


Double_t extractorMatchScale::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  // variation for every mass sample in the current bin:
  reco -= deltaRec.at(bin-1);
  bgr -= deltaBgr.at(bin-1);
  ttbgr -= deltaTtbgr.at(bin-1);

  // Call parent class signal calculation function:
  return extractor::getSignal(bin, mass, data, reco, bgr, ttbgr);
}

Double_t extractor::getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) {

  LOG(logDEBUG2) << "Calculation signal event count...";

  // Calculate the signal fraction from reconstructed events and TT background:
  Double_t fsignal = reco/(reco+ttbgr);

  // Correct the TTBgr for different TTBar Cross sections (mass dependent):
  Double_t corr_ttbgr = ttbgr*getTtbarXsec(mass)/getTtbarXsec(nominalmass);

  // Calculate signal by subtracting backround from data, multiplied by signal fraction.
  Double_t signal = (data - (bgr - corr_ttbgr))*fsignal;

  LOG(logDEBUG2) << "Bin #" << bin << ": data=" << data << " fsignal=" << fsignal << " sig=" << signal << " mc=" << reco;

  return signal;
}

Double_t extractor::getReco(Int_t bin, Double_t mass, Double_t reco) {

  LOG(logDEBUG2) << "Calculation reco event count...";

  Double_t corr_reco = reco*getTtbarXsec(mass)/getTtbarXsec(nominalmass);
  LOG(logDEBUG2) << "Bin #" << bin << ": reco=" << reco << " corr=" << corr_reco;

  // Return the reco event count corrected by the ttbar Xsec at given mass:
  return corr_reco;
}

Double_t extractorMatchScale::getReco(Int_t bin, Double_t mass, Double_t reco) {

  // Subtract the difference in event count for the nominal mass bin and systematics
  //variation for every bin:
  reco -= deltaRec.at(bin-1);

  // Call parent class signal calculation function:
  return extractor::getReco(bin, mass, reco);
}

TH1D * extractor::getSignalHistogram(Double_t mass, TFile * histos) {

  // Histogram containing data:
  TH1D * aDataHist = static_cast<TH1D*>(histos->Get("aDataHist"));
  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist"));
  // Histograms containing the background:
  TH1D * aTtBgrHist = static_cast<TH1D*>(histos->Get("aTtBgrHist"));
  TH1D * aBgrHist = static_cast<TH1D*>(histos->Get("aBgrHist"));

  Int_t nbins = aDataHist->GetNbinsX();
  LOG(logDEBUG) << "Data hist has " << nbins << " bins.";

  // Iterate over all bins:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    // Get signal corrected by ttbar background:
    Double_t signal = getSignal(bin, mass, aDataHist->GetBinContent(bin),
				aRecHist->GetBinContent(bin),
				aBgrHist->GetBinContent(bin),
				aTtBgrHist->GetBinContent(bin));
    
    // Write background subtrated signal:
    aDataHist->SetBinContent(bin,signal);
  }

  // Return signal-only histogram:
  return aDataHist;
}

TH1D * extractor::getSimulationHistogram(Double_t mass, TFile * histos) {

  // Histogram containing reconstructed events:
  TH1D * aRecHist = static_cast<TH1D*>(histos->Get("aRecHist"));

  Int_t nbins = aRecHist->GetNbinsX();
  LOG(logDEBUG) << "Reco hist has " << nbins << " bins.";

  // Iterate over all bins:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    // Correct the Reco events for different TTBar Cross sections (mass dependent):
    Double_t corr_reco = getReco(bin, mass,aRecHist->GetBinContent(bin));

    LOG(logDEBUG2) << "Bin #" << bin << ": reco=" << aRecHist->GetBinContent(bin) << " corr=" << corr_reco;
    
    // Write corrected Reco:
    aRecHist->SetBinContent(bin,corr_reco);
  }

  // Return reco histogram:
  return aRecHist;
}

std::vector<std::vector<Double_t> > extractor::splitBins(std::vector<TH1D*> histograms) {

  std::vector<std::vector<Double_t> > separated_bins;
  Int_t nbins = histograms.at(0)->GetNbinsX();

  // For every bin, prepare a new histogram containing all mass points:
  for(Int_t bin = 1; bin <= nbins; bin++) {

    LOG(logDEBUG) << "Filling bin " << bin << " mass points ";
 
    // Prepare new vector for this bin:
    std::vector<Double_t> thisbin;

    for(std::vector<TH1D*>::iterator hist = histograms.begin(); hist != histograms.end(); ++hist) {
      // Set the corresponding bin in output vector:
      thisbin.push_back((*hist)->GetBinContent(bin));
    }

    separated_bins.push_back(thisbin);
  }
  
  return separated_bins;
}

std::vector<TF1*> extractor::fitMassBins(TString channel, Int_t bin, std::vector<Double_t> masses, std::vector<Double_t> data, std::vector<Double_t> mc) {

  std::vector<TF1*> allfits;

  TString mname, dname, cname;
  mname.Form("mc_%i_",bin);
  dname.Form("dat_%i_",bin);
  cname.Form("dat_mc_%i_",bin);

  TCanvas* c;
  if(storeHistograms) {
    c = new TCanvas(cname+channel,cname+channel);
    c->cd();
  }

  TGraphErrors * graph_mc = new TGraphErrors();
  graph_mc->SetTitle(mname+channel);
  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph_mc->SetPoint(point, masses.at(point), mc.at(point));
    //graph_mc->SetPointError(point,0,sqrt(mc.at(point)));
  }

  graph_mc->SetLineWidth(0);
  graph_mc->SetMarkerStyle(20);
  graph_mc->Fit("pol2","Q");

  TF1 * mcfit = graph_mc->GetFunction("pol2");
  allfits.push_back(mcfit);

  if(storeHistograms) {
    graph_mc->Draw("A P E1");
    graph_mc->Write(mname+channel);
  }

  TGraphErrors * graph = new TGraphErrors();
  graph->SetTitle(dname+channel);
  for(UInt_t point = 0; point < masses.size(); ++point) {
    graph->SetPoint(point, masses.at(point), data.at(point));
    graph->SetPointError(point,0,sqrt(data.at(point)));
  }

  graph->SetLineWidth(0);
  graph->SetMarkerColor(kRed);
  graph->SetMarkerStyle(20);
  graph->Fit("pol2","Q");

  TF1 * datfit = graph->GetFunction("pol2");
  allfits.push_back(datfit);

  if(storeHistograms) {
    graph->Draw("same P E1");
    graph->Write(dname+channel);

    c->Write();
  }

  return allfits;
}

TF1 * extractor::getChiSquare(TString channel, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc) {

  TString name = "chi2_";

  TGraphErrors * chisquare = new TGraphErrors();
  chisquare->SetTitle(name+channel);

  for(UInt_t point = 0; point < masses.size(); ++point) {

    Double_t chi2 = 0;

    // Iterate over all bins:
    for(Int_t bin = 1; bin <= data.at(point)->GetNbinsX(); bin++) {
      Double_t chi = (data.at(point)->GetBinContent(bin) - mc.at(point)->GetBinContent(bin))/sqrt(data.at(point)->GetBinContent(bin));
      chi2 += chi*chi;
    }
    chisquare->SetPoint(point, masses.at(point), chi2);
  }

  chisquare->Fit("pol2","Q");
  if(storeHistograms) {
    chisquare->SetMarkerStyle(20);
    chisquare->Draw("AP");
    chisquare->Write(name+channel);
  }

  return chisquare->GetFunction("pol2");
}

Double_t extractor::getMinimum(TF1 * fit) {
  
  // For now, just return the function's minimum:
  return fit->GetMinimumX(0,330);
}

Double_t extractor::getTopMass() {

  std::vector<TH1D*> data_hists;
  std::vector<TH1D*> mc_hists;
  std::vector<Double_t> masses;

  for(std::vector<TString>::iterator sys = samples.begin(); sys != samples.end(); ++sys) {
    // Input files:
    TString filename = "preunfolded/" + (*sys) + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
    TFile * datafile = new TFile(filename);

    if(!datafile->IsOpen()) {
      LOG(logINFO) << "Failed to access data file " << filename;
      throw 1;
    }

    Double_t topmass = nominalmass;
    if(sys->Contains("UP")) {
      if(sys->Contains("1GEV")) topmass += 1;
      else if(sys->Contains("3GEV")) topmass += 3;
      else if(sys->Contains("6GEV")) topmass += 6;
    }
    else if(sys->Contains("DOWN")) {
      if(sys->Contains("1GEV")) topmass -= 1;
      else if(sys->Contains("3GEV")) topmass -= 3;
      else if(sys->Contains("6GEV")) topmass -= 6;
    }

    LOG(logDEBUG) << "Top Mass for Sample " << (*sys) << " m_t=" << topmass;
    
    // Subtract the estimated background from the data:
    TH1D * data = getSignalHistogram(topmass,datafile);
    TH1D * mc = getSimulationHistogram(topmass,datafile);

    masses.push_back(topmass);
    data_hists.push_back(data);
    mc_hists.push_back(mc);

  }

  std::vector<std::vector<Double_t> > separated_data = splitBins(data_hists);
  std::vector<std::vector<Double_t> > separated_mc = splitBins(mc_hists);

  TFile output;
  if(storeHistograms) {
    output.Open("massfit_bins.root","update");
    gDirectory->pwd();
  }

  std::vector<std::vector<TF1*> > fits;
  for(UInt_t bin = 0; bin < separated_data.size(); ++bin) {
    fits.push_back(fitMassBins(channel,bin+1,masses,separated_data.at(bin),separated_mc.at(bin)));
  }

  TF1 * fit = getChiSquare(channel,masses,data_hists,mc_hists);
  extracted_mass = getMinimum(fit);

  //getChiSquareFitted(channel,masses,data_hists,mc_hists);

  return extracted_mass;
}

void extractorMatchScale::calcDifferenceToNominal(TString nominal, TString systematics) {

  // Input files:
  TString nfilename = "preunfolded/" + nominal + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";
  TString sfilename = "preunfolded/" + systematics + "/" + channel + "/HypTTBar1stJetMass_UnfoldingHistos.root";

  TFile * nominalfile = new TFile(nfilename);
  TFile * systematicsfile = new TFile(sfilename);

  LOG(logDEBUG2) << nfilename;
  LOG(logDEBUG2) << sfilename;

  if(!nominalfile->IsOpen() || !systematicsfile->IsOpen()) {
    LOG(logINFO) << "Failed to access file";
    throw 1;
  }

  // Calculate (NOMINAL MASS - SYS_UP/DOWN) difference for every bin:
  TH1D * nominalReco = static_cast<TH1D*>(nominalfile->Get("aRecHist"));
  TH1D * varReco = static_cast<TH1D*>(systematicsfile->Get("aRecHist"));

  TH1D * nominalBgr = static_cast<TH1D*>(nominalfile->Get("aBgrHist"));
  TH1D * varBgr = static_cast<TH1D*>(systematicsfile->Get("aBgrHist"));

  TH1D * nominalTtbgr = static_cast<TH1D*>(nominalfile->Get("aTtBgrHist"));
  TH1D * varTtbgr = static_cast<TH1D*>(systematicsfile->Get("aTtBgrHist"));

  LOG(logDEBUG2) << nfilename;
  LOG(logDEBUG2) << sfilename;

  for(Int_t bin = 1; bin <= nominalReco->GetNbinsX(); bin++) {
    Double_t rec = nominalReco->GetBinContent(bin) - varReco->GetBinContent(bin);
    deltaRec.push_back(rec);
    
    Double_t bgr = nominalBgr->GetBinContent(bin) - varBgr->GetBinContent(bin);
    deltaBgr.push_back(bgr);

    Double_t ttbgr = nominalTtbgr->GetBinContent(bin) - varTtbgr->GetBinContent(bin);
    deltaTtbgr.push_back(ttbgr);

    LOG(logDEBUG2) << "Diffs bin " << bin << " drec=" << rec << " dbgr=" << bgr << " dttbgr=" << ttbgr;
  }

}


extractor::extractor(TString ch, std::vector<TString> samp, bool storeHistos) : channel(ch), samples(samp), storeHistograms(storeHistos) {
  LOG(logDEBUG) << "Initialized.";
}


void extract() {

  Log::ReportingLevel() = Log::FromString("INFO");

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");

  std::vector<TString> samples;
  samples.push_back("MASS_DOWN_6GEV");
  samples.push_back("MASS_DOWN_3GEV");
  samples.push_back("MASS_DOWN_1GEV");
  samples.push_back("Nominal");
  samples.push_back("MASS_UP_1GEV");
  samples.push_back("MASS_UP_3GEV");
  samples.push_back("MASS_UP_6GEV");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {
    extractor * mass_samples = new extractor(*ch,samples,true);
    Double_t topmass = mass_samples->getTopMass();
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass;
  }

  // Do the same (or similar things) for the systematics uncertainties:
  std::vector<TString> uncertainties;
  uncertainties.push_back("MATCH_UP");
  uncertainties.push_back("MATCH_DOWN");
  uncertainties.push_back("SCALE_UP");
  uncertainties.push_back("SCALE_DOWN");

  for(std::vector<TString>::iterator syst = uncertainties.begin(); syst != uncertainties.end(); ++syst) {

    LOG(logINFO) << "Getting " << (*syst) << " variation...";

    for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

      extractorMatchScale * extractor = new extractorMatchScale(*ch,samples,false,"Nominal",(*syst));

      Double_t topmass = extractor->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass;
    }

  }
}





template<class t>
bool extractor::isApprox(t a, t b, double eps) {
  if (fabs(a - b) < eps) { return true; }
  else { return false; }
}

float extractor::getTtbarXsec(float topmass, float energy, float* scaleerr, float * pdferr) {
    /*
     * all numbers following arxiv 1303.6254
     *
     */
    float mref=173.3;
    float referencexsec=0;
    float deltam=topmass-mref;


    float a1=0,a2=0;

    if(isApprox(energy,8.f,0.01)){
        a1=-1.1125;
        a2=0.070778;
        referencexsec=245.8;
        if(scaleerr)
            *scaleerr=0.034;
        if(pdferr)
            *pdferr=0.026;
    }
    else if(isApprox(energy,7.f,0.01)){
        a1=-1.24243;
        a2=0.890776;
        referencexsec=172.0;
        if(scaleerr)
            *scaleerr=0.034;
        if(pdferr)
            *pdferr=0.028;
    }

    float reldm=mref/(mref+deltam);

    float out= referencexsec* (reldm*reldm*reldm*reldm) * (1+ a1*(deltam)/mref + a2*(deltam/mref)*(deltam/mref));

    return out;
}

