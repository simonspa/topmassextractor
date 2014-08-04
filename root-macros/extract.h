#include <vector>
#include <string>
#include "Riostream.h"
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#ifndef EXTRACT_H
#define EXTRACT_H

class extractTopMass {

private:
  virtual TH1D * getSignalHistogram(Double_t mass, TFile * histos);
  virtual TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
  
  virtual Double_t getSignal(Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

  virtual Double_t getReco(Double_t mass, Double_t reco);

  std::vector<std::vector<Double_t> > splitBins(std::vector<TH1D*> histograms);
  std::vector<TF1*> fitMassBins(TString channel, Int_t bin, std::vector<Double_t> masses, std::vector<Double_t> data, std::vector<Double_t> mc);
  TF1 * getChiSquare(TString channel, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc);
  Double_t getMinimum(TF1 * fit);

  Double_t extracted_mass;
  TString channel;
  std::vector<TString> samples;

  bool storeHistograms;
public:
  Double_t getTopMass();
  extractTopMass(TString channel, std::vector<TString> systematics, bool storeHistos);
};


class extractTopMassMatchScale : public extractTopMass {

private:
  Double_t getSignal(Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
  Double_t getReco(Double_t mass, Double_t reco);

  Double_t deltaNevents;
public:
  extractTopMassMatchScale(TString channel, std::vector<TString> systematics, bool storeHistos, Double_t deltaN) : deltaNevents(deltaN), extractTopMass(channel, systematics, storeHistos) {
    LOG(unilog::logINFO) << "Running for Match/Scale systematics.";
  };

};

#endif /* EXTRACT_H */
