#include <vector>
#include <string>
#include "Riostream.h"
#include <TH1D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#ifndef EXTRACT_H
#define EXTRACT_H

class extractor {

private:
  virtual TH1D * getSignalHistogram(Double_t mass, TFile * histos);
  virtual TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
  
  virtual Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

  virtual Double_t getReco(Int_t bin, Double_t mass, Double_t reco);

  std::vector<std::vector<Double_t> > splitBins(std::vector<TH1D*> histograms);
  std::vector<TF1*> fitMassBins(TString channel, Int_t bin, std::vector<Double_t> masses, std::vector<Double_t> data, std::vector<Double_t> mc);
  TF1 * getChiSquare(TString channel, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc);
  Double_t getMinimum(TF1 * fit);

  Double_t extractedMass;
  Double_t statError;

  TString channel;
  std::vector<TString> samples;

  // Whether to produce and write out histograms or not
  bool storeHistograms;

  // Do Closure test (i.e. replace signal sample with nominal mc)
  bool doClosure;
  TH1D * pseudoData;

  // Helper functions:
  template<class t>
    bool isApprox(t a, t b, double eps = 0.01);
  /**
   * following numbers and mass dependence provided in NNLO paper arXiv:1303.6254
   * errors are NOT returned in % (so e.g. 0.026)
   */
  float getTtbarXsec(float topmass, float energy=8, float* scaleerr=0, float * pdferr=0);
  void setHHStyle(TStyle& HHStyle);
  void DrawDecayChLabel(TString decaychannel="", double textSize=0.04);
  void DrawCMSLabels(int cmsprelim=true, double energy=8, double textSize=0.04);
  void setStyle(TGraphErrors *hist);
  void setStyleAndFillLegend(TGraphErrors* histo, TString name, TLegend *leg);
  void setLegendStyle(TLegend *leg);

  Double_t getMassFromSample(TString sample);
  TString getChannelLabel(TString channel);

public:
  Double_t getTopMass();
  Double_t getStatError();
  void setClosureSample(TString closure);
  extractor(TString channel, TString sample, bool storeHistos);
};


class extractorMatchScale : public extractor {

private:
  Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
  Double_t getReco(Int_t bin, Double_t mass, Double_t reco);

  void calcDifferenceToNominal(TString nominal, TString systematics);

  std::vector<Double_t> deltaRec;
  std::vector<Double_t> deltaBgr;
  std::vector<Double_t> deltaTtbgr;

public:
 extractorMatchScale(TString channel, TString sample, bool storeHistos, TString nominal, TString systematic) : deltaRec(), deltaBgr(), deltaTtbgr(), extractor(channel, sample, storeHistos) {
    LOG(unilog::logDEBUG) << "Running for Match/Scale systematics: " << systematic;
    calcDifferenceToNominal(nominal,systematic);
  };

};

class extractorBackground : public extractor {

private:
  Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

  void prepareScaleFactor(TString systematic);
  Double_t scaleFactor;

public:
 extractorBackground(TString channel, TString sample, bool storeHistos, TString systematic) : scaleFactor(1), extractor(channel, sample, storeHistos) {
    LOG(unilog::logDEBUG) << "Running for BG/DY systematics: " << systematic;
    prepareScaleFactor(systematic);
  };
};

#endif /* EXTRACT_H */
