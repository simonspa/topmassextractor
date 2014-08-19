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

// Ship this flag in order to create and store histograms:
#define FLAG_STORE_HISTOGRAMS 0x01

// Set this flag to normalize the yield distributions:
#define FLAG_NORMALIZE_YIELD 0x02

// Allow restriction to last bin of distributions:
#define FLAG_LASTBIN_EXTRACTION 0x04

// Use data DiffXSec distributions unfolded with different mass samples:
#define FLAG_UNFOLD_ALLMASSES 0x08

// Get Chi2 distribution from already fitted bin distributions:
#define FLAG_CHISQUARE_FROM_FITS 0x10

// Do not shift the data & MC graphs to (0,0) before fitting:
#define FLAG_DONT_SHIFT_GRAPHS 0x20


class extractor {

private:
  virtual TH1D * getSignalHistogram(Double_t mass, TFile * histos);
  virtual TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
  
  virtual Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
  virtual Double_t getReco(Int_t bin, Double_t mass, Double_t reco);

  void getControlPlots(std::vector<TH1D*> histograms);

  std::vector<TH1D* > splitBins(TString type, std::vector<TH1D*> histograms);
  std::pair<TGraphErrors*,TGraphErrors*> fitMassBins(TString channel, Int_t bin, std::vector<Double_t> masses, TH1D* data, TH1D* mc);

  TGraphErrors * getShiftedGraph(TGraphErrors* ingraph, Double_t xshift, Double_t yshift);
  Double_t chiSquare(const Double_t center, const Double_t widthsquared, const Double_t eval);
  TF1 * getChiSquare(TString channel, std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc);

  TF1 * getFittedChiSquare(TString channel, std::vector<Double_t> masses, std::vector<std::pair<TGraphErrors*,TGraphErrors*> > fits);
  TGraphErrors * createIntersectionChiSquare(std::pair<TGraphErrors*,TGraphErrors*> fits);

  Double_t getMinimum(TF1 * fit);

  Double_t statError;
  Double_t extractedMass;

  TString channel;
  std::vector<TString> samples;

  // Stroing the settings flags:
  uint32_t flags;

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

 protected:
  virtual inline TString getQuantity() { return "Events"; }
  virtual TFile * selectInputFile(TString sample, TString channel);

 public:
  Double_t getTopMass();
  Double_t getStatError();
  TString getSampleLabel(TString systematic);
  virtual void setClosureSample(TString closure);
  extractor(TString channel, TString sample, uint32_t steeringFlags);
};


class extractorOtherSamples : public extractor {

private:
  Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
  Double_t getReco(Int_t bin, Double_t mass, Double_t reco);

  void calcDifferenceToNominal(TString nominal, TString systematics);

  std::vector<Double_t> deltaRec;
  std::vector<Double_t> deltaBgr;
  std::vector<Double_t> deltaTtbgr;

public:
 extractorOtherSamples(TString ch, TString sample, uint32_t steeringFlags, TString systematic) : extractor(ch, sample, steeringFlags), deltaRec(), deltaBgr(), deltaTtbgr() {
    LOG(unilog::logDEBUG) << "Running for Match/Scale systematics: " << systematic;
    calcDifferenceToNominal(sample,systematic);
  };
};

class extractorBackground : public extractor {

private:
  Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

  void prepareScaleFactor(TString systematic);
  Double_t scaleFactor;

public:
 extractorBackground(TString ch, TString sample, uint32_t steeringFlags, TString systematic) : extractor(ch, sample, steeringFlags), scaleFactor(1) {
    LOG(unilog::logDEBUG) << "Running for BG/DY systematics: " << systematic;
    prepareScaleFactor(systematic);
  };
};

class extractorDiffXSec : public extractor {

 private:
  TH1D * getSignalHistogram(Double_t mass, TFile * histos);
  TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
  TFile * selectInputFile(TString sample, TString channel);

 protected:
  inline TString getQuantity() { return "#frac{1}{#sigma} #frac{d#sigma}{d#rho_{s}} #left[GeV^{-1}#right]"; }

public:
 extractorDiffXSec(TString ch, TString sample, uint32_t steeringFlags) : extractor(ch, sample, steeringFlags) {
    LOG(unilog::logDEBUG) << "Extracting from Differential Cross Section.";
  };
  void setClosureSample(TString closure);

};
#endif /* EXTRACT_H */
