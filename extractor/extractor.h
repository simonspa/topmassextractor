#ifndef MASSEXTRACTOR_H
#define MASSEXTRACTOR_H

#include <TH1D.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include "log.h"

namespace massextractor {

  // Some ugly global defines:
#define nominalmass 172.5
#define granularity 500
#define confidenceLevel 0.95
#define chi2significance 1.5

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

  // In addition to the root file, also write PDFs (requires FLAG_STORE_HISTOGRAMS):
#define FLAG_STORE_PDFS 0x40

  // Either deliver minimum of the fit to the final ChiSquare or just the TGraph's minimum:
#define FLAG_RETURN_FITMIN 0x80

  // Do not take full covariance matrix into account, but just use stat. errors:
#define FLAG_DONT_USE_COVARIANCE 0x100

  // Do not subtract the background from data but compare all with background:
#define FLAG_DONT_SUBTRACT_BACKGROUND 0x200

  // Exclude theory predicition errors from bin-by-bin MC fit:
#define FLAG_NO_THEORYPREDICTION_ERRORS 0x400

  // Flag to explicitly exclude the statistical error on data in the chi2 calculation for systematic variation samples. If set, just the MC statistical errors are taken into account. This should only be used to evaluate the statistical errors of systematic variations, not to extract the systematic uncertainties.
#define FLAG_EXCLUDE_DATA_IN_VARIATION_STATERR 0x800

  class extractor {

  private:
    virtual TH1D * getSignalHistogram(Double_t mass, TFile * histos) = 0;
    virtual TH1D * getSimulationHistogram(Double_t mass, TFile * histos) = 0;
    virtual Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr) = 0;
    virtual Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr) = 0;

    virtual void getPredictionUncertainties() = 0;

    virtual TString getQuantity() = 0;
    virtual TString getRootFilename() = 0;
    virtual TFile * selectInputFile(TString sample) = 0;
    virtual TFile * selectInputFileTheory(TString channel, TString sample) = 0;

    virtual std::vector<Double_t> calcSampleDifference(TString nominal, TString systematic, TString histogram) = 0;

  protected:
    void getControlPlots(std::vector<TH1D*> histograms);

    // Functions for simple summed Chi2 extraction:
    Double_t chiSquare(const Double_t center, const Double_t center_widthsquared, const Double_t eval_widthsquared, const Double_t eval);
    std::pair<TGraphErrors*,TF1*> getChiSquare(std::vector<Double_t> masses, std::vector<TH1D*> data, std::vector<TH1D*> mc);

    // Functions for more involved fitted Chi2 extraction:
    std::vector<TGraphErrors*> splitBins(TString type, std::vector<Double_t> masses, std::vector<TH1D*> histograms);
    TGraphErrors * createIntersectionChiSquare(TGraphErrors* data, TGraphErrors* mc, Int_t bin, TGraphErrors* firstFit = NULL, TGraphErrors* secondFit = NULL);
    virtual std::pair<TGraphErrors*,TF1*> getFittedChiSquare(std::vector<Double_t> masses, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc);

    void getTrueConfidenceIntervals(Int_t n, Int_t ndim, const Double_t* x, Double_t* ci, Double_t cl = 0.95);
    void getTrueConfidenceIntervals(TGraphErrors * gr, Double_t cl = 0.95);

    // Minimization of the global Chi2 for extraction of the final mass value:
    Double_t getMinimum(std::pair<TGraphErrors*,TF1*> finalChiSquare);

    // Positive and negative statistical uncertainty on the extracted mass
    //calculated from varying the minimal Chi2 by +1 in either direction
    Double_t statErrorPos;
    Double_t statErrorNeg;
    // The extracted mass of the top quark in GeV
    Double_t extractedMass;

    // Path to parent folder of all input files
    TString m_inputpath;
    // Output path where histograms and PDFs as well as systematics tables will be stored
    TString m_outputpath;

    // The channel this extraction is performed on (ee, emu, mumu, combined)
    TString m_channel;
    // The sample we are using for this extraxtion
    TString m_sample;
    // Vector of all variations of the sample we are using, representing the 7 mass values
    std::vector<TString> samples;
    // Bin boundaries of the histograms under unvesigation. Calculated and filles in getSignalHistogram
    std::vector<Double_t> bin_boundaries;

    // Flag for separating runs in systematic variations from the actual mass extraction with nominal MC sample
    bool m_isSystematicVariation;

    // Storing the settings flags:
    uint32_t flags;

    // Do Closure test (i.e. replace signal sample with nominal mc)
    bool doClosure;
    TH1D * pseudoData;

    // Helper functions:
    void shiftGraph(TGraphErrors* ingraph, Double_t xshift, Double_t yshift);
    template<class t>
      bool isApprox(t a, t b, double eps = 0.01);
    /**
     * following numbers and mass dependence provided in NNLO paper arXiv:1303.6254
     * errors are NOT returned in % (so e.g. 0.026)
     */
    float getTtbarXsec(float topmass, float energy=8, float* scaleerr=0, float * pdferr=0);

    Double_t getMassFromSample(TString sample);
    TString getSampleFromMass(TString sample, Double_t mass, bool nominal);
    TFile * OpenFile(TString name, TString mode, bool output = false);

  public:
    // Return the extracted top mass - starts the extraction procedure.
    Double_t getTopMass();

    // Get both statistical error values: up and down, separately. The actual return value is the symmetrized statistical error
    virtual Double_t getStatError(Double_t &statPos, Double_t &statNeg);

    extractor(TString channel, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags);
    virtual ~extractor() {};
  };


  class extractorYield : public extractor {

  protected:
    Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);
    Double_t getPseudoData(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    Double_t m_stat_ndata;
    Double_t m_stat_nmc;

    TFile * selectInputFile(TString sample);
    TFile * selectInputFileTheory(TString /*channel*/, TString sample) { return selectInputFile(sample); }

    std::vector<Double_t> calcSampleDifference(TString nominal, TString systematic, TString histogram);

  private:
    TH1D * getSignalHistogram(Double_t mass, TFile * histos);
    TH1D * getSimulationHistogram(Double_t mass, TFile * histos);

    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_rec;
    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_bgr;
    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_ttbgr;
    void getPredictionUncertainties();

    inline TString getQuantity() { 
      if((flags &FLAG_NORMALIZE_YIELD) != 0) return "Events";
      else return "Events (norm.)";
    }
    inline TString getRootFilename() { return "MassFitRates.root"; }

    Double_t mcScalingFactor;

  public:
  extractorYield(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags) 
    : extractor(ch, sample, inputpath, outputpath, steeringFlags),
      m_prediction_errors_rec(),
      m_prediction_errors_bgr(),
      m_prediction_errors_ttbgr(),
      mcScalingFactor(0) {

      if((flags & FLAG_NORMALIZE_YIELD) != 0
	 && (flags & FLAG_LASTBIN_EXTRACTION) != 0) {
	LOG(unilog::logERROR) << "Normalization of a single bin doesn't make any sense. Dropping "
		      << "NORMALIZE_YIELD in favor for extracting from last bin only.";
	flags &= ~FLAG_NORMALIZE_YIELD;
      }
    };
    void setClosureSample(TString closure);
  };

  class extractorYieldOtherSamples : public extractorYield {

  private:
    Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    void calcDifferenceToNominal(TString nominal, TString systematics);

    std::vector<Double_t> deltaRec;
    std::vector<Double_t> deltaBgr;
    std::vector<Double_t> deltaTtbgr;

    TString m_systematic;

  public:
    Double_t getStatError(Double_t &statPos, Double_t &statNeg);

  extractorYieldOtherSamples(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString systematic) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), deltaRec(), deltaBgr(), deltaTtbgr(), m_systematic(systematic) {

      // This is a systematic variation run:
      m_isSystematicVariation = true;

      LOG(unilog::logDEBUG) << "Running for Match/Scale systematics: " << systematic;
      calcDifferenceToNominal(sample,systematic);
    };
  };

  class extractorYieldBackground : public extractorYield {

  private:
    Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareScaleFactor(TString systematic);
    Double_t scaleFactor;

  public:
  extractorYieldBackground(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString systematic) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), scaleFactor(1) {

      // This is a systematic variation run:
      m_isSystematicVariation = true;

      LOG(unilog::logDEBUG) << "Running for BG/DY systematics: " << systematic;
      prepareScaleFactor(systematic);
    };
  };

  class extractorDiffXSec : public extractor {

  private:
    TH1D * getSignalHistogram(Double_t mass, TFile * histos);
    TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
    TFile * selectInputFile(TString sample);
    TFile * selectInputFileTheory(TString channel, TString sample);

    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors;
    void getPredictionUncertainties();

    virtual Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco=0, Double_t bgr=0, Double_t ttbgr=0);
    Double_t getReco(Int_t /*bin*/, Double_t /*mass*/, Double_t /*reco*/, Double_t /*bgr*/, Double_t /*ttbgr*/) { return 0; };

    std::pair<TGraphErrors*,TF1*> getFittedChiSquare(std::vector<Double_t> masses, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc);
    // Function for fetching covariance matrix and inverting it:
    TMatrixD * getInverseCovMatrix(TString sample);

    std::vector<Double_t> calcSampleDifference(TString nominal, TString systematic, TString histogram);

    Double_t unfoldingMass;

    inline TString getQuantity() { return "#frac{1}{#sigma} #frac{d#sigma}{d#rho_{s}}"; }
    inline TString getRootFilename() { return "MassFitDiffXSec.root"; }

  public:
  extractorDiffXSec(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags) : extractor(ch, sample, inputpath, outputpath, steeringFlags), unfoldingMass(nominalmass) {
      LOG(unilog::logDEBUG) << "Extracting from Differential Cross Section.";

      // For HAD_UP, HAD_DOWN variations don't add the prediction uncertainties:
      if(sample.Contains("POWHEG") || sample.Contains("MCATNLO")) {
	flags |= FLAG_NO_THEORYPREDICTION_ERRORS;
	LOG(unilog::logDEBUG) << "Removed theory prediction errors for HAD variations: FLAG_NO_THEORYPREDICTION_ERRORS";
      }
    };
    void setUnfoldingMass(Double_t mass);
  };

  class extractorDiffXSecScaled : public extractorDiffXSec {

  private:
    Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareScaleFactors(TString systematic);
    std::vector<Double_t> scaleFactors;

  public:
  extractorDiffXSecScaled(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString scale) : extractorDiffXSec(ch, sample, inputpath, outputpath, steeringFlags), scaleFactors() {

      // This is a systematic variation run:
      m_isSystematicVariation = true;

      LOG(unilog::logDEBUG) << "Running sample " << sample << " with scale factors " << scale;
      prepareScaleFactors(scale);
    };
  };
}
#endif /* MASSEXTRACTOR_H */
