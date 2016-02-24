#ifndef MASSEXTRACTOR_H
#define MASSEXTRACTOR_H

#include <TH1D.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TRandom.h>

#include "log.h"

namespace massextractor {

  // Some ugly global defines:
#define nominalmass 172.5
#define granularity 5000
#define confidenceLevel 0.68
#define chi2significance 1.5

#define PLOT_LOWER_LIMIT 168.5
#define PLOT_UPPER_LIMIT 177.0
#define PLOT_LOWER_LIMIT_DIFFXS 162.5
#define PLOT_UPPER_LIMIT_DIFFXS 176.0

  // Ship this flag in order to create and store histograms:
#define FLAG_STORE_HISTOGRAMS 0x01

  // Set this flag to normalize the distributions:
#define FLAG_NORMALIZE_DISTRIBUTIONS 0x02

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

  // Flag to explicitly exclude/ignore the statistical error on the MC sample in the chi2 calculation. If set, just the data statistical errors (convoluted with whatever has been used for unfolding) are taken into account. This should only be used to evaluate the statistical errors of systematic variations, not to extract the mass or any systematic uncertainties.
#define FLAG_IGNORE_MC_STATERR 0x800

  // Flag to explicitly set the statistical uncertainty on data (signal histogram) to zero (or something super tiny) to pretend having infinite statistics in data
#define FLAG_INFINITE_DATA_STATISTICS 0x1000

  // Do not plot channel labels:
#define FLAG_DONT_PLOT_CHANNELLABELS 0x2000

  // Use the Powheg NLO predictions instead of LO MadGraph:
#define FLAG_USE_NLO 0x4000


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
    // The sample we are using for this extraction
    TString m_sample;
    // Vector of all variations of the sample we are using, representing the 7 mass values
    std::vector<TString> samples;
    // Bin boundaries of the histograms under unvesigation. Calculated and filles in getSignalHistogram
    std::vector<Double_t> bin_boundaries;

    // Flag for separating runs in systematic variations from the actual mass extraction with nominal MC sample
    bool m_isSystematicVariation;

    // Flag for requesting calculation of prediction errors. Is set by constructor if FLAG_NO_THEORYPREDICTION_ERRORS is not present. Can also be set by derived classes that need prediction errors:
    bool m_requestPredictionErrors;

    // Storing the settings flags:
    uint32_t flags;

    // Do Closure test (i.e. replace signal sample with nominal mc)
    bool doClosure;
    TH1D * pseudoData;

    // File pointer for output histogram ROOT file:
    TFile * m_root_output;

    // Helper functions:
    void shiftGraph(TGraphErrors* ingraph, Double_t xshift, Double_t yshift);
    template<class t>
      bool isApprox(t a, t b, double eps = 0.01);
    /**
     * following numbers and mass dependence provided in NNLO paper arXiv:1303.6254
     * errors are NOT returned in % (so e.g. 0.026)
     */
    float getTtbarXsec(float topmass, float energy=8, float* scaleerr=0, float * pdferr=0);

    TString getSampleFromMass(TString sample, Double_t mass, bool nominal);
    TFile * OpenFile(TString name, TString mode);

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

    TH1D * getSignalHistogram(Double_t mass, TFile * histos);

  private:
    TH1D * getSimulationHistogram(Double_t mass, TFile * histos);

    void getPredictionUncertainties();

    inline TString getQuantity() { 
      if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) == 0) return "events";
      else return "events (norm.)";
    }
    inline TString getRootFilename() { return "MassFitRates.root"; }

  protected:
    Double_t mcScalingFactor;

    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_rec;
    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_bgr;
    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors_ttbgr;

  public:
  extractorYield(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags) 
    : extractor(ch, sample, inputpath, outputpath, steeringFlags),
      mcScalingFactor(0),
      m_prediction_errors_rec(),
      m_prediction_errors_bgr(),
      m_prediction_errors_ttbgr() {

      if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0
	 && (flags & FLAG_LASTBIN_EXTRACTION) != 0) {
	LOG(unilog::logERROR) << "Normalization of a single bin doesn't make any sense. Dropping "
			      << "NORMALIZE_DISTRIBUTIONS in favor for extracting from last bin only.";
	flags &= ~FLAG_NORMALIZE_DISTRIBUTIONS;
      }

      // Remove potential "ignore MC error" flag:
      if((flags & FLAG_IGNORE_MC_STATERR) != 0) {
	LOG(unilog::logWARNING) << "Removed FLAG_IGNORE_MC_STATERR, this should never be used for Yield extraction!";
	flags &= ~FLAG_IGNORE_MC_STATERR;
      }
    };
    void setClosureSample(TString closure);
  };

  class extractorYieldPseudoExp : public extractorYield {

  private:
    TH1D * getSignalHistogram(Double_t mass, TFile * histos);
    TRandom * myrnd;

  public:
  extractorYieldPseudoExp(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TRandom * random) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), myrnd(random) {

      LOG(unilog::logDEBUG) << "Running pseudo experiments on sample " << sample;
    };
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
      m_outputpath += "/" + m_systematic;

      LOG(unilog::logDEBUG) << "Running for Match/Scale systematics: " << systematic;
      calcDifferenceToNominal(sample,systematic);
    };
  };

  class extractorYieldPrediction : public extractorYield {

  private:
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);
    TString m_systematic;
  public:
  extractorYieldPrediction(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString systematic) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), m_systematic(systematic) {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      m_outputpath += "/" + m_systematic;

      // Request prediction errors to be calculated:
      m_requestPredictionErrors = true;
      LOG(unilog::logDEBUG) << "Running for Prediction systematics: " << systematic;
    };
  };

  class extractorYieldBackground : public extractorYield {

  private:
    Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco, Double_t bgr, Double_t ttbgr);
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareScaleFactor(TString systematic);
    Double_t scaleFactor;

  public:
  extractorYieldBackground(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString systematic) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), scaleFactor(1) {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      m_outputpath += "/" + systematic;

      LOG(unilog::logDEBUG) << "Running for BG/DY systematics: " << systematic;
      prepareScaleFactor(systematic);
    };
  };

  class extractorYieldScaled : public extractorYield {

  private:
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareScaleFactors(TString systematic);
    std::vector<Double_t> getPDFScaleFactors(Int_t sign, TString channel);
    std::vector<Double_t> scaleFactors;

  public:
  extractorYieldScaled(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString scale) : extractorYield(ch, sample, inputpath, outputpath, steeringFlags), scaleFactors() {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      m_outputpath += "/" + scale;

      LOG(unilog::logDEBUG) << "Running sample " << sample << " with scale factors " << scale;
      prepareScaleFactors(scale);
    };
  };

  class extractorDiffXSec : public extractor {

  private:
    TH1D * getSimulationHistogram(Double_t mass, TFile * histos);
    TFile * selectInputFile(TString sample);

    std::vector<std::pair<Double_t,Double_t> > m_prediction_errors;
    void getPredictionUncertainties();

    virtual Double_t getSignal(Int_t bin, Double_t mass, Double_t data, Double_t reco=0, Double_t bgr=0, Double_t ttbgr=0);
    virtual Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr=0, Double_t ttbgr=0);

    std::pair<TGraphErrors*,TF1*> getFittedChiSquare(std::vector<Double_t> masses, std::vector<TGraphErrors*> data, std::vector<TGraphErrors*> mc);
    // Helper function to read matrix from file and divide bins by their width:
    TMatrixD * readMatrix(TString sample, TString channel);
    // Function for fetching covariance matrix and preparing it:
    TMatrixD * getCovMatrix(TString sample);
    TMatrixD * invertCovMatrix(TMatrixD * cov, std::vector<Double_t> mcstats, Int_t drop_bin);

    // Mass to be used for unfolding procedure:
    Double_t unfoldingMass;
    // Total number of events in the input signal histogram, to be used for COV matrix normalization:
    Int_t m_nevents;

    inline TString getQuantity() { 
      // Non-normalized diff. cross section:
      if((flags & FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { return "#frac{1}{#sigma} #frac{d#sigma}{d#rho_{s}}"; }
      // Normalized diff. cross section:
      else { return "#frac{d#sigma}{d#rho_{s}}"; }
    }
    inline TString getRootFilename() { return "MassFitDiffXSec.root"; }

  protected:
    virtual TH1D * getSignalHistogram(Double_t mass, TFile * histos);
    TFile * selectInputFileTheory(TString channel, TString sample);
    std::vector<Double_t> calcSampleDifference(TString nominal, TString systematic, TString histogram);

  public:
  extractorDiffXSec(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags) : extractor(ch, sample, inputpath, outputpath, steeringFlags), unfoldingMass(nominalmass) {
      LOG(unilog::logDEBUG) << "Extracting from Differential Cross Section.";

      // For HAD_UP, HAD_DOWN variations don't add the prediction uncertainties:
      if(sample.Contains("POWHEG") || sample.Contains("MCATNLO")) {
	flags |= FLAG_NO_THEORYPREDICTION_ERRORS;
	LOG(unilog::logDEBUG) << "Removed theory prediction errors for HAD variations: FLAG_NO_THEORYPREDICTION_ERRORS";
      }

      if((flags&FLAG_IGNORE_MC_STATERR) != 0) { LOG(unilog::logWARNING) << "This run will ignore all statistical (and theory prediction) errors assigned to the MC sample!"; }
    };
    void setUnfoldingMass(Double_t mass);
  };

  class extractorDiffXSecPseudoExp : public extractorDiffXSec {

  private:
    TH1D * getSignalHistogram(Double_t mass, TFile * histos);
    TRandom * myrnd;

  public:
  extractorDiffXSecPseudoExp(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TRandom * random) : extractorDiffXSec(ch, sample, inputpath, outputpath, steeringFlags), myrnd(random) {

      LOG(unilog::logDEBUG) << "Running pseudo experiments on sample " << sample;
    };
  };

  class extractorDiffXSecScaled : public extractorDiffXSec {

  private:
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareScaleFactors(TString systematic);
    std::vector<Double_t> getPDFScaleFactors(Int_t sign, TString channel);
    std::vector<Double_t> scaleFactors;

  public:
  extractorDiffXSecScaled(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString scale) : extractorDiffXSec(ch, sample, inputpath, outputpath, steeringFlags), scaleFactors() {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      m_outputpath += "/" + scale;

      LOG(unilog::logDEBUG) << "Running sample " << sample << " with scale factors " << scale;
      prepareScaleFactors(scale);
    };
  };

  class extractorDiffXSecPrediction : public extractorDiffXSec {

  private:
    Double_t getReco(Int_t bin, Double_t mass, Double_t reco, Double_t bgr, Double_t ttbgr);

    void prepareShiftFactors(TString systematic);
    std::vector<Double_t> m_shiftFactors;

  public:
  extractorDiffXSecPrediction(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString scale) : extractorDiffXSec(ch, sample, inputpath, outputpath, steeringFlags), m_shiftFactors() {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      if(sample == "Nominal") m_outputpath += "/" + scale;

      LOG(unilog::logDEBUG) << "Running Theory Prediction error " << sample << ", shifting MC using " << scale;
      prepareShiftFactors(scale);
    };
  };

  class extractorDiffXSecGenLevelPrediction : public extractorDiffXSec {

  private:
    TFile * selectInputFileTheory(TString channel, TString sample);
    TString m_systematic;

  public:
  extractorDiffXSecGenLevelPrediction(TString ch, TString sample, TString inputpath, TString outputpath, uint32_t steeringFlags, TString systematic) : extractorDiffXSec(ch, sample, inputpath, outputpath, steeringFlags) {

      // This is a systematic variation run:
      m_isSystematicVariation = true;
      m_systematic = systematic;
      if(sample == "Nominal") m_outputpath += "/" + systematic;

      LOG(unilog::logDEBUG) << "Running Theory Prediction uncertainty " << systematic << " using gen-level histograms, data from sample " << sample;
    };
  };

}
#endif /* MASSEXTRACTOR_H */
