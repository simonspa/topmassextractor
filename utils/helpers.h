#ifndef MASSEXTRACTOR_HELPERS_H
#define MASSEXTRACTOR_HELPERS_H

#include <sstream>
#include "extractor.h"

namespace massextractor {
/** Helper function to return a printed list of flags
 */
  std::string inline listFlags(uint32_t flags) {
    std::stringstream os;

    // No flags given:
    if(flags == 0) return "(0) ";

    if((flags&FLAG_STORE_HISTOGRAMS) != 0) { os << "FLAG_STORE_HISTOGRAMS, "; flags -= FLAG_STORE_HISTOGRAMS; }
    if((flags&FLAG_NORMALIZE_DISTRIBUTIONS) != 0) { os << "FLAG_NORMALIZE_DISTRIBUTIONS, "; flags -= FLAG_NORMALIZE_DISTRIBUTIONS; }
    if((flags&FLAG_LASTBIN_EXTRACTION) != 0) { os << "FLAG_LASTBIN_EXTRACTION, "; flags -= FLAG_LASTBIN_EXTRACTION; }
    if((flags&FLAG_UNFOLD_ALLMASSES) != 0) { os << "FLAG_UNFOLD_ALLMASSES, "; flags -= FLAG_UNFOLD_ALLMASSES; }
    if((flags&FLAG_CHISQUARE_FROM_FITS) != 0) { os << "FLAG_CHISQUARE_FROM_FITS, "; flags -= FLAG_CHISQUARE_FROM_FITS; }
    if((flags&FLAG_DONT_SHIFT_GRAPHS) != 0) { os << "FLAG_DONT_SHIFT_GRAPHS, "; flags -= FLAG_DONT_SHIFT_GRAPHS; }
    if((flags&FLAG_STORE_PDFS) != 0) { os << "FLAG_STORE_PDFS, "; flags -= FLAG_STORE_PDFS; }
    if((flags&FLAG_RETURN_FITMIN) != 0) { os << "FLAG_RETURN_FITMIN, "; flags -= FLAG_RETURN_FITMIN; }
    if((flags&FLAG_DONT_USE_COVARIANCE) != 0) { os << "FLAG_DONT_USE_COVARIANCE, "; flags -= FLAG_DONT_USE_COVARIANCE; }
    if((flags&FLAG_DONT_SUBTRACT_BACKGROUND) != 0) { os << "FLAG_DONT_SUBTRACT_BACKGROUND, "; flags -= FLAG_DONT_SUBTRACT_BACKGROUND; }
    if((flags&FLAG_NO_THEORYPREDICTION_ERRORS) != 0) { os << "FLAG_NO_THEORYPREDICTION_ERRORS, "; flags -= FLAG_NO_THEORYPREDICTION_ERRORS; }
    if((flags&FLAG_IGNORE_MC_STATERR) != 0) { os << "FLAG_IGNORE_MC_STATERR, "; flags -= FLAG_IGNORE_MC_STATERR; }
    if((flags&FLAG_INFINITE_DATA_STATISTICS) != 0) { os << "FLAG_INFINITE_DATA_STATISTICS, "; flags -= FLAG_INFINITE_DATA_STATISTICS; }
    if((flags&FLAG_DONT_PLOT_CHANNELLABELS) != 0) { os << "FLAG_DONT_PLOT_CHANNELLABELS, "; flags -= FLAG_DONT_PLOT_CHANNELLABELS; }
    if((flags&FLAG_USE_NLO) != 0) { os << "FLAG_USE_NLO, "; flags -= FLAG_USE_NLO; }
    if(flags != 0) os << "Unknown flag: " << flags;
    return os.str();
  }

  std::vector<Double_t> getBinningFromHistogram(TH1D * histo, Int_t startbin=1, Int_t nbins=0);

  // Splitting a string at tokens:
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);

  // Helper to calculate the stat. error on systematics for DiffXSec:
  Double_t systStatErr(Double_t nominalStatErr, Double_t systStatErr);

  TString getSampleLabel(TString systematic);
  Double_t getMassFromSample(TString sample);

  // Helper to check if both systematic variations (UP and DOWN) produce errors in the same
  // direction. If so, only add the maximum of both, otherwise count up and down.
  // Return value is decision (true: same direction, false: different directions)
  bool getSystematicUpDownError(Double_t delta_up, Double_t delta_down, Double_t & total_syst_pos, Double_t & total_syst_neg);

}

#endif /* MASSEXTRACTOR_HELPERS_H */
