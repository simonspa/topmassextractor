#ifndef EXTRACT_H
#define EXTRACT_H

#include <vector>
#include <string>

#include <TMath.h>
#include <TSystem.h>

namespace bin {

  void extract_yield(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake);
  void extract_yield_stats(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake);

  void extract_diffxsec(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake);
  void extract_diffxsec_stats(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake);

  // Splitting a string at tokens:
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);

  // Helper to calculate the stat. error on systematics for DiffXSec:
  Double_t systStatErr(Double_t nominalStatErr, Double_t systStatErr);

  TString getSampleLabel(TString systematic);

  // Helper to check if both systematic variations (UP and DOWN) produce errors in the same
  // direction. If so, only add the maximum of both, otherwise count up and down.
  // Return value is decision (true: same direction, false: different directions)
  bool getSystematicUpDownError(Double_t delta_up, Double_t delta_down, Double_t & total_syst_pos, Double_t & total_syst_neg);

}

#endif /*EXTRACT_H*/
