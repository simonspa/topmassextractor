#ifndef EXTRACT_H
#define EXTRACT_H

#include <vector>
#include <string>

#include <TMath.h>
#include <TSystem.h>

namespace massextractor {

  void extract_yield(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake);
  void extract_yield_stats(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake);

  void extract_diffxsec(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake);
  void extract_diffxsec_stats(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake);

}

#endif /*EXTRACT_H*/
