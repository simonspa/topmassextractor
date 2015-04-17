#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdint.h>

#include "extract.h"
#include "extractor.h"
#include "log.h"
#include "latex.h"
#include "helpers.h"

#include <TMath.h>
#include <TSystem.h>

using namespace std;
using namespace massextractor;
using namespace unilog;

int main(int argc, char* argv[]) {

  Log::ReportingLevel() = Log::FromString("INFO");
  FILE* logfile;

  TString outputpath = "ExtractionResults";
  TString inputpath = "";
  std::string type, channel;

  std::vector<std::string> flagtokens;
  std::vector<std::string> systlist;
  uint32_t flags = 0;
  bool syst = false;

  // If this bool is set, PDFs for *all* systematics are produced (many!)
  bool fulltake = false;

  bool closure = false;
  TString closure_sample = "Nominal";

  for (int i = 1; i < argc; i++) {
    // select to either extract from yield of differential cross section:
    if (!strcmp(argv[i],"-t")) {
      type = string(argv[++i]);
      if(type !=  "diffxs" && type != "yield" && type != "yieldstats" && type != "diffxsstats") {
	LOG(logERROR) << "Unknown extraction type \"" << type << "\".";
	return 0;
      }
    }
    // Set the verbosity level:
    else if (!strcmp(argv[i],"-v")) { Log::ReportingLevel() = Log::FromString(argv[++i]); }
    // Set output path option
    else if(!strcmp(argv[i],"-i")) { inputpath = string(argv[++i]); }
    // Set output path option
    else if(!strcmp(argv[i],"-o")) { outputpath = string(argv[++i]); }
    // Set "data" flag to run on real data instead of closure test (just yield):
    else if(!strcmp(argv[i],"-d")) { closure = false; }
    // Set "pseudo" flag to run on pseudo data (just yield):
    else if(!strcmp(argv[i],"--pseudo")) { closure = true; }
    // Set systematics flag to run on systematic variation sample
    else if(!strcmp(argv[i],"-s")) { syst = true; systlist = split(string(argv[++i]), ','); }
    // Select the channel to run on:
    else if(!strcmp(argv[i],"-c")) { channel = string(argv[++i]); }
    // Allow additional logging to file:
    else if(!strcmp(argv[i],"-l")) {
      logfile = fopen(argv[++i], "a");
      unilog::SetLogOutput::Stream() = logfile;
      unilog::SetLogOutput::Duplicate() = true;
    }
    // Mass sample to be used for closure test:
    else if(!strcmp(argv[i],"-m")) { closure_sample = string(argv[++i]); }
    // Read and tokeinze the flags:
    else if(!strcmp(argv[i],"-f")) { flagtokens = split(string(argv[++i]), ','); }
    else { LOG(logERROR) << "Unrecognized command line argument \"" << argv[i] << "\".";}
  }

  // Check that the channel we selected is fine:
  std::vector<TString> channels;
  if(channel == "ee" || channel.empty()) { channels.push_back("ee"); }
  if(channel == "emu" || channel.empty()) { channels.push_back("emu"); }
  if(channel == "mumu" || channel.empty()) { channels.push_back("mumu"); }
  if(channel == "combined" || channel.empty()) { channels.push_back("combined");}
  if(channels.empty()) {
    LOG(logERROR) << "No valid channel selected.";
    return 0;
  }
  else {
    std::stringstream ch;
    for(std::vector<TString>::iterator c = channels.begin(); c != channels.end(); ++c) { ch << *c << ", "; }
    LOG(logINFO) << "Running on channels " << ch.str();
  }

  if(syst) {
    std::stringstream ch;
    for(std::vector<std::string>::iterator c = systlist.begin(); c != systlist.end(); ++c) { ch << *c << ", "; }
    LOG(logINFO) << "Running on systematics " << ch.str();
  }

  if(gSystem->AccessPathName(inputpath)) {
    LOG(logERROR) << "Input path \"" << inputpath << "\" does not exist!";
    return 0;
  }
  // Check and create output path if necessary:
  if(gSystem->AccessPathName(outputpath)) {
    gSystem->mkdir(outputpath, true);
    LOG(logDEBUG) << "Created output directory \"" << outputpath << "\"";
  }

  // Standard flag settings: Chi2 fit, no Covariance, normalized.
  flags = FLAG_CHISQUARE_FROM_FITS | FLAG_NORMALIZE_DISTRIBUTIONS | FLAG_DONT_SUBTRACT_BACKGROUND | FLAG_NO_THEORYPREDICTION_ERRORS;
  // Check and assign the flags:
  for(std::vector<std::string>::iterator tok = flagtokens.begin(); tok != flagtokens.end(); ++tok) {
    if(*tok == "fit") { flags |= FLAG_CHISQUARE_FROM_FITS; }
    if(*tok == "nofit") { flags &= ~FLAG_CHISQUARE_FROM_FITS; }
    else if(*tok == "pdf") { flags |= FLAG_STORE_PDFS; }
    else if(*tok == "pdfall") { flags |= FLAG_STORE_PDFS; fulltake = true; }
    else if(*tok == "root") { flags |= FLAG_STORE_HISTOGRAMS; }
    else if(*tok == "cov") { flags &= ~FLAG_DONT_USE_COVARIANCE; }
    else if(*tok == "nocov") { flags |= FLAG_DONT_USE_COVARIANCE; }
    else if(*tok == "norm") { flags |= FLAG_NORMALIZE_DISTRIBUTIONS; }
    else if(*tok == "nonorm") { flags &= ~FLAG_NORMALIZE_DISTRIBUTIONS; }
    else if(*tok == "lastbin") { flags |= FLAG_LASTBIN_EXTRACTION; }
    else if(*tok == "bgr") { flags |= FLAG_DONT_SUBTRACT_BACKGROUND; }
    else if(*tok == "nobgr") { flags &= ~FLAG_DONT_SUBTRACT_BACKGROUND; }
    else if(*tok == "pred") { flags &= ~FLAG_NO_THEORYPREDICTION_ERRORS; }
    else if(*tok == "nopred") { flags |= FLAG_NO_THEORYPREDICTION_ERRORS; }
    else if(*tok == "mcstat") { flags |= FLAG_IGNORE_MC_STATERR; }
    else { LOG(logERROR) << "Unrecognized flag \"" << *tok << "\"."; }
  }
  LOG(logINFO) << "Flags: " << massextractor::listFlags(flags);

  Double_t unfolding_mass = nominalmass;

  try {
    if(type == "yield") extract_yield(inputpath,outputpath,channels,closure,closure_sample,flags,syst,systlist,fulltake);
    else if(type == "yieldstats") extract_yield_stats(inputpath,outputpath,channels,closure,closure_sample,flags,syst,systlist, fulltake);
    else if(type == "diffxsstats") extract_diffxsec_stats(inputpath,outputpath,channels,unfolding_mass,flags,syst,systlist,fulltake);
    else extract_diffxsec(inputpath,outputpath,channels,unfolding_mass,flags,syst,systlist,fulltake);
  }
  catch(...) {
    LOG(logCRITICAL) << "Mass extraction failed.";
    return -1;
  }

  return 0;
}
