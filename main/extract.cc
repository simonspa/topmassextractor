#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdint.h>

#include "extractor.h"
#include "log.h"
#include "helpers.h"

#include <TMath.h>
#include <TSystem.h>

using namespace std;
using namespace massextractor;
using namespace unilog;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

void extract_yield(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags) {

  // #######################################
  // ###              YIELD              ###
  // #######################################

  // Do the same (or similar things) for the systematic uncertainties:
  std::vector<TString> syst_on_nominal;
  syst_on_nominal.push_back("MATCH_UP");
  syst_on_nominal.push_back("MATCH_DOWN");
  syst_on_nominal.push_back("SCALE_UP");
  syst_on_nominal.push_back("SCALE_DOWN");
  syst_on_nominal.push_back("HAD_UP");
  syst_on_nominal.push_back("HAD_DOWN");
  syst_on_nominal.push_back("CR_UP");
  syst_on_nominal.push_back("CR_DOWN");
  syst_on_nominal.push_back("UE_UP");
  syst_on_nominal.push_back("UE_DOWN");
  
  std::vector<TString> syst_bg;
  syst_bg.push_back("BG_UP");
  syst_bg.push_back("BG_DOWN");

  std::vector<TString> systematics;
  systematics.push_back("JES_UP"); systematics.push_back("JES_DOWN");
  systematics.push_back("JER_UP"); systematics.push_back("JER_DOWN");
  systematics.push_back("PU_UP"); systematics.push_back("PU_DOWN");
  systematics.push_back("TRIG_UP"); systematics.push_back("TRIG_DOWN");
  systematics.push_back("KIN_UP"); systematics.push_back("KIN_DOWN");
  systematics.push_back("LEPT_UP"); systematics.push_back("LEPT_DOWN");
  systematics.push_back("BTAG_UP"); systematics.push_back("BTAG_DOWN");
  systematics.push_back("BTAG_LJET_UP"); systematics.push_back("BTAG_LJET_DOWN");
  systematics.push_back("BTAG_PT_UP"); systematics.push_back("BTAG_PT_DOWN");
  systematics.push_back("BTAG_ETA_UP"); systematics.push_back("BTAG_ETA_DOWN");
  systematics.push_back("BTAG_LJET_PT_UP"); systematics.push_back("BTAG_LJET_PT_DOWN");
  systematics.push_back("BTAG_LJET_ETA_UP"); systematics.push_back("BTAG_LJET_ETA_DOWN");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream SystOutputFile(outputpath + "/MassFitRatesSystematics_" + *ch + ".txt", std::ofstream::trunc);
    SystOutputFile << "Top Mass, Channel: " << *ch << endl;
    SystOutputFile << "Systematic & Syst. error on m_t & [GeV] \\\\" << endl;
    SystOutputFile << "\\hline" << std::endl;


    extractorYield * mass_samples = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_STORE_HISTOGRAMS);
    if(closure) mass_samples->setClosureSample(closure_sample);

    Double_t topmass = mass_samples->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_samples->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;

    // Systematic Variations with own samples:
    for(std::vector<TString>::iterator syst = syst_on_nominal.begin(); syst != syst_on_nominal.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorYieldOtherSamples * matchscale_samples = new extractorYieldOtherSamples(*ch,"Nominal", inputpath, outputpath, flags,(*syst));
      if(closure) matchscale_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = matchscale_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("UP")) SystOutputFile << matchscale_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(2) << std::fixed << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      delete matchscale_samples;
    }

    // Getting DrellYan/Background variations:
    for(std::vector<TString>::iterator syst = syst_bg.begin(); syst != syst_bg.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorYieldBackground * bg_samples = new extractorYieldBackground(*ch,"Nominal",inputpath, outputpath, flags,(*syst));
      if(closure) bg_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = bg_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;
      
      if(syst->Contains("UP")) SystOutputFile << bg_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(2) << std::fixed << (delta > 0 ? "+" : "" ) << delta << "}_{";
      else SystOutputFile << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
    }

    // Systematic Variations produced by varying nominal samples:
    Double_t btag_syst_pos = 0, btag_syst_neg = 0;
    for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorYield * variation_samples = new extractorYield(*ch,*syst,inputpath, outputpath, flags);
      if(closure) variation_samples->setClosureSample(closure_sample);

      Double_t topmass_variation = variation_samples->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("BTAG")) {
	if(delta > 0) btag_syst_pos += delta*delta;
	else btag_syst_neg += delta*delta;
      }
      else {
	if(syst->Contains("UP")) SystOutputFile << variation_samples->getSampleLabel((*syst)) << " & $^{" << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}_{";
	else SystOutputFile << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      }
    }
    SystOutputFile << mass_samples->getSampleLabel("BTAG") << " & $^{+" << setprecision(2) << std::fixed <<  TMath::Sqrt(btag_syst_pos) << "}_{"
		   << setprecision(2) << std::fixed << TMath::Sqrt(btag_syst_neg) << "}$ \\\\" << endl;

    total_syst_pos = TMath::Sqrt(total_syst_pos);
    total_syst_neg = TMath::Sqrt(total_syst_neg);
    LOG(logRESULT) << "Channel " << *ch << ": m_t = " << setprecision(2) << std::fixed << topmass << setprecision(2) << std::fixed << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV";

    SystOutputFile << "Channel " << *ch << ": m_t = " << setprecision(2) << std::fixed << topmass << setprecision(2) << std::fixed << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV" << endl;
    SystOutputFile.close();
    delete mass_samples;
  }
  return;
}

void extract_diffxsec(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags) {

  std::vector<TString> systematics;
  systematics.push_back("MATCH_UP"); systematics.push_back("MATCH_DOWN");
  systematics.push_back("SCALE_UP"); systematics.push_back("SCALE_DOWN");
  // If we only take one mass for unfolding, we need to add this as systematic error:
  if((flags & FLAG_UNFOLD_ALLMASSES) == 0) {
    systematics.push_back("MASS_UP"); systematics.push_back("MASS_DOWN");
  }
  systematics.push_back("BG_UP"); systematics.push_back("BG_DOWN");
  systematics.push_back("JES_UP"); systematics.push_back("JES_DOWN");
  systematics.push_back("JER_UP"); systematics.push_back("JER_DOWN");
  systematics.push_back("PU_UP"); systematics.push_back("PU_DOWN");
  systematics.push_back("TRIG_UP"); systematics.push_back("TRIG_DOWN");
  systematics.push_back("KIN_UP"); systematics.push_back("KIN_DOWN");
  systematics.push_back("LEPT_UP"); systematics.push_back("LEPT_DOWN");
  systematics.push_back("BTAG_UP"); systematics.push_back("BTAG_DOWN");
  systematics.push_back("BTAG_LJET_UP"); systematics.push_back("BTAG_LJET_DOWN");
  systematics.push_back("BTAG_PT_UP"); systematics.push_back("BTAG_PT_DOWN");
  systematics.push_back("BTAG_ETA_UP"); systematics.push_back("BTAG_ETA_DOWN");
  systematics.push_back("BTAG_LJET_PT_UP"); systematics.push_back("BTAG_LJET_PT_DOWN");
  systematics.push_back("BTAG_LJET_ETA_UP"); systematics.push_back("BTAG_LJET_ETA_DOWN");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream DiffSystOutputFile(outputpath + "/MassFitDiffXSecSystematics_" + *ch + ".txt", std::ofstream::trunc);
    DiffSystOutputFile << "Top Mass, Channel: " << *ch << endl;
    DiffSystOutputFile << "Systematic & Syst. error on m_t & [GeV] \\\\" << endl;
    DiffSystOutputFile << "\\hline" << std::endl;

    extractorDiffXSec * mass_diffxs = new extractorDiffXSec(*ch,"Nominal", inputpath, outputpath, flags | FLAG_STORE_HISTOGRAMS);
    mass_diffxs->setUnfoldingMass(unfoldingMass);
    Double_t topmass = mass_diffxs->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_diffxs->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;
    
    // Systematic Variations with own samples:
    extractorDiffXSec * var, * var2;
    Double_t diff = 0;

    LOG(logDEBUG) << "Getting HAD variation...";
    var = new extractorDiffXSec(*ch,"POWHEG", inputpath, outputpath, flags);
    var2 = new extractorDiffXSec(*ch,"MCATNLO", inputpath, outputpath, flags);
    diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("HAD_UP") << " & $\\pm " << setprecision(2) << std::fixed << diff << "$ \\\\" << endl;
    LOG(logINFO) << "HAD - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2;

    LOG(logDEBUG) << "Getting CR variation...";
    var = new extractorDiffXSec(*ch,"PERUGIA11NoCR", inputpath, outputpath, flags);
    var2 = new extractorDiffXSec(*ch,"PERUGIA11", inputpath, outputpath, flags);
    diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("CR_UP") << " & $\\pm " << setprecision(2) << std::fixed << diff << "$ \\\\" << endl;
    LOG(logINFO) << "CR - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2;

    LOG(logDEBUG) << "Getting UE variation...";
    extractorDiffXSec * var3;
    var = new extractorDiffXSec(*ch,"PERUGIA11mpiHi", inputpath, outputpath, flags);
    var2 = new extractorDiffXSec(*ch,"PERUGIA11TeV", inputpath, outputpath, flags);
    var3 = new extractorDiffXSec(*ch,"PERUGIA11", inputpath, outputpath, flags);
    diff = (TMath::Abs(var->getTopMass() - var3->getTopMass()) + TMath::Abs(var2->getTopMass() - var3->getTopMass()))/2;
    DiffSystOutputFile << mass_diffxs->getSampleLabel("UE_UP") << " & $\\pm " << setprecision(2) << std::fixed << diff << "$ \\\\" << endl;
    LOG(logINFO) << "UE - " << *ch << ": delta = " << diff;
    total_syst_pos += diff*diff;
    total_syst_neg += diff*diff;
    delete var; delete var2; delete var3;

    LOG(logDEBUG) << "Getting PDF variation...";
    extractorDiffXSecScaled * pdf;
    pdf = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags, "PDF_UP");
    diff = (Double_t)topmass - pdf->getTopMass();
    LOG(logINFO) << "PDF_UP - " << *ch << ": delta = " << diff;
    if(diff > 0) total_syst_pos += diff*diff;
    else total_syst_neg += diff*diff;
    DiffSystOutputFile << pdf->getSampleLabel("PDF_UP") << " & $^{" << setprecision(2) << std::fixed << (diff > 0 ? "+" : "" ) << diff << "}_{";

    pdf = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags, "PDF_DOWN");
    diff = (Double_t)topmass - pdf->getTopMass();
    LOG(logINFO) << "PDF_DOWN - " << *ch << ": delta = " << diff;
    if(diff > 0) total_syst_pos += diff*diff;
    else total_syst_neg += diff*diff;
    DiffSystOutputFile << setprecision(2) << std::fixed <<  (diff > 0 ? "+" : "" ) << diff << "}$ \\\\" << endl;

    // Systematic Variations produced by varying nominal samples:
    Double_t btag_syst_pos = 0, btag_syst_neg = 0;
    for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
      LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
      extractorDiffXSec * variation_diffxs = new extractorDiffXSec(*ch,*syst,inputpath, outputpath, flags);
      variation_diffxs->setUnfoldingMass(unfoldingMass);

      Double_t topmass_variation = variation_diffxs->getTopMass();
      LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
      Double_t delta = (Double_t)topmass-topmass_variation;
      LOG(logINFO) << *syst << ": delta = " << delta;
      LOG(logINFO) << *syst << ": stat. unc. = " << variation_diffxs->getStatError();
      if(delta > 0) total_syst_pos += delta*delta;
      else total_syst_neg += delta*delta;

      if(syst->Contains("BTAG")) {
	if(delta > 0) btag_syst_pos += delta*delta;
	else btag_syst_neg += delta*delta;
      }
      else {
	if(syst->Contains("UP")) DiffSystOutputFile << variation_diffxs->getSampleLabel((*syst)) << " & $^{" << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}_{";
	else DiffSystOutputFile << setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta << "}$ \\\\" << endl;
      }
    }
    DiffSystOutputFile << mass_diffxs->getSampleLabel("BTAG") << " & $^{+" << setprecision(2) << std::fixed <<  TMath::Sqrt(btag_syst_pos) << "}_{"
		       << setprecision(2) << std::fixed << TMath::Sqrt(btag_syst_neg) << "}$ \\\\" << endl;

    total_syst_pos = TMath::Sqrt(total_syst_pos);
    total_syst_neg = TMath::Sqrt(total_syst_neg);
    LOG(logRESULT) << "Channel " << *ch << ": m_t = " << setprecision(2) << std::fixed << topmass << setprecision(2) << std::fixed << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV";

    DiffSystOutputFile << "Channel " << *ch << ": m_t = " << setprecision(2) << std::fixed << topmass << setprecision(2) << std::fixed << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV" << endl;
    DiffSystOutputFile.close();
    delete mass_diffxs;
  }
  return;
}

int main(int argc, char* argv[]) {

  Log::ReportingLevel() = Log::FromString("INFO");
  TString outputpath = "ExtractionResults";
  TString inputpath = "";
  bool closure = true;
  std::string type, channel;

  std::vector<std::string> flagtokens;
  uint32_t flags = 0;

  for (int i = 1; i < argc; i++) {
    // select to either extract from yield of differential cross section:
    if (!strcmp(argv[i],"-t")) {
      type = string(argv[++i]);
      if(type !=  "diffxs" && type != "yield") {
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
    // Select the channel to run on:
    else if(!strcmp(argv[i],"-c")) { channel = string(argv[++i]); }
    // Read and tokeinze the flags::
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
  flags = FLAG_CHISQUARE_FROM_FITS | FLAG_DONT_USE_COVARIANCE | FLAG_NORMALIZE_YIELD;
  // Check and assign the flags:
  for(std::vector<std::string>::iterator tok = flagtokens.begin(); tok != flagtokens.end(); ++tok) {
    if(*tok == "fit") { flags |= FLAG_CHISQUARE_FROM_FITS; }
    else if(*tok == "pdf") { flags |= FLAG_STORE_PDFS; }
    else if(*tok == "root") { flags |= FLAG_STORE_HISTOGRAMS; }
    else if(*tok == "cov") { flags &= ~FLAG_DONT_USE_COVARIANCE; }
    else if(*tok == "nocov") { flags |= FLAG_DONT_USE_COVARIANCE; }
    else if(*tok == "norm") { flags |= FLAG_NORMALIZE_YIELD; }
    else if(*tok == "nonorm") { flags &= ~FLAG_NORMALIZE_YIELD; }
    else if(*tok == "lastbin") { flags |= FLAG_LASTBIN_EXTRACTION; }
    else if(*tok == "bgr") { flags |= FLAG_DONT_SUBTRACT_BACKGROUND; }
    else { LOG(logERROR) << "Unrecognized flag \"" << *tok << "\"."; }
  }
  LOG(logINFO) << "Flags: " << massextractor::listFlags(flags);

  TString closure_sample = "Nominal";
  Double_t unfolding_mass = nominalmass;

  try {
    if(type == "yield") extract_yield(inputpath,outputpath,channels,closure,closure_sample,flags);
    else extract_diffxsec(inputpath,outputpath,channels,unfolding_mass,flags);
  }
  catch(...) {
    LOG(logCRITICAL) << "Mass extraction failed.";
    return -1;
  }

  return 0;
}
