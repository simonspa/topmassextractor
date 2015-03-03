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


void extract_yield(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake) {

  // #######################################
  // ###              YIELD              ###
  // #######################################

  // Do the same (or similar things) for the systematic uncertainties:
  std::vector<TString> syst_on_nominal;
  std::vector<TString> syst_bg;
  std::vector<TString> systematics;

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

  syst_bg.push_back("BG_UP");
  syst_bg.push_back("BG_DOWN");

  systematics.push_back("JES_UP"); systematics.push_back("JES_DOWN");

  /*systematics.push_back("JES_MPF_UP"); systematics.push_back("JES_MPF_DOWN");
  systematics.push_back("JES_INTERCAL_UP"); systematics.push_back("JES_INTERCAL_DOWN");
  systematics.push_back("JES_UNCORR_UP"); systematics.push_back("JES_UNCORR_DOWN");
  systematics.push_back("JES_BJES_UP"); systematics.push_back("JES_BJES_DOWN");

  systematics.push_back("JES_FLAVOR_GLUON_UP"); systematics.push_back("JES_FLAVOR_GLUON_DOWN");
  systematics.push_back("JES_FLAVOR_QUARK_UP"); systematics.push_back("JES_FLAVOR_QUARK_DOWN");
  systematics.push_back("JES_FLAVOR_CHARM_UP"); systematics.push_back("JES_FLAVOR_CHARM_DOWN");
  systematics.push_back("JES_FLAVOR_BOTTOM_UP"); systematics.push_back("JES_FLAVOR_BOTTOM_DOWN");*/

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
    latex::table * systab = new latex::table();
    SystOutputFile << "Top Mass, Channel: " << *ch << endl;
    SystOutputFile << "Systematic & Syst. error on m_t [GeV] & Stat. error on Variation \\\\" << endl;
    SystOutputFile << "\\hline" << std::endl;

    extractorYield * mass_samples = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_STORE_HISTOGRAMS);
    if(closure) mass_samples->setClosureSample(closure_sample);

    Double_t topmass = mass_samples->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_samples->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;
    LOG(logRESULT) << "Nominal: mass = " << topmass << " GeV +" << total_stat_pos << "-" << total_stat_neg;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;
    Double_t total_theo_pos = 0;
    Double_t total_theo_neg = 0;
    Double_t var_stat_pos, var_stat_neg;
    if(syst) {
      // We are processing systematic variations. Maybe user requested histograms for all of them?
      if(fulltake) { flags |= FLAG_STORE_HISTOGRAMS; }

      // Systematic Variations with own samples:
      for(std::vector<TString>::iterator syst = syst_on_nominal.begin(); syst != syst_on_nominal.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";

	extractorYieldOtherSamples * matchscale_samples = new extractorYieldOtherSamples(*ch,"Nominal", inputpath, outputpath, flags,(*syst));
	if(closure) matchscale_samples->setClosureSample(closure_sample);

	Double_t topmass_variation = matchscale_samples->getTopMass();
	LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
	Double_t delta = (Double_t)topmass-topmass_variation;
	matchscale_samples->getStatError(var_stat_pos,var_stat_neg);
	LOG(logRESULT) << *syst << ": delta = " << delta << " GeV +" << var_stat_pos << "-" << var_stat_neg;
	
	if(syst->Contains("MATCH") || syst->Contains("SCALE")) {
	  if(delta > 0) total_theo_pos += delta*delta;
	  else total_theo_neg += delta*delta;
	}
	else {
	  if(delta > 0) total_syst_pos += delta*delta;
	  else total_syst_neg += delta*delta;
	}

	if(syst->Contains("UP")) { SystOutputFile << systab->writeSystematicsTableUp(*syst, delta, var_stat_pos, var_stat_neg);	}
	else { SystOutputFile << systab->writeSystematicsTableDown(delta, var_stat_pos, var_stat_neg); }
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
	bg_samples->getStatError(var_stat_pos,var_stat_neg);
	LOG(logRESULT) << *syst << ": delta = " << delta << " GeV +" << var_stat_pos << "-" << var_stat_neg;
	if(delta > 0) total_syst_pos += delta*delta;
	else total_syst_neg += delta*delta;
      
	if(syst->Contains("UP")) SystOutputFile << systab->writeSystematicsTableUp(*syst, delta, var_stat_pos, var_stat_neg);
	else SystOutputFile << systab->writeSystematicsTableDown(delta, var_stat_pos, var_stat_neg);
      }

      // Systematic Variations produced by varying nominal samples:
      Double_t btag_syst_pos = 0, btag_syst_neg = 0;
      Double_t btag_stat_pos = 0, btag_stat_neg = 0;
      Double_t jes_flavor_pos = 0, jes_flavor_neg = 0;
      Double_t jes_flavor_stat_pos = 0, jes_flavor_stat_neg = 0;

      Double_t total_jes_pos = 0, total_jes_neg = 0;
      for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";
	extractorYield * variation_samples = new extractorYield(*ch,*syst,inputpath, outputpath, flags);
	if(closure) variation_samples->setClosureSample(closure_sample);

	Double_t topmass_variation = variation_samples->getTopMass();
	LOG(logINFO) << *syst << " - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation;
	Double_t delta = (Double_t)topmass-topmass_variation;
	variation_samples->getStatError(var_stat_pos,var_stat_neg);
	LOG(logRESULT) << *syst << ": delta = " << delta << " GeV +" << var_stat_pos << "-" << var_stat_neg;

	if(syst->Contains("FLAVOR")) {
	  // Adding up linearly, preserving the sign:
	  if(delta > 0) jes_flavor_pos += delta;
	  else jes_flavor_neg += delta;
	  jes_flavor_stat_pos += var_stat_pos;
	  jes_flavor_stat_neg += var_stat_neg;
	}
	else {
	  // All others are summed quadratically directly to the total uncertainty:
	  if(delta > 0) total_syst_pos += delta*delta;
	  else total_syst_neg += delta*delta;
	}

	if(syst->Contains("BTAG")) {
	  if(delta > 0) btag_syst_pos += delta*delta;
	  else btag_syst_neg += delta*delta;
	  btag_stat_pos += var_stat_pos;
	  btag_stat_neg += var_stat_neg;
	}
	else if(syst->Contains("FLAVOR")) {}
	else {
	  if(syst->Contains("JES") && *syst != "JES_UP" && *syst != "JES_DOWN") {
	    LOG(logCRITICAL) << "Adding " << *syst;
	    if(delta > 0) total_jes_pos += delta*delta;
	    else total_jes_neg += delta*delta;
	  }
	  if(syst->Contains("UP")) SystOutputFile << systab->writeSystematicsTableUp(*syst, delta, var_stat_pos, var_stat_neg);
	  else SystOutputFile << systab->writeSystematicsTableDown(delta, var_stat_pos, var_stat_neg);
	}
      }
      SystOutputFile << systab->writeSystematicsTableUp("JES_FLAVOR", TMath::Sqrt(jes_flavor_pos), jes_flavor_stat_pos/4, jes_flavor_stat_neg/4)
		     << systab->writeSystematicsTableDown(TMath::Sqrt(jes_flavor_neg), jes_flavor_stat_pos/4, jes_flavor_stat_neg/4);

      // Add JES_FLAVOR to the total uncertainty:
      total_syst_pos += jes_flavor_pos*jes_flavor_pos;
      total_syst_neg += jes_flavor_neg*jes_flavor_neg;
      total_jes_pos += jes_flavor_pos*jes_flavor_pos;
      total_jes_neg += jes_flavor_neg*jes_flavor_neg;

      // Also writing the summed JES uncertainty for comparison, not added again (is already...)
      SystOutputFile << systab->writeSystematicsTableUp("JES", TMath::Sqrt(total_jes_pos), 0, 0)
		     << systab->writeSystematicsTableDown(TMath::Sqrt(total_jes_neg), 0, 0);

      SystOutputFile << systab->writeSystematicsTableUp("BTAG", TMath::Sqrt(btag_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
		     << systab->writeSystematicsTableDown(TMath::Sqrt(btag_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
    }

    total_syst_pos = TMath::Sqrt(total_syst_pos);
    total_syst_neg = TMath::Sqrt(total_syst_neg);
    total_theo_pos = TMath::Sqrt(total_theo_pos);
    total_theo_neg = TMath::Sqrt(total_theo_neg);

    LOG(logRESULT) << systab->writeSystematicsTableSummary(*ch, topmass, total_stat_pos, total_stat_neg, total_syst_pos, total_syst_neg, total_theo_pos, total_theo_neg);
    SystOutputFile << systab->writeSystematicsTableSummary(*ch, topmass, total_stat_pos, total_stat_neg, total_syst_pos, total_syst_neg, total_theo_pos, total_theo_neg);

    SystOutputFile.close();
    delete mass_samples;
  }
  return;
}

void extract_yield_stats(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, bool fulltake) {

  std::vector<TString> syst_bg;
  //syst_bg.push_back("BG_UP"); syst_bg.push_back("BG_DOWN");

  std::vector<TString> systematic;
  systematic.push_back("MATCH_UP"); systematic.push_back("MATCH_DOWN");
  systematic.push_back("SCALE_UP"); systematic.push_back("SCALE_DOWN");
  /*systematic.push_back("HAD_UP"); systematic.push_back("HAD_DOWN");*/
  systematic.push_back("JES_UP"); systematic.push_back("JES_DOWN");
  systematic.push_back("JER_UP"); systematic.push_back("JER_DOWN");
  systematic.push_back("PU_UP"); systematic.push_back("PU_DOWN");
  systematic.push_back("TRIG_UP"); systematic.push_back("TRIG_DOWN");
  systematic.push_back("KIN_UP"); systematic.push_back("KIN_DOWN");
  systematic.push_back("LEPT_UP"); systematic.push_back("LEPT_DOWN");
  systematic.push_back("BTAG_UP"); systematic.push_back("BTAG_DOWN");
  systematic.push_back("BTAG_LJET_UP"); systematic.push_back("BTAG_LJET_DOWN");
  systematic.push_back("BTAG_PT_UP"); systematic.push_back("BTAG_PT_DOWN");
  systematic.push_back("BTAG_ETA_UP"); systematic.push_back("BTAG_ETA_DOWN");
  systematic.push_back("BTAG_LJET_PT_UP"); systematic.push_back("BTAG_LJET_PT_DOWN");
  systematic.push_back("BTAG_LJET_ETA_UP"); systematic.push_back("BTAG_LJET_ETA_DOWN");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream SystOutputFile(outputpath + "/MassFitRatesSystematicStats_" + *ch + ".txt", std::ofstream::trunc);
    latex::table * systab = new latex::table();
    SystOutputFile << "Top Mass, Channel: " << *ch << endl;
    SystOutputFile << "Systematic & & Stat. error on Systematic \\\\" << endl;
    SystOutputFile << "\\hline" << std::endl;

    Double_t normal_stat_pos, normal_stat_neg;
    Double_t inf_stat_pos, inf_stat_neg;
    Double_t btag_stat_pos = 0, btag_stat_neg = 0;

    if(syst) {
      // Color Reconnection:
      extractorYield * normal, * infstat;
      normal = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags);
      if(closure) normal->setClosureSample("PERUGIA11");
      normal->getTopMass(); normal->getStatError(normal_stat_pos,normal_stat_neg);
      infstat = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS);
      if(closure) infstat->setClosureSample("PERUGIA11");
      infstat->getTopMass(); infstat->getStatError(inf_stat_pos,inf_stat_neg);
      LOG(logRESULT) << "CR_UP - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
      SystOutputFile << systab->writeSystematicsTableUp("CR_UP", 0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg));
      delete normal; delete infstat;
      normal = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags);
      if(closure) normal->setClosureSample("PERUGIA11NoCR");
      normal->getTopMass(); normal->getStatError(normal_stat_pos,normal_stat_neg);
      infstat = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS);
      if(closure) infstat->setClosureSample("PERUGIA11NoCR");
      infstat->getTopMass(); infstat->getStatError(inf_stat_pos,inf_stat_neg);
      LOG(logRESULT) << "CR_DOWN - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
      SystOutputFile << systab->writeSystematicsTableDown(0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg));
      delete normal; delete infstat;

      // Underlying Event
      normal = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags);
      if(closure) normal->setClosureSample("PERUGIA11");
      normal->getTopMass(); normal->getStatError(normal_stat_pos,normal_stat_neg);
      infstat = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS);
      if(closure) infstat->setClosureSample("PERUGIA11");
      infstat->getTopMass(); infstat->getStatError(inf_stat_pos,inf_stat_neg);
      LOG(logRESULT) << "UE_UP - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
      SystOutputFile << systab->writeSystematicsTableUp("UE_UP", 0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg));
      delete normal; delete infstat;
      normal = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags);
      if(closure) normal->setClosureSample("PERUGIA11TeV");
      normal->getTopMass(); normal->getStatError(normal_stat_pos,normal_stat_neg);
      infstat = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS);
      if(closure) infstat->setClosureSample("PERUGIA11TeV");
      infstat->getTopMass(); infstat->getStatError(inf_stat_pos,inf_stat_neg);
      LOG(logRESULT) << "UE_DOWN - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
      SystOutputFile << systab->writeSystematicsTableDown(0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg));
      delete normal; delete infstat;

      // Backgrounds first:
      for(std::vector<TString>::iterator syst = syst_bg.begin(); syst != syst_bg.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";

	extractorYieldBackground * bg_normal = new extractorYieldBackground(*ch,"Nominal",inputpath, outputpath, flags,(*syst));
	if(closure) bg_normal->setClosureSample("Nominal");
	bg_normal->getTopMass();
	bg_normal->getStatError(normal_stat_pos,normal_stat_neg);

	extractorYieldBackground * bg_inf = new extractorYieldBackground(*ch,"Nominal",inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS,(*syst));
	if(closure) bg_inf->setClosureSample("Nominal");
	bg_inf->getTopMass();
	bg_inf->getStatError(inf_stat_pos,inf_stat_neg);

	LOG(logRESULT) << *syst << " - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
	if(syst->Contains("UP")) { SystOutputFile << systab->writeSystematicsTableUp(*syst, 0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg)); }
	else { SystOutputFile << systab->writeSystematicsTableDown(0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg)); }

	delete bg_normal;
	delete bg_inf;
      }

      // Other Systematic Variations
      for(std::vector<TString>::iterator syst = systematic.begin(); syst != systematic.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";

	extractorYield * systematic_normal = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags);
	if(closure) systematic_normal->setClosureSample((*syst));
	systematic_normal->getTopMass();
	systematic_normal->getStatError(normal_stat_pos,normal_stat_neg);

	extractorYield * systematic_infstat = new extractorYield(*ch,"Nominal", inputpath, outputpath, flags | FLAG_INFINITE_DATA_STATISTICS);
	if(closure) systematic_infstat->setClosureSample((*syst));
	systematic_infstat->getTopMass();
	systematic_infstat->getStatError(inf_stat_pos,inf_stat_neg);

	if(syst->Contains("BTAG")) {
	  LOG(logINFO) << *syst << " - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
	  btag_stat_pos += systStatErr(normal_stat_pos,inf_stat_pos);
	  btag_stat_neg += systStatErr(normal_stat_neg,inf_stat_neg);
	}
	else {
	  LOG(logRESULT) << *syst << " - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
	  if(syst->Contains("UP")) { SystOutputFile << systab->writeSystematicsTableUp(*syst, 0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg)); }
	  else { SystOutputFile << systab->writeSystematicsTableDown(0, systStatErr(normal_stat_pos,inf_stat_pos), systStatErr(normal_stat_neg,inf_stat_neg)); }
	}
	delete systematic_normal;
	delete systematic_infstat;
      }

      LOG(logRESULT) << "BTAG - " << *ch << ": stat. unc: +" << btag_stat_pos/12 << "-" << btag_stat_neg/12;

      SystOutputFile << systab->writeSystematicsTableUp("BTAG", 0, btag_stat_pos/12, btag_stat_neg/12)
		     << systab->writeSystematicsTableDown(0, btag_stat_pos/12, btag_stat_neg/12);

    }

    SystOutputFile.close();
  }
  return;
}

void extract_diffxsec(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake) {

  std::vector<TString> predictions;
  predictions.push_back("MATCH_UP"); predictions.push_back("MATCH_DOWN");
  predictions.push_back("SCALE_UP"); predictions.push_back("SCALE_DOWN");

  std::vector<TString> systematics;
  systematics.push_back("MATCH");
  systematics.push_back("SCALE");
  systematics.push_back("BG");
  systematics.push_back("JES");
  systematics.push_back("JER");
  systematics.push_back("PU");
  systematics.push_back("TRIG");
  systematics.push_back("KIN");
  systematics.push_back("LEPT");
  systematics.push_back("BTAG");
  systematics.push_back("BTAG_LJET");
  systematics.push_back("BTAG_PT");
  systematics.push_back("BTAG_ETA");
  systematics.push_back("BTAG_LJET_PT");
  systematics.push_back("BTAG_LJET_ETA");

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    std::ofstream DiffSystOutputFile(outputpath + "/MassFitDiffXSecSystematics_" + *ch + ".txt", std::ofstream::trunc);
    latex::table * systab = new latex::table();
    DiffSystOutputFile << "Top Mass, Channel: " << *ch << endl;
    DiffSystOutputFile << "Systematic & Syst. error on m_t [GeV] & Stat. error on Variation \\\\" << endl;
    DiffSystOutputFile << "\\hline" << std::endl;

    extractorDiffXSec * mass_diffxs = new extractorDiffXSec(*ch,"Nominal", inputpath, outputpath, flags | FLAG_STORE_HISTOGRAMS);
    mass_diffxs->setUnfoldingMass(unfoldingMass);
    Double_t topmass = mass_diffxs->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_diffxs->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;
    LOG(logRESULT) << "Nominal: mass = " << topmass << " GeV +" << total_stat_pos << "-" << total_stat_neg;

    Double_t total_syst_pos = 0;
    Double_t total_syst_neg = 0;
    Double_t total_theo_pos = 0;
    Double_t total_theo_neg = 0;
    Double_t var_stat_pos, var_stat_neg;
    Double_t var_stat_pos2, var_stat_neg2;
    if(syst) {
      // We are processing systematic variations. Maybe user requested histograms for all of them?
      if(fulltake) { flags |= FLAG_STORE_HISTOGRAMS; }
 
      // Theory prediction errors - this uses Nominal as data and shifts MC according to the given Syst. Sample:
      for(std::vector<TString>::iterator pred = predictions.begin(); pred != predictions.end(); ++pred) {
	LOG(logDEBUG) << "Getting Theory Prediction error for " << (*pred) << "...";
	extractorDiffXSecPrediction * variation_diffxs = new extractorDiffXSecPrediction(*ch,"Nominal",inputpath, outputpath, flags,*pred);
	variation_diffxs->setUnfoldingMass(unfoldingMass);

	Double_t topmass_variation = variation_diffxs->getTopMass();
	variation_diffxs->getStatError(var_stat_pos,var_stat_neg);
	LOG(logINFO) << *pred << "_PRED - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation << " +" << var_stat_pos << " -" << var_stat_neg;
	Double_t delta = (Double_t)topmass-topmass_variation;
	LOG(logRESULT) << *pred << "_PRED: delta = " << delta << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
	if(delta > 0) total_theo_pos += delta*delta;
	else total_theo_neg += delta*delta;

	if(pred->Contains("UP")) DiffSystOutputFile << systab->writeSystematicsTableUp((*pred)+"_PRED", delta, 
										       systStatErr(total_stat_pos,var_stat_pos),
										       systStatErr(total_stat_neg,var_stat_neg));
	else DiffSystOutputFile << systab->writeSystematicsTableDown(delta,
								     systStatErr(total_stat_pos,var_stat_pos),
								     systStatErr(total_stat_neg,var_stat_neg));
      }

      // Systematic Variations with own samples:
      extractorDiffXSec * var, * var2;
      Double_t diff = 0;

      LOG(logDEBUG) << "Getting HAD variation...";
      var = new extractorDiffXSec(*ch,"POWHEG", inputpath, outputpath, flags);
      var2 = new extractorDiffXSec(*ch,"MCATNLO", inputpath, outputpath, flags);
      diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
      var->getStatError(var_stat_pos,var_stat_neg);
      var2->getStatError(var_stat_pos2,var_stat_neg2);
      LOG(logINFO) << "POWHEG - " << *ch << ": minimum Chi2 @ m_t=" << var->getTopMass() << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "MCATNLO - " << *ch << ": minimum Chi2 @ m_t=" << var2->getTopMass() << " +" << var_stat_pos2 << " -" << var_stat_neg2;
      DiffSystOutputFile << systab->writeSystematicsTableUpDown("HAD_UP", diff, 
								systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2),
								systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2));
      LOG(logRESULT) << "HAD- " << *ch << ": delta = " << diff << " GeV "
		     << "+" << systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2)
		     << "-" << systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2);
      total_syst_pos += diff*diff;
      total_syst_neg += diff*diff;
      delete var; delete var2;

      LOG(logDEBUG) << "Getting CR variation...";
      var = new extractorDiffXSec(*ch,"PERUGIA11NoCR", inputpath, outputpath, flags);
      var2 = new extractorDiffXSec(*ch,"PERUGIA11", inputpath, outputpath, flags);
      diff = TMath::Abs(var->getTopMass()-var2->getTopMass())/2;
      var->getStatError(var_stat_pos,var_stat_neg);
      var2->getStatError(var_stat_pos2,var_stat_neg2);
      LOG(logINFO) << "PERUGIA11NoCR - " << *ch << ": minimum Chi2 @ m_t=" << var->getTopMass() << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "PERUGIA11 - " << *ch << ": minimum Chi2 @ m_t=" << var2->getTopMass() << " +" << var_stat_pos2 << " -" << var_stat_neg2;
      DiffSystOutputFile << systab->writeSystematicsTableUpDown("CR_UP", diff,
								systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2),
								systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2));
      LOG(logRESULT) << "CR - " << *ch << ": delta = " << diff << " GeV "
		     << "+" << systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2)
		     << "-" << systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2);
      total_syst_pos += diff*diff;
      total_syst_neg += diff*diff;
      delete var; delete var2;

      LOG(logDEBUG) << "Getting UE variation...";
      extractorDiffXSec * var3;
      var = new extractorDiffXSec(*ch,"PERUGIA11mpiHi", inputpath, outputpath, flags);
      var2 = new extractorDiffXSec(*ch,"PERUGIA11TeV", inputpath, outputpath, flags);
      var3 = new extractorDiffXSec(*ch,"PERUGIA11", inputpath, outputpath, flags);
      diff = (TMath::Abs(var->getTopMass() - var3->getTopMass()) + TMath::Abs(var2->getTopMass() - var3->getTopMass()))/2;
      Double_t var_stat_pos3, var_stat_neg3;
      var->getStatError(var_stat_pos,var_stat_neg);
      var2->getStatError(var_stat_pos2,var_stat_neg2);
      var3->getStatError(var_stat_pos3,var_stat_neg3);
      LOG(logINFO) << "PERUGIA11mpiHi - " << *ch << ": minimum Chi2 @ m_t=" << var->getTopMass() << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "PERUGIA11TeV - " << *ch << ": minimum Chi2 @ m_t=" << var2->getTopMass() << " +" << var_stat_pos2 << " -" << var_stat_neg2;
      LOG(logINFO) << "PERUGIA11 - " << *ch << ": minimum Chi2 @ m_t=" << var3->getTopMass() << " +" << var_stat_pos3 << " -" << var_stat_neg3;
      DiffSystOutputFile << systab->writeSystematicsTableUpDown("UE_UP", diff,
								systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2)+systStatErr(total_stat_pos,var_stat_pos3),
								systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2)+systStatErr(total_stat_neg,var_stat_neg3));
      LOG(logRESULT) << "UE - " << *ch << ": delta = " << diff << " GeV "
		     << "+" << systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos2)+systStatErr(total_stat_pos,var_stat_pos3)
		     << "-" << systStatErr(total_stat_neg,var_stat_neg)+systStatErr(total_stat_neg,var_stat_neg2)+systStatErr(total_stat_neg,var_stat_neg3);
      total_syst_pos += diff*diff;
      total_syst_neg += diff*diff;
      delete var; delete var2; delete var3;

      LOG(logDEBUG) << "Getting PDF variation...";
      extractorDiffXSecScaled * pdf;
      pdf = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags, "PDF_UP");
      diff = (Double_t)topmass - pdf->getTopMass();
      pdf->getStatError(var_stat_pos,var_stat_neg);
      LOG(logINFO) << "PDF_UP - " << *ch << ": minimum Chi2 @ m_t=" << pdf->getTopMass() << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logRESULT) << "PDF_UP - " << *ch << ": delta = " << diff << " GeV "
		     << "+" << systStatErr(total_stat_pos,var_stat_pos)+systStatErr(total_stat_pos,var_stat_pos)
		     << "-" << systStatErr(total_stat_neg,var_stat_neg);
      if(diff > 0) total_syst_pos += diff*diff;
      else total_syst_neg += diff*diff;
      DiffSystOutputFile << systab->writeSystematicsTableUp("PDF_UP", diff,
							    systStatErr(total_stat_pos,var_stat_pos),
							    systStatErr(total_stat_neg,var_stat_neg));

      pdf = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags, "PDF_DOWN");
      diff = (Double_t)topmass - pdf->getTopMass();
      pdf->getStatError(var_stat_pos,var_stat_neg);
      LOG(logINFO) << "PDF_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << pdf->getTopMass() << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logRESULT) << "PDF_DOWN - " << *ch << ": delta = " << diff << " GeV "
		     << "+" << systStatErr(total_stat_pos,var_stat_pos)
		     << "-" << systStatErr(total_stat_neg,var_stat_neg);
      if(diff > 0) total_syst_pos += diff*diff;
      else total_syst_neg += diff*diff;
      DiffSystOutputFile << systab->writeSystematicsTableDown(diff,
							      systStatErr(total_stat_pos,var_stat_pos),
							      systStatErr(total_stat_neg,var_stat_neg));

      // Systematic Variations produced by varying nominal samples:
      Double_t btag_syst_pos = 0, btag_syst_neg = 0;
      Double_t btag_stat_pos = 0, btag_stat_neg = 0;
      Double_t stat_up_pos = 0, stat_up_neg = 0;
      Double_t stat_down_pos = 0, stat_down_neg = 0;
      for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << "_UP/DOWN variation...";
	
	// Prepare the extraction modules for up and down variation:
	extractorDiffXSec * variation_up = new extractorDiffXSec(*ch,*syst + "_UP",inputpath, outputpath, flags);
	extractorDiffXSec * variation_down = new extractorDiffXSec(*ch,*syst + "_DOWN",inputpath, outputpath, flags);
	variation_up->setUnfoldingMass(unfoldingMass);
	variation_down->setUnfoldingMass(unfoldingMass);

	// Extract the top quark mass from data:
	Double_t topmass_variation_up = variation_up->getTopMass();
	Double_t topmass_variation_down = variation_down->getTopMass();
	variation_up->getStatError(stat_up_pos,stat_up_neg);
	variation_down->getStatError(stat_down_pos,stat_down_neg);

	LOG(logINFO) << *syst << "_UP   - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_up << " +" << stat_up_pos << " -" << stat_up_neg;
	LOG(logINFO) << *syst << "_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_down << " +" << stat_down_pos << " -" << stat_down_neg;

	// Calculate difference to nominal extraction:
	Double_t delta_up = (Double_t)topmass-topmass_variation_up;
	Double_t delta_down = (Double_t)topmass-topmass_variation_down;

	LOG(logRESULT) << *syst << "_UP:   delta = " << delta_up << " GeV +" << systStatErr(total_stat_pos,stat_up_pos) << "-" << systStatErr(total_stat_neg,stat_up_neg);
	LOG(logRESULT) << *syst << "_DOWN: delta = " << delta_down << " GeV +" << systStatErr(total_stat_pos,stat_down_pos) << "-" << systStatErr(total_stat_neg,stat_down_neg);

	// Check if variation go in the same direction:
	if((delta_up < 0) == (delta_down < 0)) {
	  // Variations have the same sign, so just take the maximum:
	  Double_t delta = 0;
	  if(abs(delta_up) > abs(delta_down)) { delta = delta_up; }
	  else { delta = delta_down; }
	  
	  if(delta > 0) { total_syst_pos += delta*delta; }
	  else { total_syst_neg += delta*delta; }
	  LOG(logINFO) << "Counted maximum deviation only: delta = " << delta << " GeV";
	}
	else {
	  // Variations have different sign, add them both:
	  if(delta_up > 0) { total_syst_pos += delta_up*delta_up; total_syst_neg += delta_down*delta_down; }
	  else { total_syst_pos += delta_down*delta_down; total_syst_neg += delta_up*delta_up; }
	  LOG(logINFO) << "Counted both variations.";
	}

	if(syst->Contains("BTAG")) {
	  if(delta > 0) btag_syst_pos += delta*delta;
	  else btag_syst_neg += delta*delta;
	  btag_stat_pos += systStatErr(total_stat_pos,var_stat_pos);
	  btag_stat_neg += systStatErr(total_stat_neg,var_stat_neg);
	}
	else {
	  if(syst->Contains("UP")) DiffSystOutputFile << systab->writeSystematicsTableUp(*syst, delta, 
											 systStatErr(total_stat_pos,var_stat_pos),
											 systStatErr(total_stat_neg,var_stat_neg));
	  else DiffSystOutputFile << systab->writeSystematicsTableDown(delta,
								       systStatErr(total_stat_pos,var_stat_pos),
								       systStatErr(total_stat_neg,var_stat_neg));
	}
      }
      DiffSystOutputFile << systab->writeSystematicsTableUp("BTAG", TMath::Sqrt(btag_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
			 << systab->writeSystematicsTableDown(TMath::Sqrt(btag_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
    }

    total_syst_pos = TMath::Sqrt(total_syst_pos);
    total_syst_neg = TMath::Sqrt(total_syst_neg);
    total_theo_pos = TMath::Sqrt(total_theo_pos);
    total_theo_neg = TMath::Sqrt(total_theo_neg);

    LOG(logRESULT) << systab->writeSystematicsTableSummary(*ch, topmass, total_stat_pos, total_stat_neg, total_syst_pos, total_syst_neg,total_theo_pos,total_theo_neg);
    DiffSystOutputFile << systab->writeSystematicsTableSummary(*ch, topmass, total_stat_pos, total_stat_neg, total_syst_pos, total_syst_neg,total_theo_pos,total_theo_neg);

    DiffSystOutputFile.close();
    delete mass_diffxs;
  }
  return;
}

void extract_diffxsec_stats(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, bool fulltake) {

  // here we need to use two different input paths since the statistical errors have already been changed (for the inf. stat. data set)
  // before the unfolding (appending "_infstat" to the input path)

  std::vector<TString> systematics;
  systematics.push_back("MATCH_UP"); systematics.push_back("MATCH_DOWN");
  systematics.push_back("SCALE_UP"); systematics.push_back("SCALE_DOWN");
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

    std::ofstream DiffSystOutputFile(outputpath + "/MassFitDiffXSecSystematicsStats_" + *ch + ".txt", std::ofstream::trunc);
    latex::table * systab = new latex::table();
    DiffSystOutputFile << "Top Mass, Channel: " << *ch << endl;
    DiffSystOutputFile << "Systematic &  & Stat. error on Variation \\\\" << endl;
    DiffSystOutputFile << "\\hline" << std::endl;

    Double_t normal_stat_pos, normal_stat_neg;
    Double_t inf_stat_pos, inf_stat_neg;

    if(syst) {
      // Systematic Variations produced by varying nominal samples:
      Double_t btag_stat_pos = 0, btag_stat_neg = 0;
      
      for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";

	extractorDiffXSec * diffxs_normal = new extractorDiffXSec(*ch,*syst,inputpath, outputpath, flags);
	diffxs_normal->setUnfoldingMass(unfoldingMass);
	diffxs_normal->getTopMass();
	diffxs_normal->getStatError(normal_stat_pos,normal_stat_neg);

	extractorDiffXSec * diffxs_inf = new extractorDiffXSec(*ch,*syst,inputpath + "_infstat", outputpath, flags);
	diffxs_inf->setUnfoldingMass(unfoldingMass);
	diffxs_inf->getTopMass();
	diffxs_inf->getStatError(inf_stat_pos,inf_stat_neg);

	delete diffxs_normal;
	delete diffxs_inf;

	if(syst->Contains("BTAG")) {
	  LOG(logINFO) << *syst << " - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
	  btag_stat_pos += systStatErr(normal_stat_pos,inf_stat_pos);
	  btag_stat_neg += systStatErr(normal_stat_neg,inf_stat_neg);
	}
	else {
	  LOG(logRESULT) << *syst << " - " << *ch << ": stat. unc: +" << systStatErr(normal_stat_pos,inf_stat_pos) << "-" << systStatErr(normal_stat_neg,inf_stat_neg);
	  if(syst->Contains("UP")) DiffSystOutputFile << systab->writeSystematicsTableUp(*syst, 0, 
											 systStatErr(normal_stat_pos,inf_stat_pos),
											 systStatErr(normal_stat_neg,inf_stat_neg));
	  else DiffSystOutputFile << systab->writeSystematicsTableDown(0,
								       systStatErr(normal_stat_pos,inf_stat_pos),
								       systStatErr(normal_stat_neg,inf_stat_neg));
	}
      }

      LOG(logRESULT) << "BTAG - " << *ch << ": stat. unc: +" << btag_stat_pos/12 << "-" << btag_stat_neg/12;
      DiffSystOutputFile << systab->writeSystematicsTableUp("BTAG", 0, btag_stat_pos/12, btag_stat_neg/12)
			 << systab->writeSystematicsTableDown(0, btag_stat_pos/12, btag_stat_neg/12);
    }

    DiffSystOutputFile.close();
  }
  return;
}

int main(int argc, char* argv[]) {

  Log::ReportingLevel() = Log::FromString("INFO");
  FILE* logfile;

  TString outputpath = "ExtractionResults";
  TString inputpath = "";
  std::string type, channel;

  std::vector<std::string> flagtokens;
  uint32_t flags = 0;
  bool syst = false;

  // If this bool is set, PDFs for *all* systematics are produced (many!)
  bool fulltake = false;

  bool closure = true;
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
    // Set systematics flag to run on systematic variation sample
    else if(!strcmp(argv[i],"-s")) { syst = true; }
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
  flags = FLAG_CHISQUARE_FROM_FITS | FLAG_DONT_USE_COVARIANCE | FLAG_NORMALIZE_DISTRIBUTIONS | FLAG_DONT_SUBTRACT_BACKGROUND | FLAG_NO_THEORYPREDICTION_ERRORS;
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
    if(type == "yield") extract_yield(inputpath,outputpath,channels,closure,closure_sample,flags,syst,fulltake);
    else if(type == "yieldstats") extract_yield_stats(inputpath,outputpath,channels,closure,closure_sample,flags,syst,fulltake);
    else if(type == "diffxsstats") extract_diffxsec_stats(inputpath,outputpath,channels,unfolding_mass,flags,syst,fulltake);
    else extract_diffxsec(inputpath,outputpath,channels,unfolding_mass,flags,syst,fulltake);
  }
  catch(...) {
    LOG(logCRITICAL) << "Mass extraction failed.";
    return -1;
  }

  return 0;
}
