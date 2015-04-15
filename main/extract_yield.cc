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

void massextractor::extract_yield(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, std::vector<std::string> systlist, bool fulltake) {

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

void massextractor::extract_yield_stats(TString inputpath, TString outputpath, std::vector<TString> channels, bool closure, TString closure_sample, uint32_t flags, bool syst, std::vector<std::string> systlist, bool fulltake) {

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
