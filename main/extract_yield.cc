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

  syst_on_nominal.push_back("MATCH");
  syst_on_nominal.push_back("SCALE");
  syst_on_nominal.push_back("HAD");
  syst_on_nominal.push_back("CR");
  syst_on_nominal.push_back("UE");

  systematics.push_back("JES");

  /*systematics.push_back("JES_MPF");
  systematics.push_back("JES_INTERCAL");
  systematics.push_back("JES_UNCORR");
  systematics.push_back("JES_BJES");

  systematics.push_back("JES_FLAVOR_GLUON");
  systematics.push_back("JES_FLAVOR_QUARK");
  systematics.push_back("JES_FLAVOR_CHARM");
  systematics.push_back("JES_FLAVOR_BOTTOM");*/

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
    Double_t var_stat_pos2, var_stat_neg2;
    if(syst) {
      // We are processing systematic variations. Maybe user requested histograms for all of them?
      if(fulltake) { flags |= FLAG_STORE_HISTOGRAMS; }

      // Systematic Variations with own samples:
      for(std::vector<TString>::iterator syst = syst_on_nominal.begin(); syst != syst_on_nominal.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << " variation...";

	extractorYieldOtherSamples * variation_up = new extractorYieldOtherSamples(*ch,"Nominal", inputpath, outputpath, flags,(*syst)+"_UP");
	extractorYieldOtherSamples * variation_down = new extractorYieldOtherSamples(*ch,"Nominal", inputpath, outputpath, flags,(*syst)+"_DOWN");
	if(closure) {
	  variation_up->setClosureSample(closure_sample);
	  variation_down->setClosureSample(closure_sample);
	}

	Double_t topmass_variation_up = variation_up->getTopMass();
	Double_t topmass_variation_down = variation_down->getTopMass();

	variation_up->getStatError(var_stat_pos,var_stat_neg);
	variation_down->getStatError(var_stat_pos2,var_stat_neg2);

	LOG(logINFO) << *syst << "_UP - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_up << " +" << var_stat_pos << " -" << var_stat_neg;
	LOG(logINFO) << *syst << "_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_down << " +" << var_stat_pos2 << " -" << var_stat_neg2;

	Double_t delta_up = (Double_t)topmass_variation_up-topmass;
	Double_t delta_down = (Double_t)topmass_variation_down-topmass;

	LOG(logRESULT) << *syst << "_UP: delta = " << delta_up << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
	LOG(logRESULT) << *syst << "_DOWN: delta = " << delta_down << " GeV +" << systStatErr(total_stat_pos,var_stat_pos2) << "-" << systStatErr(total_stat_neg,var_stat_neg2);

	if(syst->Contains("MATCH") || syst->Contains("SCALE")) {
	  getSystematicUpDownError(delta_up,delta_down,total_theo_pos,total_theo_neg);
	  LOG(logDEBUG) << "Added to theory uncertainties.";
	}
	else {
	  getSystematicUpDownError(delta_up,delta_down,total_syst_pos,total_syst_neg);
	  LOG(logDEBUG) << "Added to experimental uncertainties.";
	}

	//if(syst->Contains("UP")) { SystOutputFile << systab->writeSystematicsTableUp(*syst, delta, var_stat_pos, var_stat_neg);	}
	//else { SystOutputFile << systab->writeSystematicsTableDown(delta, var_stat_pos, var_stat_neg); }
	delete variation_up;
	delete variation_down;
      }

      // Getting DrellYan/Background variations:
      LOG(logDEBUG) << "Getting BG variation...";
      extractorYieldBackground * bg_up = new extractorYieldBackground(*ch,"Nominal",inputpath, outputpath, flags,"BG_UP");
      extractorYieldBackground * bg_down = new extractorYieldBackground(*ch,"Nominal",inputpath, outputpath, flags,"BG_DOWN");
      if(closure) {
	bg_up->setClosureSample(closure_sample);
	bg_down->setClosureSample(closure_sample);
      }

      Double_t topmass_variation_bg_up = bg_up->getTopMass();
      Double_t topmass_variation_bg_down = bg_down->getTopMass();

      bg_up->getStatError(var_stat_pos,var_stat_neg);
      bg_down->getStatError(var_stat_pos2,var_stat_neg2);

      LOG(logINFO) << "BG_UP - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_bg_up << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "BG_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_bg_down << " +" << var_stat_pos2 << " -" << var_stat_neg2;

      Double_t delta_bg_up = (Double_t)topmass_variation_bg_up-topmass;
      Double_t delta_bg_down = (Double_t)topmass_variation_bg_down-topmass;

      LOG(logRESULT) << "BG_UP: delta = " << delta_bg_up << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
      LOG(logRESULT) << "BG_DOWN: delta = " << delta_bg_down << " GeV +" << systStatErr(total_stat_pos,var_stat_pos2) << "-" << systStatErr(total_stat_neg,var_stat_neg2);

      getSystematicUpDownError(delta_bg_up,delta_bg_down,total_syst_pos,total_syst_neg);
      
      //if(syst->Contains("UP")) SystOutputFile << systab->writeSystematicsTableUp(*syst, delta, var_stat_pos, var_stat_neg);
      //else SystOutputFile << systab->writeSystematicsTableDown(delta, var_stat_pos, var_stat_neg);

      // Getting PDF variations:
      LOG(logDEBUG) << "Getting PDF variation...";
      extractorYieldScaled * pdf_up   = new extractorYieldScaled(*ch,"Nominal",inputpath, outputpath, flags,"PDF_UP");
      extractorYieldScaled * pdf_down = new extractorYieldScaled(*ch,"Nominal",inputpath, outputpath, flags,"PDF_DOWN");
      if(closure) {
	bg_up->setClosureSample(closure_sample);
	bg_down->setClosureSample(closure_sample);
      }

      Double_t topmass_variation_pdf_up = pdf_up->getTopMass();
      Double_t topmass_variation_pdf_down = pdf_down->getTopMass();

      pdf_up->getStatError(var_stat_pos,var_stat_neg);
      pdf_down->getStatError(var_stat_pos2,var_stat_neg2);

      LOG(logINFO) << "PDF_UP - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_pdf_up << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "PDF_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_pdf_down << " +" << var_stat_pos2 << " -" << var_stat_neg2;

      Double_t delta_pdf_up = (Double_t)topmass_variation_pdf_up-topmass;
      Double_t delta_pdf_down = (Double_t)topmass_variation_pdf_down-topmass;

      LOG(logRESULT) << "PDF_UP: delta = " << delta_pdf_up << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
      LOG(logRESULT) << "PDF_DOWN: delta = " << delta_pdf_down << " GeV +" << systStatErr(total_stat_pos,var_stat_pos2) << "-" << systStatErr(total_stat_neg,var_stat_neg2);

      getSystematicUpDownError(delta_pdf_up,delta_pdf_down,total_syst_pos,total_syst_neg);


      // Systematic Variations produced by varying nominal samples:
      Double_t btag_syst_pos = 0, btag_syst_neg = 0;
      Double_t jes_syst_pos = 0, jes_syst_neg = 0;

      Double_t btag_stat_pos = 0, btag_stat_neg = 0;
      Double_t stat_up_pos = 0,   stat_up_neg = 0;
      Double_t stat_down_pos = 0, stat_down_neg = 0;

      for(std::vector<TString>::iterator syst = systematics.begin(); syst != systematics.end(); ++syst) {
	LOG(logDEBUG) << "Getting " << (*syst) << "_UP/DOWN variation...";

	// Prepare the extraction modules for up and down variation:
	extractorYield * variation_up = new extractorYield(*ch,(*syst)+"_UP",inputpath, outputpath, flags);
	extractorYield * variation_down = new extractorYield(*ch,(*syst)+"_DOWN",inputpath, outputpath, flags);
	if(closure) {
	  variation_up->setClosureSample(closure_sample);
	  variation_down->setClosureSample(closure_sample);
	}

	// Extract the top quark mass from data:
	Double_t topmass_variation_up = variation_up->getTopMass();
	Double_t topmass_variation_down = variation_down->getTopMass();
	variation_up->getStatError(stat_up_pos,stat_up_neg);
	variation_down->getStatError(stat_down_pos,stat_down_neg);

	LOG(logINFO) << *syst << "_UP   - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_up << " +" << stat_up_pos << " -" << stat_up_neg;
	LOG(logINFO) << *syst << "_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_down << " +" << stat_down_pos << " -" << stat_down_neg;

	// Calculate difference to nominal extraction:
	Double_t delta_up = (Double_t)topmass_variation_up-topmass;
	Double_t delta_down = (Double_t)topmass_variation_down-topmass;

	LOG(logRESULT) << *syst << "_UP:   delta = " << delta_up << " GeV +" << systStatErr(total_stat_pos,stat_up_pos) << "-" << systStatErr(total_stat_neg,stat_up_neg);
	LOG(logRESULT) << *syst << "_DOWN: delta = " << delta_down << " GeV +" << systStatErr(total_stat_pos,stat_down_pos) << "-" << systStatErr(total_stat_neg,stat_down_neg);

	bool systUpDown;
	if(syst->Contains("BTAG")) { getSystematicUpDownError(delta_up,delta_down,btag_syst_pos,btag_syst_neg); }
	else if(syst->Contains("JES_FLAVOR")) { getSystematicUpDownError(delta_up,delta_down,jes_syst_pos,jes_syst_neg); }
	else { systUpDown = getSystematicUpDownError(delta_up,delta_down,total_syst_pos,total_syst_neg); }

	if(!syst->Contains("BTAG") && !syst->Contains("JES_FLAVOR")) {
	  if(systUpDown) { 
	    if(abs(delta_up) > abs(delta_down)) { delta_down = 0; }
	    else { delta_up = delta_down; delta_down = 0; }
	  }

	  SystOutputFile << systab->writeSystematicsTableUp(*syst, delta_up, 
							    systStatErr(total_stat_pos,var_stat_pos),
							    systStatErr(total_stat_neg,var_stat_neg));
	  SystOutputFile << systab->writeSystematicsTableDown(delta_down,
							      systStatErr(total_stat_pos,var_stat_pos),
							      systStatErr(total_stat_neg,var_stat_neg));
	}
      }
      SystOutputFile << systab->writeSystematicsTableUp("JES_FLAVOR", TMath::Sqrt(jes_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
		     << systab->writeSystematicsTableDown(-1*TMath::Sqrt(jes_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
      SystOutputFile << systab->writeSystematicsTableUp("BTAG", TMath::Sqrt(btag_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
		     << systab->writeSystematicsTableDown(-1*TMath::Sqrt(btag_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
      total_syst_pos += jes_syst_pos + btag_syst_pos;
      total_syst_neg += jes_syst_neg + btag_syst_neg;
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
