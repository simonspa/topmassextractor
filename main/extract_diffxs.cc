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

void massextractor::extract_diffxsec(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, std::vector<std::string> systlist, bool fulltake) {

  std::vector<TString> predictions;
  if((flags & FLAG_USE_NLO) == 0) {
    predictions.push_back("MATCH");
  }
  predictions.push_back("SCALE");

  std::vector<TString> systematics;
  if((flags & FLAG_USE_NLO) != 0) {
    systematics.push_back("MATCH");
    systematics.push_back("SCALE");
  }
  systematics.push_back("BG");
  systematics.push_back("JES");
  /*systematics.push_back("JES_MPF");
  systematics.push_back("JES_INTERCAL");
  systematics.push_back("JES_UNCORR");
  systematics.push_back("JES_BJES");

  systematics.push_back("JES_FLAVOR_GLUON");
  systematics.push_back("JES_FLAVOR_QUARK");
  systematics.push_back("JES_FLAVOR_CHARM");
  systematics.push_back("JES_FLAVOR_BOTTOM");
  */
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
	extractorDiffXSec * variation_diffxs_up, * variation_diffxs_down;

	// Coherently vary MadGraph Prediction and DiffXS up and down for LO:
	if((flags & FLAG_USE_NLO) == 0) {
	  variation_diffxs_up = new extractorDiffXSecPrediction(*ch,(*pred)+"_UP",inputpath, outputpath, flags,(*pred)+"_UP");
	  variation_diffxs_down = new extractorDiffXSecPrediction(*ch,(*pred)+"_DOWN",inputpath, outputpath, flags,(*pred)+"_DOWN");
	}
	// Compare Powheg NLO prediction only agains nominal DiffXSec:
	else {
	  variation_diffxs_up = new extractorDiffXSecGenLevelPrediction(*ch,"Nominal",inputpath, outputpath, flags,(*pred)+"_UP");
	  variation_diffxs_down = new extractorDiffXSecGenLevelPrediction(*ch,(*pred)+"Nominal",inputpath, outputpath, flags,(*pred)+"_DOWN");
	}

	variation_diffxs_up->setUnfoldingMass(unfoldingMass);
	variation_diffxs_down->setUnfoldingMass(unfoldingMass);

	Double_t topmass_variation_up = variation_diffxs_up->getTopMass();
	Double_t topmass_variation_down = variation_diffxs_down->getTopMass();

	variation_diffxs_up->getStatError(var_stat_pos,var_stat_neg);
	variation_diffxs_down->getStatError(var_stat_pos2,var_stat_neg2);

	LOG(logINFO) << *pred << "_UP_PRED - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_up << " +" << var_stat_pos << " -" << var_stat_neg;
	LOG(logINFO) << *pred << "_DOWN_PRED - " << *ch << ": minimum Chi2 @ m_t=" << topmass_variation_down << " +" << var_stat_pos2 << " -" << var_stat_neg2;

	Double_t delta_up = (Double_t)topmass_variation_up-topmass;
	Double_t delta_down = (Double_t)topmass_variation_down-topmass;

	LOG(logRESULT) << *pred << "_UP_PRED: delta = " << delta_up << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
	LOG(logRESULT) << *pred << "_DOWN_PRED: delta = " << delta_down << " GeV +" << systStatErr(total_stat_pos,var_stat_pos2) << "-" << systStatErr(total_stat_neg,var_stat_neg2);

	getSystematicUpDownError(delta_up,delta_down,total_theo_pos,total_theo_neg);

	/*
	  if(pred->Contains("UP")) DiffSystOutputFile << systab->writeSystematicsTableUp((*pred)+"_PRED", delta, 
										       systStatErr(total_stat_pos,var_stat_pos),
										       systStatErr(total_stat_neg,var_stat_neg));
	else DiffSystOutputFile << systab->writeSystematicsTableDown(delta,
								     systStatErr(total_stat_pos,var_stat_pos),
								     systStatErr(total_stat_neg,var_stat_neg));
	*/
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
      getSystematicUpDownError(diff,-1*diff,total_syst_pos,total_syst_neg);
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
      getSystematicUpDownError(diff,-1*diff,total_syst_pos,total_syst_neg);
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
      getSystematicUpDownError(diff,-1*diff,total_syst_pos,total_syst_neg);
      delete var; delete var2; delete var3;

      LOG(logDEBUG) << "Getting PDF variation...";
      extractorDiffXSecScaled * pdf_up   = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags | FLAG_DONT_USE_COVARIANCE, "PDF_UP");
      extractorDiffXSecScaled * pdf_down = new extractorDiffXSecScaled(*ch,"Nominal", inputpath, outputpath, flags | FLAG_DONT_USE_COVARIANCE, "PDF_DOWN");

      Double_t topmass_pdf_up = pdf_up->getTopMass();
      Double_t topmass_pdf_down = pdf_down->getTopMass();

      pdf_up->getStatError(var_stat_pos,var_stat_neg);
      pdf_up->getStatError(var_stat_pos2,var_stat_neg2);

      LOG(logINFO) << "PDF_UP   - " << *ch << ": minimum Chi2 @ m_t=" << topmass_pdf_up << " +" << var_stat_pos << " -" << var_stat_neg;
      LOG(logINFO) << "PDF_DOWN - " << *ch << ": minimum Chi2 @ m_t=" << topmass_pdf_down << " +" << var_stat_pos2 << " -" << var_stat_neg2;

      // Calculate difference to nominal extraction:
      Double_t delta_pdf_up = (Double_t)topmass_pdf_up-topmass;
      Double_t delta_pdf_down = (Double_t)topmass_pdf_down-topmass;

      LOG(logRESULT) << "PDF_UP:   delta = " << delta_pdf_up << " GeV +" << systStatErr(total_stat_pos,var_stat_pos) << "-" << systStatErr(total_stat_neg,var_stat_neg);
      LOG(logRESULT) << "PDF_DOWN: delta = " << delta_pdf_down << " GeV +" << systStatErr(total_stat_pos,var_stat_pos2) << "-" << systStatErr(total_stat_neg,var_stat_neg2);

      bool systPdfUpDown = getSystematicUpDownError(delta_pdf_up,delta_pdf_down,total_syst_pos,total_syst_neg);
      if(systPdfUpDown) {
	if(abs(delta_pdf_up) > abs(delta_pdf_down)) { delta_pdf_down = 0; }
	else { delta_pdf_up = delta_pdf_down; delta_pdf_down = 0; }
      }

      DiffSystOutputFile << systab->writeSystematicsTableUp("PDF_UP", delta_pdf_up, 
							    systStatErr(total_stat_pos,var_stat_pos),
							    systStatErr(total_stat_neg,var_stat_neg));
      DiffSystOutputFile << systab->writeSystematicsTableDown(delta_pdf_down,
							      systStatErr(total_stat_pos,var_stat_pos),
							      systStatErr(total_stat_neg,var_stat_neg));

      // Systematic Variations produced by varying nominal samples:
      Double_t btag_syst_pos = 0, btag_syst_neg = 0;
      Double_t jes_syst_pos = 0, jes_syst_neg = 0;

      Double_t btag_stat_pos = 0, btag_stat_neg = 0;
      Double_t stat_up_pos = 0,   stat_up_neg = 0;
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

	  DiffSystOutputFile << systab->writeSystematicsTableUp(*syst, delta_up, 
								systStatErr(total_stat_pos,var_stat_pos),
								systStatErr(total_stat_neg,var_stat_neg));
	  DiffSystOutputFile << systab->writeSystematicsTableDown(delta_down,
								  systStatErr(total_stat_pos,var_stat_pos),
								  systStatErr(total_stat_neg,var_stat_neg));
	}
      }
      DiffSystOutputFile << systab->writeSystematicsTableUp("JES_FLAVOR", TMath::Sqrt(jes_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
			 << systab->writeSystematicsTableDown(TMath::Sqrt(jes_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
      DiffSystOutputFile << systab->writeSystematicsTableUp("BTAG", TMath::Sqrt(btag_syst_pos), btag_stat_pos/12, btag_stat_neg/12)
			 << systab->writeSystematicsTableDown(-1*TMath::Sqrt(btag_syst_neg), btag_stat_pos/12, btag_stat_neg/12);
      total_syst_pos += jes_syst_pos + btag_syst_pos;
      total_syst_neg += jes_syst_neg + btag_syst_neg;
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

void massextractor::extract_diffxsec_stats(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags, bool syst, std::vector<std::string> systlist, bool fulltake) {

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
