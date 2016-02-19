#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdint.h>
#include <TRandom3.h>
#include <TF1.h>

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

#define MINMASS 163.5

void massextractor::extract_diffxsec_pseudo(TString inputpath, TString outputpath, std::vector<TString> channels, Double_t unfoldingMass, uint32_t flags) {

  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    // Set number of experiments:
    Int_t nexp = 10000;

    // Get random number generator:
    TRandom3 *myrnd = new TRandom3();

    TH1D * pseudoMasses = new TH1D("pseudomasses_" + *ch, "pseudomasses_" + *ch, 250,160,175);
    
    for(Int_t exp = 0; exp < nexp; exp++) {
      extractorDiffXSecPseudoExp * mass_diffxs = new extractorDiffXSecPseudoExp(*ch,"Nominal", inputpath, outputpath, flags, myrnd);
      mass_diffxs->setUnfoldingMass(unfoldingMass);
      Double_t topmass = mass_diffxs->getTopMass();
      if(topmass < 100) continue;

      Double_t total_stat_pos;
      Double_t total_stat_neg;
      mass_diffxs->getStatError(total_stat_pos,total_stat_neg);
      LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;
      LOG(logRESULT) << "(" << exp << "/" << nexp << ") Nominal: mass = " << topmass << " GeV +" << total_stat_pos << "-" << total_stat_neg;
      pseudoMasses->Fill(topmass);
      delete mass_diffxs;
    }

    TString fileName = "pseudomasses_" + *ch + ".root";
    TFile * out = TFile::Open(fileName,"RECREATE");
    gDirectory->pwd();

    pseudoMasses->Fit("gaus");
    TF1 *fit = pseudoMasses->GetFunction("gaus");
    fit->SetLineWidth(3);
    fit->SetLineColor(kRed+1);
    LOG(logRESULT) << "Pseudo Experiments: mass = " << fit->GetParameter(1) << " GeV +-" << fit->GetParameter(2);
    pseudoMasses->Write();


    extractorDiffXSec * mass_xs = new extractorDiffXSec(*ch,"Nominal", inputpath, outputpath, flags);
    mass_xs->setUnfoldingMass(unfoldingMass);
    Double_t topmass = mass_xs->getTopMass();
    Double_t total_stat_pos;
    Double_t total_stat_neg;
    mass_xs->getStatError(total_stat_pos,total_stat_neg);
    LOG(logINFO) << *ch << ": minimum Chi2 @ m_t=" << topmass << " +" << total_stat_pos << " -" << total_stat_neg;
    LOG(logRESULT) << "Reference, Nominal: mass = " << topmass << " GeV +" << total_stat_pos << "-" << total_stat_neg;
  }

  return;
}
