#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include "latex.h"
#include "helpers.h"

using namespace latex;

std::string table::writeSystematicsTableUp(TString syst, Double_t delta, Double_t stat_up, Double_t stat_down) {
  std::stringstream out;
  temp_up = stat_up; temp_down = stat_down;
  out << massextractor::getSampleLabel(syst) << " & $^{";
  if(TMath::Abs(delta) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed << (delta > 0 ? "+" : "" ) << delta; }
  out << "}_{";
  return out.str();
}

std::string table::writeSystematicsTableDown(Double_t delta, Double_t stat_up, Double_t stat_down) {
  std::stringstream out;
  Double_t stat_pos = (temp_up+stat_up)/2;
  Double_t stat_neg = (temp_down+stat_down)/2;
  if(TMath::Abs(delta) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed <<  (delta > 0 ? "+" : "" ) << delta; }
  out << "}$ & ";

  if(stat_pos > 0.015 && stat_neg > 0.015) { }
  out << "$^{+";
  if(TMath::Abs(stat_pos) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed << stat_pos; }
  out << "}_{-";
  if(TMath::Abs(stat_neg) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed << stat_neg; }
  out << "}$";
 
  out << " \\\\" << std::endl;
  temp_up = 0; temp_down = 0;
  return out.str();
}

std::string table::writeSystematicsTableUpDown(TString syst, Double_t delta, Double_t stat_pos, Double_t stat_neg) {
  std::stringstream out;
  out << massextractor::getSampleLabel(syst) << " & $";
  if(TMath::Abs(delta) < 0.01) { out << "< 0.01"; }
  else { out << "\\pm " << std::setprecision(2) << std::fixed << TMath::Abs(delta); }
  out << "$ & ";  

  if(stat_pos > 0.015 && stat_neg > 0.015) { }
  out << "$^{+";
  if(TMath::Abs(stat_pos) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed << stat_pos; }
  out << "}_{-";
  if(TMath::Abs(stat_neg) < 0.01) { out << "< 0.01"; }
  else { out << std::setprecision(2) << std::fixed << stat_neg; }
  out << "}$";
 
  out << " \\\\" << std::endl;
  return out.str();
}

std::string table::writeSystematicsTableSummary(TString channel, Double_t topmass, Double_t total_stat_pos, Double_t total_stat_neg, Double_t total_syst_pos, Double_t total_syst_neg) {
  std::stringstream out;
  out << "\\hline" << std::endl;
  out << "Stat. & $^{" << std::setprecision(2) << std::fixed << "+" << total_stat_pos << "}_{-" << total_stat_neg << "}$ & \\\\" << std::endl;
  out << "Total syst. & $^{+" << total_syst_pos << "}_{-" << total_syst_neg << "}$ & \\\\" << std::endl;
  out << "%Channel " << channel << ": m_t = " << std::setprecision(2) << std::fixed << topmass << std::setprecision(2) << std::fixed << " +" << total_stat_pos << " -" << total_stat_neg << " (stat) +" << total_syst_pos << " -" << total_syst_neg << " (syst) GeV" << std::endl;
  out << "\\hline" << std::endl;
  return out.str();
}
