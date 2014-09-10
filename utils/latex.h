#ifndef LATEX_H
#define LATEX_H

#include "TROOT.h"
namespace latex {

  class table {
  public:
    table() {};
    ~table() {};

    std::string writeSystematicsTableUp(TString syst, Double_t delta, Double_t stat_up = 0, Double_t stat_down = 0);
    std::string writeSystematicsTableDown(Double_t delta, Double_t stat_up, Double_t stat_down);
    std::string writeSystematicsTableUpDown(TString syst, Double_t delta, Double_t stat_pos = 0, Double_t stat_neg = 0);
    std::string writeSystematicsTableSummary(TString channel, Double_t topmass, Double_t total_stat_pos, Double_t total_stat_neg, Double_t total_syst_pos, Double_t total_syst_neg);

  private:
    Double_t temp_up;
    Double_t temp_down;
  };
}

#endif /*EXTRACT_H*/
