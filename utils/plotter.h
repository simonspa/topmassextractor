#include <TStyle.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLegend.h>

namespace massextractor {

  void setHHStyle(TStyle& HHStyle);

  // Draw label for Decay Channel in upper left corner of plot
  void DrawDecayChLabel(TString decaychannel="", Int_t bin=0, std::vector<Double_t> boundaries=std::vector<Double_t>(), int cmsprelim=1, double textSize=0.04);

  // Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
  void DrawCMSLabels(double lumi=19712, double energy=8, double textSize=0.045);

  void setStyle(TGraphErrors *hist, TString name="");
  void setStyleAndFillLegend(TGraphErrors* hist, TString name, TLegend *leg, bool closure=false);
  void setLegendStyle(TLegend *leg);

  void rescaleGraph(TGraphErrors * g, Double_t up=1.2, Double_t down=0.95);
  TString getChannelLabel(TString channel);
}
