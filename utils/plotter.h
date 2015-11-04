#include <TStyle.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLegend.h>
#include <TH2D.h>

namespace massextractor {

  void setHHStyle(TStyle& HHStyle);

  // Draw label for Decay Channel in upper left corner of plot
  void DrawDecayChLabel(TString decaychannel="", bool drawchannel=true, Int_t bin=0, std::vector<Double_t> boundaries=std::vector<Double_t>(), int cmsprelim=1, double textSize=0.04);

  // Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
  void DrawCMSLabels(double lumi=19712, double energy=8, double textSize=0.045);
  void DrawFreeCMSLabels(char* text, double textSize = 0.045);

  void setStyle(TGraphErrors *hist, TString name="");
  void setStyle(TH2D *hist, TString name="");
  void setStyle(TH1D *hist, TString name="");
  void setStyleAndFillLegend(TGraphErrors* hist, TString name, TLegend *leg, bool closure=false);
  void setTheoryStyleAndFillLegend(TH1* histo, TString theoryName, TLegend *leg);
  void setLegendStyle(TLegend *leg);

  void rescaleGraph(TGraphErrors * g, Double_t up=1.2, Double_t down=0.95);
  TString getChannelLabel(TString channel);

  void drawRatio(const TH1* histNumerator, const TH1* histDenominator1, const TH1* uncband,
		 TGraphAsymmErrors *ratio_stat = 0, TGraphAsymmErrors *ratio_total = 0, 
		 const TH1* histDenominator2 = 0, const TH1* histDenominator3 = 0, 
		 const TH1* histDenominator4 = 0, const TH1* histDenominator5 = 0, 
		 const TH1* histDenominator6 = 0, const TH1* histDenominator7 = 0, 
		 const Double_t& ratioMin = 0.5, const Double_t& ratioMax = 1.5, 
		 const bool data = true,
		 TStyle myStyle = *gStyle);

}
