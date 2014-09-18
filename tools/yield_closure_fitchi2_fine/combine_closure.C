#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

Double_t nominalmass = 172.5;
Double_t lumi = 19712;

std::string getNextToLastLine(std::ifstream& in)
{
  std::string line;
  std::string oldline;
  std::string oldoldline;
    while (in >> std::ws && std::getline(in, line)) {
      // skip empty lines
      oldoldline = oldline;
      oldline = line
        ;
    }

    return oldoldline;
}

void combine_closure() {

  // Compile the plotter utils:
  gROOT->ProcessLine(".L ../../utils/plotter.cc+");

  bool drawSys = false;

  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");
  channels.push_back("combined");

  std::vector<TString> masses;
  masses.push_back("166");
  masses.push_back("169");
  masses.push_back("171");
  masses.push_back("172");
  masses.push_back("173");
  masses.push_back("175");
  masses.push_back("178");
  
  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ++ch) {

    massextractor::setHHStyle(*gStyle);

    TGraphAsymmErrors * graph = new TGraphAsymmErrors();
    TGraphAsymmErrors * graphtotal = new TGraphAsymmErrors();

    massextractor::setStyle(dynamic_cast<TGraphErrors*>(graph));
    massextractor::setStyle(dynamic_cast<TGraphErrors*>(graphtotal));
    
    for(std::vector<TString>::iterator m = masses.begin(); m != masses.end(); ++m) {
      TString filename = "closure" + *m + "/MassFitRatesSystematics_" + *ch + ".txt";
      std::ifstream infile(filename);

      if (infile) {
	std::stringstream stream;
	std::string dummy;
  	std::string line = getNextToLastLine(infile);
	std::cout << line << '\n';
	stream.clear();
	stream << line;
	for(size_t i = 0; i < 4; ++i ) { stream >> dummy; }
	double outmass, statp, statn, systp, systn;
	stream >> outmass >> statp >> statn >> dummy >> systp >> systn;
	//stat = atof(dummy.erase(0,2).c_str());
	double inmass = (atof(*m)+0.5);

	double tot_pos = sqrt(statp*statp + systp*systp);
	double tot_neg = sqrt(statn*statn + systn*systn);

	std::cout << outmass << "@" << inmass << " | " 
		  << statp << " " << statn << " | " 
		  << systp << " " << systn << " | "
		  << tot_pos << " " << tot_neg << std::endl;

	// Add point:
	graph->SetPoint(m-masses.begin(),inmass,outmass);
	graphtotal->SetPoint(m-masses.begin(),inmass,outmass);

	graphtotal->SetPointError(m-masses.begin(),0,0,tot_neg,tot_pos);
	graph->SetPointError(m-masses.begin(),0,0,TMath::Abs(statn),TMath::Abs(statp));
      }
      else
        std::cout << "Unable to open file " << filename << ".\n";

      infile.close();
    }

    graphtotal->GetXaxis()->SetTitle("input m_{t} #left[GeV#right]");
    graphtotal->GetYaxis()->SetTitle("output m_{t} #left[GeV#right]");
    graph->GetXaxis()->SetTitle("input m_{t} #left[GeV#right]");
    graph->GetYaxis()->SetTitle("output m_{t} #left[GeV#right]");

    Double_t m_min = 163, m_max = 183;
    graphtotal->GetXaxis()->SetLimits(m_min,m_max);
    graphtotal->GetYaxis()->SetRangeUser(m_min,m_max);
    graph->GetXaxis()->SetLimits(m_min,m_max);
    graph->GetYaxis()->SetRangeUser(m_min,m_max);

    graph->SetLineColor(kGray+1);
    graph->SetMarkerStyle(20);
    graphtotal->SetLineColor(kOrange-4);

    TLine *l = new TLine(m_min,m_min,m_max,m_max);
    l->SetLineColor(kGray+1);
    l->SetLineStyle(7);
    l->SetLineWidth(2);

    TCanvas * c = new TCanvas("yield_closure_" + *ch,"yield_closure_" + *ch);
    c->cd();
    if(drawSys) {
      graphtotal->Draw("a p e1");
      graph->Draw("same e1");
    }
    else {
      graph->Draw("a p e1");
    }

    l->Draw();
    massextractor::DrawCMSLabels();
    massextractor::DrawDecayChLabel(*ch);
    c->Print("yield_closure_" + *ch + ".pdf");
  }

  return;
}
