#include "TFile.h"
#include <iostream>
#include <vector>
#include "TH1D.h"
#include <exception>
#include <fstream>

using namespace std;

void calcYieldPDFUncertainty() {

  TString path = "preunfolded/PDF_";
  TString path2 = "/Nominal/";
  TString file = "/HypTTBar1stJetMass_UnfoldingHistos.root";

  TString histoname = "aRecHist";

  // All channels:
  std::vector<TString> channels;
  channels.push_back("ee");
  channels.push_back("emu");
  channels.push_back("mumu");

  // For combined we need to sum all channels:
  std::vector<Double_t> combined_central;
  std::vector< vector< Double_t > > combined_variations_up;
  std::vector< vector< Double_t > > combined_variations_down;


  for(std::vector<TString>::iterator ch = channels.begin(); ch != channels.end(); ch++) {

    std::cout << "Processing channel " << *ch << std::endl;

    // File for the central value:
    TFile * pdf_0_central = TFile::Open(path+"0_CENTRAL"+path2+*ch+file,"read");
    if(!pdf_0_central->IsOpen()) {
      std::cout << "Problem opening file \"" << path << "0_CENTRAL" << path2 << *ch << file << "!" << endl;
      throw 1;    
    }

    TH1D *central_pdf = (TH1D*)pdf_0_central->Get(histoname);
    //central_pdf->Scale(1./central_pdf->Integral("width"));

    // Vector for central values:
    std::vector<Double_t> central;

    // Vector for all bins containing vectors for all variations:
    std::vector< vector< Double_t > > variations_up;
    std::vector< vector< Double_t > > variations_down;
    
    // Get binning:
    std::vector<Double_t> binning;

    Int_t nbins = central_pdf->GetNbinsX();

    for (Int_t bin = 1; bin <= nbins; bin++) { 
      binning.push_back(central_pdf->GetBinLowEdge(bin));
      central.push_back(central_pdf->GetBinContent(bin));

      // combined:
      if(bin-1 < combined_central.size()) { combined_central.at(bin-1) += central_pdf->GetBinContent(bin); }
      else { combined_central.push_back(central_pdf->GetBinContent(bin)); }

      variations_up.push_back(std::vector<Double_t>());
      variations_down.push_back(std::vector<Double_t>());

      // combined:
      if(bin-1 >= combined_variations_up.size()) { combined_variations_up.push_back(std::vector<Double_t>()); }
      if(bin-1 >= combined_variations_down.size()) { combined_variations_down.push_back(std::vector<Double_t>()); }

    }
    binning.push_back(central_pdf->GetBinLowEdge(nbins) + central_pdf->GetBinWidth(nbins));
    
    std::cout << "Central PDF has " << central.size() << " bins." << endl;
    std::cout << "UP Variations prepared to have " << variations_up.size() << " bins." << endl;
    std::cout << "DOWN Variations prepared to have " << variations_down.size() << " bins." << endl;

    // Now iterate over the other 26 UP and DOWN uncertainties:
    for(Int_t i = 1; i < 27; i++) {
      TString name;
      name.Form("%s%d",path.Data(),i);
      TFile * pdf_i_up = TFile::Open(name+"_UP"+path2+*ch+file,"read");
      TFile * pdf_i_down = TFile::Open(name+"_DOWN"+path2+*ch+file,"read");
      if(!pdf_i_up->IsOpen() || !pdf_i_down->IsOpen()) {
	std::cout << "either problem opening file \"" << name << "_UP" << path2 << *ch << file << "!" << endl;
	std::cout << "or problem opening file \"" << path << "_DOWN" << path2 << *ch << file << "!" << endl;
	throw 1;
      }

      TH1D *pdf_up = (TH1D*)pdf_i_up->Get(histoname);
      TH1D *pdf_down = (TH1D*)pdf_i_down->Get(histoname);
      //pdf_up->Scale(1./pdf_up->Integral("width"));
      //pdf_down->Scale(1./pdf_down->Integral("width"));

      for (Int_t bin = 0; bin < nbins; bin++) { 
	//std::cout << "Bin " << bin << " got " << pdf_up->GetBinContent(bin) << endl;
	variations_up.at(bin).push_back(pdf_up->GetBinContent(bin+1));
	variations_down.at(bin).push_back(pdf_down->GetBinContent(bin+1));

	// combined:
	if(i <= combined_variations_up.at(bin).size()) { combined_variations_up.at(bin).at(i-1) += pdf_up->GetBinContent(bin+1); }
	else { combined_variations_up.at(bin).push_back(pdf_up->GetBinContent(bin+1)); }
	if(i <= combined_variations_down.at(bin).size()) { combined_variations_down.at(bin).at(i-1) += pdf_down->GetBinContent(bin+1); }
	else { combined_variations_down.at(bin).push_back(pdf_down->GetBinContent(bin+1)); }
      }
    }

    // Analyse data, sum:
    std::cout << " Have " << variations_up.at(0).size() << " variations (UP) in " << variations_up.size() << " bins." << endl;
    std::cout << " Have " << variations_down.at(0).size() << " variations (DOWN) in " << variations_down.size() << " bins." << endl;

    ofstream outfile ("pdf_uncert_"+*ch+".txt");

    // For all bins:
    for(size_t bin = 0; bin < central.size(); bin++) {
      // sum the difference between central and all variations:
      Double_t sum_up = 0;
      for(std::vector<Double_t>::iterator var = variations_up.at(bin).begin(); var != variations_up.at(bin).end(); var++) {
	//std::cout << "Bin " << bin << " up add " << central.at(bin) << "-" << *var << "=" << (central.at(bin)-*var) << endl;
	sum_up += fabs(central.at(bin)-*var);
      }
      // sum the difference between central and all variations:
      Double_t sum_down = 0;
      for(std::vector<Double_t>::iterator var = variations_down.at(bin).begin(); var != variations_down.at(bin).end(); var++) {
	//std::cout << "Bin " << bin << " down add " << central.at(bin) << "-" << *var << "=" << (central.at(bin)-*var) << endl;
	sum_down += fabs(central.at(bin)-*var);
      }

      std::cout << "Bin " << bin << " up   central " << central.at(bin) << " sum " << sum_up;
      sum_up /= variations_up.at(bin).size();
      std::cout << " mean " << sum_up << " rel " << (100*sum_up/central.at(bin)) << "%" << endl;
      
      std::cout << "Bin " << bin << " down   central " << central.at(bin) << " sum " << sum_down;
      sum_down /= variations_down.at(bin).size();
      std::cout << " mean " << sum_down << " rel " << (100*sum_down/central.at(bin)) << "%" << endl;

      // Symmetrize:
      Double_t scalefactor = (sum_up + sum_down)/(2*central.at(bin));
      outfile << scalefactor << endl;
    }

    outfile.close();
  }

  // And finally: combined:
  // Analyse data, sum:
  std::cout << " Have " << combined_variations_up.at(0).size() << " variations (UP) in " << combined_variations_up.size() << " bins." << endl;
  std::cout << " Have " << combined_variations_down.at(0).size() << " variations (DOWN) in " << combined_variations_down.size() << " bins." << endl;
  
  ofstream coutfile ("pdf_uncert_combined.txt");

  // For all bins:
  for(size_t bin = 0; bin < combined_central.size(); bin++) {
    // sum the difference between central and all variations:
    Double_t sum_up = 0;
    for(std::vector<Double_t>::iterator var = combined_variations_up.at(bin).begin(); var != combined_variations_up.at(bin).end(); var++) {
      //std::cout << "Bin " << bin << " up add " << central.at(bin) << "-" << *var << "=" << (central.at(bin)-*var) << endl;
      sum_up += fabs(combined_central.at(bin)-*var);
    }
    // sum the difference between central and all variations:
    Double_t sum_down = 0;
    for(std::vector<Double_t>::iterator var = combined_variations_down.at(bin).begin(); var != combined_variations_down.at(bin).end(); var++) {
      //std::cout << "Bin " << bin << " down add " << central.at(bin) << "-" << *var << "=" << (central.at(bin)-*var) << endl;
      sum_down += fabs(combined_central.at(bin)-*var);
    }
    
    std::cout << "Bin " << bin << " up   central " << combined_central.at(bin) << " sum " << sum_up;
    sum_up /= combined_variations_up.at(bin).size();
    std::cout << " mean " << sum_up << " rel " << (100*sum_up/combined_central.at(bin)) << "%" << endl;
    
    std::cout << "Bin " << bin << " down   central " << combined_central.at(bin) << " sum " << sum_down;
    sum_down /= combined_variations_down.at(bin).size();
    std::cout << " mean " << sum_down << " rel " << (100*sum_down/combined_central.at(bin)) << "%" << endl;

    // Symmetrize:
    Double_t scalefactor = (sum_up + sum_down)/(2*combined_central.at(bin));
    coutfile << scalefactor << endl;
  }

  coutfile.close();
}
