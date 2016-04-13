
#include <iostream>
#include <stdlib.h> 

//#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TAxis.h"
#include "TMath.h"
#include "TChain.h"

#include "DataNEUT.h"

// const int nToys = 1000;
// const int nEvents = 100000;
// const int nbinsX = 8, nbinsY = 3;
// double binningX[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // Target type CH, O, Al, Cu, Fe, Zn, Pb
// double binningY[4] = {0, 200, 600, 10000}; // The binning in momentum

//void FillChi2Cov(TH2D *rates[9][4], TH2D *ratesMean[9][4], TH2D *ratesUncertainties[9][4], TMatrixD *Correlations[9][7]);

int main(int argc, char **argv){
  std::cout << "111" << std::endl;
  // TChain *input = new TChain("neut");
  // // All the files that we generated patiently with CalculateOOFVTheory.exe,
  // // note that NEUT cannot generate on H only target, this is why the output_3 doesnt exist (and we will essentially assign CH error to C and to H target,
  // // generous over estimation)
  // for(int i = 1; i < argc; i++)
  //   input->Add(Form("%s", argv[i]));
  
  int nToys = 1000;
  int nEvents = 10000;
  int nbinsX = 8, nbinsY = 3;
  double binningX[9] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5}; // Target type CH, O, Al, Cu, Fe, Zn, Pb
  double binningY[4] = {0, 200, 600, 10000}; // The binning in momentum
 
  TTree *input;
  TFile *files[7];
  files[0] = new TFile("previous/output_1.root", "READ");
  files[1] = new TFile("previous/output_2.root", "READ");
  files[2] = new TFile("previous/output_4.root", "READ");
  files[3] = new TFile("previous/output_5.root", "READ");
  files[4] = new TFile("previous/output_6.root", "READ");
  files[5] = new TFile("previous/output_7.root", "READ");
  files[6] = new TFile("previous/output_8.root", "READ");
  std::cout << "111" << std::endl;

  // What we want
  TFile *OutputFile = new TFile("output.root", "RECREATE");
  OutputFile->cd();
  OutputFile->SetOption("COLZ");
  //  TTree* input;
  TH2D *rates[9][4];
  TH2D *ratesMean[9][4];
  TH2D *ratesUncertainties[9][4];

  // The source of uncertainties are split, just as a check that everything is working like we expect
  std::string source[4];
  source[0] = "Flux";
  source[1] = "FSI";
  source[2] = "XSec";
  source[3] = "Total";
 
  // A map of all the pdg to the histgram array index
  std::map<int,int> PDGmap;
  const int nPDGs = 9;
  int PDGs[nPDGs] = {13, -13, 11, -11, 22, 211, -211, 2212, 2112};
  for(int iPDG = 0; iPDG < nPDGs; iPDG++)
    PDGmap[PDGs[iPDG]] = iPDG;
  
  
  // TMatrixD *Correlations[9][7];

  // This is part of what we want:
  // The correlations in momentum of the particles exiting nucleus
  // (this might very well turn out to be useless but why not calculating it since we are here...)
  // For now, I do only the total contribution, one could check for each individual source (especially the FSI should be quite interesting)
  // for(int iParticle; iParticle < 9; iParticle++)
  //   for(int iTarget; iTarget < 7; iTarget++)
  //     Correlations[iParticle][iTarget] = new TMatrixD(nbinsY, nbinsY);

  for(int iSource = 0; iSource < 4; iSource++){
    // This is what we want:
    // histograms with the uncertainty of every particle exiting all sort of nucleus, all this binned in momentum bins

    ratesUncertainties[0][iSource] = new TH2D(Form("%sRU_13", source[iSource].c_str()),   Form("%s Rate Uncertainty PDG = 13", source[iSource].c_str()),   
					      nbinsX, binningX, nbinsY, binningY); // muon
    ratesUncertainties[1][iSource] = new TH2D(Form("%sRU_-13", source[iSource].c_str()),  Form("%s Rate Uncertainty PDG = -13", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // antimuon
    ratesUncertainties[2][iSource] = new TH2D(Form("%sRU_11", source[iSource].c_str()),   Form("%s Rate Uncertainty PDG = 11", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // electron
    ratesUncertainties[3][iSource] = new TH2D(Form("%sRU_-11", source[iSource].c_str()),  Form("%s Rate Uncertainty PDG = -11", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // positron
    ratesUncertainties[4][iSource] = new TH2D(Form("%sRU_22", source[iSource].c_str()),   Form("%s Rate Uncertainty PDG = 22", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // gamma
    ratesUncertainties[5][iSource] = new TH2D(Form("%sRU_211", source[iSource].c_str()),  Form("%s Rate Uncertainty PDG = 211", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // pion+
    ratesUncertainties[6][iSource] = new TH2D(Form("%sRU_-211", source[iSource].c_str()), Form("%s Rate Uncertainty PDG = -211", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // pion-
    ratesUncertainties[7][iSource] = new TH2D(Form("%sRU_2212", source[iSource].c_str()), Form("%s Rate Uncertainty PDG = 2212", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // proton
    ratesUncertainties[8][iSource] = new TH2D(Form("%sRU_2112", source[iSource].c_str()), Form("%s Rate Uncertainty PDG = 2112", source[iSource].c_str()),
					      nbinsX, binningX, nbinsY, binningY); // neutron
    

    // These are the mean of the rates of all the toys
    ratesMean[0][iSource] = new TH2D(Form("%sRM_13", source[iSource].c_str()),   Form("%s Rate Mean PDG = 13", source[iSource].c_str()),   
				     nbinsX, binningX, nbinsY, binningY); // muon
    ratesMean[1][iSource] = new TH2D(Form("%sRM_-13", source[iSource].c_str()),  Form("%s Rate Mean PDG = -13", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // antimuon
    ratesMean[2][iSource] = new TH2D(Form("%sRM_11", source[iSource].c_str()),   Form("%s Rate Mean PDG = 11", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // electron
    ratesMean[3][iSource] = new TH2D(Form("%sRM_-11", source[iSource].c_str()),  Form("%s Rate Mean PDG = -11", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // positron
    ratesMean[4][iSource] = new TH2D(Form("%sRM_22", source[iSource].c_str()),   Form("%s Rate Mean PDG = 22", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // gamma
    ratesMean[5][iSource] = new TH2D(Form("%sRM_211", source[iSource].c_str()),  Form("%s Rate Mean PDG = 211", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // pion+
    ratesMean[6][iSource] = new TH2D(Form("%sRM_-211", source[iSource].c_str()), Form("%s Rate Mean PDG = -211", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // pion-
    ratesMean[7][iSource] = new TH2D(Form("%sRM_2212", source[iSource].c_str()), Form("%s Rate Mean PDG = 2212", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // proton
    ratesMean[8][iSource] = new TH2D(Form("%sRM_2112", source[iSource].c_str()), Form("%s Rate Mean PDG = 2112", source[iSource].c_str()),
				     nbinsX, binningX, nbinsY, binningY); // neutron


    // These are the rate histograms for each toy, they are going to be filled everytime
    rates[0][iSource] = new TH2D(Form("%sR_13", source[iSource].c_str()),   Form("%sRate Error PDG = 13", source[iSource].c_str()),   
				 nbinsX, binningX, nbinsY, binningY); // muon
    rates[1][iSource] = new TH2D(Form("%sR_-13", source[iSource].c_str()),  Form("%s Rate Error PDG = -13", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // antimuon
    rates[2][iSource] = new TH2D(Form("%sR_11", source[iSource].c_str()),   Form("%s Rate Error PDG = 11", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // electron
    rates[3][iSource] = new TH2D(Form("%sR_-11", source[iSource].c_str()),  Form("%s Rate Error PDG = -11", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // positron
    rates[4][iSource] = new TH2D(Form("%sR_22", source[iSource].c_str()),   Form("%s Rate Error PDG = 22", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // gamma
    rates[5][iSource] = new TH2D(Form("%sR_211", source[iSource].c_str()),  Form("%s Rate Error PDG = 211", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // pion+
    rates[6][iSource] = new TH2D(Form("%sR_-211", source[iSource].c_str()), Form("%s Rate Error PDG = -211", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // pion-
    rates[7][iSource] = new TH2D(Form("%sR_2212", source[iSource].c_str()), Form("%s Rate Error PDG = 2212", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // proton
    rates[8][iSource] = new TH2D(Form("%sR_2112", source[iSource].c_str()), Form("%s Rate Error PDG = 2112", source[iSource].c_str()),
				 nbinsX, binningX, nbinsY, binningY); // neutron
  }

  const char *TargetString[7]  = {"CH", "O", "Al", "Cu", "Fe", "Zn", "Pb"};
  for(int iSource = 0; iSource < 4; iSource++){
    for(int iPDG = 0; iPDG < 9; iPDG++){
      for(int iTarget = 1; iTarget < 8; iTarget++){
	ratesUncertainties[iPDG][iSource]->GetXaxis()->SetBinLabel(iTarget, TargetString[iTarget-1]);
	ratesMean         [iPDG][iSource]->GetXaxis()->SetBinLabel(iTarget, TargetString[iTarget-1]);
	rates             [iPDG][iSource]->GetXaxis()->SetBinLabel(iTarget, TargetString[iTarget-1]);
      }
      ratesUncertainties[iPDG][iSource]->GetXaxis()->SetTitle("Target");
      ratesMean         [iPDG][iSource]->GetXaxis()->SetTitle("Target");
      rates             [iPDG][iSource]->GetXaxis()->SetTitle("Target");

      ratesUncertainties[iPDG][iSource]->GetYaxis()->SetTitle("Momentum [MeV]");
      ratesMean         [iPDG][iSource]->GetYaxis()->SetTitle("Momentum [MeV]");
      rates             [iPDG][iSource]->GetYaxis()->SetTitle("Momentum [MeV]");

      ratesUncertainties[iPDG][iSource]->GetXaxis()->SetDrawOption("COLZ");
      ratesMean         [iPDG][iSource]->GetXaxis()->SetDrawOption("COLZ");
      rates             [iPDG][iSource]->GetXaxis()->SetDrawOption("COLZ");

      ratesUncertainties[iPDG][iSource]->SetStats(false);
      ratesMean         [iPDG][iSource]->SetStats(false);
      rates             [iPDG][iSource]->SetStats(false);

      // ratesUncertainties[iPDG][iSource]->SetOptLogy();
      // ratesMean         [iPDG][iSource]->SetOptLogy();
      // rates             [iPDG][iSource]->SetOptLogy();
    }
  }
  // TStyle *style;
  // style->SetOptLogy();
  std::cout << "eelele" <<std::endl;
  // Define the data handler
  Data data;

  double weight[9][4];
  std::cout << "All the histos are instanciated!" << std::endl;


  for(int iFile = 0; iFile < 7; iFile++){

    input = (TTree*)files[iFile]->Get("neut");

    if(nEvents > (int)input->GetEntries())
      nEvents = (int)input->GetEntries();

    input->Print();  
    data.SetBranchesAddress(input);
    // Loop over all the entries, except for the nominal (it could be done with but it makes number look weird and I don't want to bother)
    // for(int iEntry = 0; iEntry < input->GetEntries(); iEntry++){
    for(int iEntry = 0; iEntry < nEvents; iEntry++){

      if(data.ThrowId == -1) // Remove the nominal
	continue;    
      
      input->GetEntry(iEntry);
      if(iEntry%1000 == 0) 
	std::cout << iEntry << " / " << nEvents << std::endl;
    
      for(int iPart = 0; iPart < data.nparticle; iPart++){ // Loop over all the particles exiting the nucleus
	int pdg = data.part_id[iPart];
	if(abs(pdg) == 11 || abs(pdg) == 13   ||
	   pdg == 22      || abs(pdg) == 221  ||
	   pdg == 111     || abs(pdg) == 2212 ||
	   pdg == 2112){
	  for(int iThrow = 0; iThrow < nToys; iThrow++){
	    if(data.FluxWeightByEvent [iThrow] > 100. ||
	       data.FSIWeightByEvent  [iThrow] > 100. ||
	       data.XSecWeightByEvent [iThrow] > 100. ||
	       data.TotalWeightByEvent[iThrow] > 100. ||
	       data.FluxWeightByEvent [iThrow] < 0.   ||
	       data.FSIWeightByEvent  [iThrow] < 0.   ||
	       data.XSecWeightByEvent [iThrow] < 0.   ||
	       data.TotalWeightByEvent[iThrow] < 0.)
	      continue;
	    if(pdg == 111){ // in case it is a neutral pion fill twice the gamma histogram
	      ratesMean[PDGmap[22]][0]->Fill(iFile, data.mom[iPart], 2. * data.FluxWeightByEvent [iThrow] / (double)nToys);
	      ratesMean[PDGmap[22]][1]->Fill(iFile, data.mom[iPart], 2. * data.FSIWeightByEvent  [iThrow] / (double)nToys);
	      ratesMean[PDGmap[22]][2]->Fill(iFile, data.mom[iPart], 2. * data.XSecWeightByEvent [iThrow] / (double)nToys);
	      ratesMean[PDGmap[22]][3]->Fill(iFile, data.mom[iPart], 2. * data.TotalWeightByEvent[iThrow] / (double)nToys);
	    }else{
	      ratesMean[PDGmap[pdg]][0]->Fill(iFile, data.mom[iPart], data.FluxWeightByEvent [iThrow] / (double)nToys);
	      ratesMean[PDGmap[pdg]][1]->Fill(iFile, data.mom[iPart], data.FSIWeightByEvent  [iThrow] / (double)nToys);
	      ratesMean[PDGmap[pdg]][2]->Fill(iFile, data.mom[iPart], data.XSecWeightByEvent [iThrow] / (double)nToys);
	      ratesMean[PDGmap[pdg]][3]->Fill(iFile, data.mom[iPart], data.TotalWeightByEvent[iThrow] / (double)nToys);
	    }
	  }
	}
      }
    }

    for(int iEntry = 0; iEntry < nEvents; iEntry++){
      if(data.ThrowId == -1) // Remove the nominal
	continue;    
    
      input->GetEntry(iEntry);
      if(iEntry%1000 == 0) 
	std::cout << iEntry << " / " << nEvents << std::endl;


      for(int iPart = 0; iPart < data.nparticle; iPart++){ // Loop over all the particles exiting the nucleus
	for(int iSource = 0; iSource < 4; iSource++)
	  for(int iPDG = 0; iPDG < 9; iPDG++)
	    weight[iPDG][iSource] = 0;

	int pdg = data.part_id[iPart];
	if(abs(pdg) == 11 || abs(pdg) == 13   ||
	   pdg == 22      || abs(pdg) == 221  ||
	   pdg == 111     || abs(pdg) == 2212 ||
	   pdg == 2112){
	  for(int iThrow = 0; iThrow < 1000; iThrow++){
	    if(data.FluxWeightByEvent [iThrow] > 100. ||
	       data.FSIWeightByEvent  [iThrow] > 100. ||
	       data.XSecWeightByEvent [iThrow] > 100. ||
	       data.TotalWeightByEvent[iThrow] > 100. ||
	       data.FluxWeightByEvent [iThrow] < 0.   ||
	       data.FSIWeightByEvent  [iThrow] < 0.   ||
	       data.XSecWeightByEvent [iThrow] < 0.   ||
	       data.TotalWeightByEvent[iThrow] < 0.)
	      continue;
	    if(pdg == 111){ // in case it is a neutral pion fill twice the gamma histogram
	      weight[PDGmap[22]][0] = weight[PDGmap[22]][0] + 2. * data.FluxWeightByEvent [iThrow];
	      weight[PDGmap[22]][1] = weight[PDGmap[22]][1] + 2. * data.FSIWeightByEvent  [iThrow];
	      weight[PDGmap[22]][2] = weight[PDGmap[22]][2] + 2. * data.XSecWeightByEvent [iThrow];
	      weight[PDGmap[22]][3] = weight[PDGmap[22]][3] + 2. * data.TotalWeightByEvent[iThrow];
	    }else{
	      weight[PDGmap[pdg]][0] = weight[PDGmap[pdg]][0] + data.FluxWeightByEvent [iThrow];
	      weight[PDGmap[pdg]][1] = weight[PDGmap[pdg]][1] + data.FSIWeightByEvent  [iThrow];
	      weight[PDGmap[pdg]][2] = weight[PDGmap[pdg]][2] + data.XSecWeightByEvent [iThrow];
	      weight[PDGmap[pdg]][3] = weight[PDGmap[pdg]][3] + data.TotalWeightByEvent[iThrow];
	    }
	  }
	}
      }
      for(int iPart = 0; iPart < data.nparticle; iPart++){ // Loop over all the particles exiting the nucleus
	int bin = rates[PDGmap[22]][0]->GetBin(iFile, data.mom[iPart]);
	for(int iPDG = 0; iPDG < 9; iPDG++){
	  rates[iPDG][0]->Fill(iFile, data.mom[iPart], TMath::Power(ratesMean[iPDG][0]->GetBinContent(bin) - weight[iPDG][0], 2.));
	  rates[iPDG][1]->Fill(iFile, data.mom[iPart], TMath::Power(ratesMean[iPDG][1]->GetBinContent(bin) - weight[iPDG][1], 2.));
	  rates[iPDG][2]->Fill(iFile, data.mom[iPart], TMath::Power(ratesMean[iPDG][2]->GetBinContent(bin) - weight[iPDG][2], 2.));
	  rates[iPDG][3]->Fill(iFile, data.mom[iPart], TMath::Power(ratesMean[iPDG][3]->GetBinContent(bin) - weight[iPDG][3], 2.));
	}
      }
    }
  }
  int nbinX = rates[0][0]->GetXaxis()->GetNbins();
  int nbinY = rates[0][0]->GetYaxis()->GetNbins();
  for(int iPDG = 0; iPDG < 9; iPDG++){
    for(int iBinX = 0; iBinX < nbinX+1; iBinX++){
      for(int iBinY = 0; iBinY < nbinY+1; iBinY++){
	int iBin = rates[iPDG][0]->GetBin(iBinX, iBinY);
	for(int iSource = 0; iSource < 4; iSource++){
	  if(rates[iPDG][iSource]->GetBinContent(iBin) != 0 && ratesMean[iPDG][iSource]->GetBinContent(iBin)!= 0)
	    ratesUncertainties[iPDG][iSource]->SetBinContent(iBin, TMath::Sqrt(rates[iPDG][0]->GetBinContent(iBin)) / ratesMean[iPDG][0]->GetBinContent(iBin));
	}

      }
    }
  }
 
  OutputFile->cd();
  for(int iParticle = 0; iParticle < 9; iParticle++){
    for(int iSource = 0; iSource < 4;  iSource++){
      rates[iParticle][iSource]->Write();
      ratesMean[iParticle][iSource]->Write();
      ratesUncertainties[iParticle][iSource]->Write();
    }
  }
  
  return 1;

}


// void FillChi2Cov(TH2D *rates[9][4], TH2D *ratesMean[9][4], TH2D *ratesUncertainties[9][4], TMatrixD *Correlations[9][7]){
  
//   TAxis *MomentumBins = new TAxis(nbinsY, binningY);
  
//   for(int iParticle; iParticle < 9; iParticle++){
//     for(int iTarget; iTarget < 7; iTarget++){
//       double mom_toy[3], mom_mean[3];
//       for(int iSource = 0; iSource < 4;  iSource++){
// 	for(int iBinY = 0; iBinY < nbinsY; iBinY++){ //momentum bins
// 	  int GlobalBin = rates[iParticle][iSource]->GetBin(iTarget+1, iBinY+1);
// 	  double binContentRateToy = rates[iParticle][iSource]->GetBinContent(GlobalBin);
// 	  double binContentRateMean = ratesMean[iParticle][iSource]->GetBinContent(GlobalBin);
// 	  double binContentRateUncertainties = ratesUncertainties[iParticle][iSource]->GetBinContent(GlobalBin);
// 	  if(iSource ==3){
// 	    mom_toy[iBinY] = binContentRateToy;
// 	    mom_mean[iBinY] = binContentRateMean;
// 	  }
// 	  ratesUncertainties[iParticle][iSource]->SetBinContent(GlobalBin, (binContentRateUncertainties +
// 									    TMath::Power(binContentRateToy - binContentRateMean, 2) / nToys));
// 	}
//       }
//       for(int iBinY1 = 0; iBinY1 < nbinsY; iBinY1++)
// 	for(int iBinY2 = 0; iBinY2 < nbinsY; iBinY2++)
// 	  (*Correlations[iParticle][iTarget])[iBinY1][iBinY2] += (mom_toy[iBinY1] - mom_mean[iBinY1]) * (mom_toy[iBinY2] - mom_mean[iBinY2]) / nToys;
//     }
//   }
  
  
//   return;
// }
