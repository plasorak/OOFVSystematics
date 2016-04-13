#include <string>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <vector>

#include "neutvect.h"
#include "neutpart.h"

#include "TFile.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TRandom3.h>

#include "T2KReWeight.h"
#include "T2KSyst.h"
#include "T2KNIWGUtils.h"
#include "T2KNIWGReWeight.h"
#include "T2KNeutUtils.h"
#include "T2KNeutReWeight.h"
#include "T2KJNuBeamReWeight.h"
// #include "JReWeightEnu2013a.h"
// #define JReWeightEnu jnubeam::rew::JReWeightEnu2013a
#define NIWGReWeight niwg::rew::NIWGReWeight2014a
#include "T2KWeightsStorer.h"

#include "ThrowParms.h"

#include "DataNEUT.h"

Double_t CalcMECWeight(const NeutVect *nvect, double rand[6]);
Double_t CalcFluxWeight(const NeutVect* nvect, Data& d, TVectorT<Double_t>* dial);
void CalculateOOFVTheory(const std::string infilename, const std::string outfilename);
std::string to_lower(const std::string& input);
void AddToSystematics(t2krew::T2KReWeight& rw, const t2krew::T2KSyst_t& key);
void DoThrowsAndSaveData(TTree& NEUTinputTree,
			 TVectorT<Double_t> * dial,
			 t2krew::T2KReWeight& rwXSec, std::vector<t2krew::T2KSyst_t> XSecSystematics,
			 t2krew::T2KReWeight& rwFSI,  std::vector<t2krew::T2KSyst_t> FSISystematics,
			 Data& data,  TTree* outputTree);
void SaveData(TTree& NEUTinputTree,
	      TVectorT<Double_t> * dial,
	      t2krew::T2KReWeight& rwXSec, std::vector<t2krew::T2KSyst_t> XSecSystematics,
	      t2krew::T2KReWeight& rwFSI,  std::vector<t2krew::T2KSyst_t> FSISystematics,
	      Data& data,  TTree* outputTree);
void FillNeutData(const NeutVect* nvect, Data& d);

const int nThrows = 10;
const int nEvents = 100;
const std::string BANFF_FILE_NAME = "../../data/postfit_banff_2015_data_20150417_allparams.root";

const int nFluxSystematics = 10;
const int nFSISystematics = 6;
const int nXSecSystematics = 26;

int main(int argc, char* argv[]) {
  if(argc < 2) {
    std::cout << "Usagex: ./CalculateOOFVTheory.exe <input.root> <output.root> <particle>" << std::endl;    
    std::cout << "Order matters!" << std::endl;
    std::cout << "Particle is the PDG of the particle you are interested in (by default, calculate everything)." << std::endl;
    std::cout << "You need to have the BANFF parameter in the same folder." << std::endl;
    std::cout << "This is going to make " << nThrows << " trows." << std::endl;
    return 1;
  }

  const std::string infilename = argv[1];
  const std::string outfilename = argv[2];

  CalculateOOFVTheory(infilename, outfilename);
  return 0;
}

std::string to_lower(const std::string& input) {
  std::string output = input;
  std::transform(input.begin(), input.end(), output.begin(), ::tolower);
  return output;
}

void CalculateOOFVTheory(const std::string infilename, const std::string outfilename) {
  /// --------------- Output Config --------------- ///
  TFile* out = new TFile(outfilename.c_str(), "recreate");
  TTree* out_tree = new TTree("neut", "neut");
  TTree* out_tree_nominal = new TTree("neut_nominal", "neut_nominal");
  
  Data data;
  data.SetBranches(out_tree);
  data.SetBranches(out_tree_nominal);

  std::vector<t2krew::T2KSyst_t> FSISystematics;
  FSISystematics.push_back(t2krew::kNCasc_FrPiProd_pi)  ;
  FSISystematics.push_back(t2krew::kNCasc_FrCExHigh_pi) ;
  FSISystematics.push_back(t2krew::kNCasc_FrAbs_pi)     ;
  FSISystematics.push_back(t2krew::kNCasc_FrInelLow_pi) ;
  FSISystematics.push_back(t2krew::kNCasc_FrInelHigh_pi);
  FSISystematics.push_back(t2krew::kNCasc_FrCExLow_pi)  ;

  //XSec
  std::vector<t2krew::T2KSyst_t> XSecSystematics;
  std::vector<t2krew::T2KSyst_t> XSecTune;

  // NIWG tuning (not varying)
  XSecTune.push_back(t2krew::kNXSec_VecFFCCQE);        // To be set at: 2. (NEUT reweight for RFG MaQE (Smith-Moniz))
  XSecTune.push_back(t2krew::kNIWG2014a_SF_RFG);       // To be set at: 1. (SF->RFG tuning) that was generated with SF on
  XSecTune.push_back(t2krew::kNIWG_rpaCCQE_norm);      // To be set at: 1. (nominal RPA correction)
  XSecTune.push_back(t2krew::kNIWG_rpaCCQE_shape);     // To be set at: 1. (shape RPA correction)
  // XSecTune.push_back(t2krew::kNIWG_protonFSI_bugfix);  // Proton bug fix (?)

  // To be varied		                         
  XSecSystematics.push_back(t2krew::kNXSec_MaCCQE);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_C12);
  XSecSystematics.push_back(t2krew::kNIWGMEC_Norm_C12);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_C12);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_O16);
  XSecSystematics.push_back(t2krew::kNIWGMEC_Norm_O16);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_O16);

  XSecSystematics.push_back(t2krew::kNXSec_CA5RES);
  XSecSystematics.push_back(t2krew::kNXSec_MaNFFRES);
  XSecSystematics.push_back(t2krew::kNXSec_BgSclRES);

  XSecSystematics.push_back(t2krew::kNIWG2012a_ccnueE0);
  XSecSystematics.push_back(t2krew::kNIWG2012a_dismpishp);
  XSecSystematics.push_back(t2krew::kNIWG2012a_cccohE0);
  XSecSystematics.push_back(t2krew::kNIWG2012a_cccohE0);
  XSecSystematics.push_back(t2krew::kNIWG2012a_nccohE0);
  XSecSystematics.push_back(t2krew::kNIWG2012a_ncotherE0);

  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_Al27);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_Fe56);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_Cu63);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_Zn64);
  XSecSystematics.push_back(t2krew::kNIWG2014a_pF_Pb208);
				                        
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_Al27);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_Fe56);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_Cu63);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_Zn64);
  XSecSystematics.push_back(t2krew::kNIWG2014a_Eb_Pb208);
  
  TFile *BANFF_file = new TFile(Form("%s", BANFF_FILE_NAME.c_str()), "READ");
  BANFF_file->ls();
  
  if(!BANFF_file || BANFF_file->IsZombie()){
    std::cout << "The BANFF file is Zombie, are you sure you have the correct name in the begging of the file??" << std::endl;
    return;
  }
  
  TMatrixDSym *BANFF_prefit = (TMatrixDSym*)BANFF_file->Get("prefit_cov");
  
  std::cout << "BANFF_prefit->GetNrows() " << BANFF_prefit->GetNrows() << std::endl;
  TMatrixDSym *BANFF_prefit_copy = new TMatrixDSym(nFluxSystematics + nFSISystematics + nXSecSystematics);

  TVectorD *BANFF_nominal = (TVectorD*)BANFF_file->Get("prefit_params");
  BANFF_nominal->Print();
  TVectorT<Double_t> *BANFF_nominal_copy = new TVectorT<Double_t>(nFluxSystematics + nFSISystematics + nXSecSystematics);

  // We have to play a little bit to trim the covariance matrix
  // (by the way, it is not positive definite... hum...)
  for(Int_t iParam = 0; iParam < BANFF_prefit->GetNrows(); iParam++){
    if(iParam < nFluxSystematics){
      std::cout << iParam << " vector = " << (*BANFF_nominal)[iParam] << std::endl;
      (*BANFF_nominal_copy)[iParam] = (*BANFF_nominal)[iParam];
      for(Int_t iParam2 = 0; iParam2 < nFluxSystematics; iParam2++){
	(*BANFF_prefit_copy)[iParam][iParam2] = (*BANFF_prefit)[iParam][iParam2];
	//std::cout << "     2: "<< iParam2 << "   " << (*BANFF_prefit)[iParam][iParam2] << "  " << (*BANFF_prefit_copy)[iParam][iParam2] << std::endl;
      }
    }else if(iParam >= nFluxSystematics && iParam < 100){
      continue;
    }else if(iParam >= 100 && iParam < 122){
      (*BANFF_nominal_copy)[iParam - 90] = (*BANFF_nominal)[iParam];
      
      std::cout << iParam - 90 << " vector = " << (*BANFF_nominal)[iParam] << std::endl;
      for(Int_t iParam2 = 100; iParam2 < 122; iParam2++){
	(*BANFF_prefit_copy)[iParam - 90][iParam2 - 90] = (*BANFF_prefit)[iParam][iParam2];
	//std::cout << "     2: "<< iParam2 - 75 << "   " << (*BANFF_prefit_copy)[iParam - 75][iParam2 - 75] << "  " << (*BANFF_prefit)[iParam][iParam2] << std::endl;
      }
    }
  }

  // My own prescriptions for all the heavy-targets-specific dials:
  // A-scaled from Carbon (uncorrelated), hopefully this doesn't become crazy for Pb
  for(int i = 0; i < 10; i++){
    (*BANFF_nominal_copy)[i + 122 - 90] = 1.;
  }

  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(27./12.,  2.) * (*BANFF_nominal)[107];   // kNIWG2014a_pF_Al27
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(56./12.,  2.) * (*BANFF_nominal)[107];   // kNIWG2014a_pF_Fe56
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(63./12.,  2.) * (*BANFF_nominal)[107];   // kNIWG2014a_pF_Cu63
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(64./12.,  2.) * (*BANFF_nominal)[107];   // kNIWG2014a_pF_Zn64
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(208./12., 2.) * (*BANFF_nominal)[107];   // kNIWG2014a_pF_Pb208
  
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(27./12.,  2.) * (*BANFF_nominal)[109];   // kNIWG2014a_Eb_Al27
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(56./12.,  2.) * (*BANFF_nominal)[109];   // kNIWG2014a_Eb_Fe56
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(63./12.,  2.) * (*BANFF_nominal)[109];   // kNIWG2014a_Eb_Cu63
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(64./12.,  2.) * (*BANFF_nominal)[109];   // kNIWG2014a_Eb_Zn64
  (*BANFF_prefit_copy)[i + 122 - 90][i + 122 - 90] = TMath::Power(208./12., 2.) * (*BANFF_nominal)[109];   // kNIWG2014a_Eb_Pb208

  BANFF_nominal_copy->Print();
  BANFF_prefit_copy->Print();

  data.BANFFNominal = BANFF_nominal_copy;
  std::cout << "PRINTING NOMINAL" << std::endl;
  data.BANFFNominal->Print();
  ThrowParms *throwParams = new ThrowParms(*BANFF_nominal_copy, *BANFF_prefit_copy);
  throwParams->SetSeed(63432946);
  data.throwParams = throwParams;
  data.BANFFMatrix = BANFF_prefit_copy;
  data.FluxAxis.clear();
  data.FluxAxis.push_back((TAxis*)BANFF_file->Get("nd5_numode_numu_bins"));
  data.FluxAxis.push_back((TAxis*)BANFF_file->Get("nd5_numode_numub_bins"));
  data.FluxAxis.push_back((TAxis*)BANFF_file->Get("nd5_numode_nue_bins"));
  data.FluxAxis.push_back((TAxis*)BANFF_file->Get("nd5_numode_nueb_bins"));
  
  
  // -----------------------
  // Set the nominal first
  // -----------------------
  
  TVectorT<Double_t> *DialVector = new TVectorT<Double_t>(nFluxSystematics + nFSISystematics + nXSecSystematics);
  
  for(int iSyst = 0; iSyst < nFluxSystematics; iSyst++){
    std::cout << iSyst << " Filling Flux" << std::endl;
    (*DialVector)[iSyst] = 1;
  }

  t2krew::T2KReWeight rwFSI;
  rwFSI.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  rwFSI.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());
  
  for(int iSyst = 0; iSyst < nFSISystematics; iSyst++){
    rwFSI.Systematics().Include(FSISystematics[iSyst]);
  }
  for(int iSyst = 0; iSyst < nFSISystematics; iSyst++){
    rwFSI.Systematics().SetAbsTwk(FSISystematics[iSyst]);
  }
  for(int iSyst = 0; iSyst < nFSISystematics; iSyst++){
    rwFSI.Systematics().SetTwkDial(FSISystematics[iSyst], (*BANFF_nominal_copy)(iSyst+nFluxSystematics));
    (*DialVector)[iSyst+nFluxSystematics] = (*BANFF_nominal_copy)(iSyst+nFluxSystematics);
    std::cout << iSyst+nFluxSystematics << " Filling FSI" << std::endl;
  }
  rwFSI.Reconfigure();


  t2krew::T2KReWeight rwXSec;
  rwXSec.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  rwXSec.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());

  // Define the dials that are constant (NIWG tunning for the first 5)
  for(int iTune = 0; iTune < XSecTune.size(); iTune++){
    rwXSec.Systematics().Include(XSecTune[iTune]);
    rwXSec.Systematics().SetAbsTwk(XSecTune[iTune]);
  }

  // Define the real systematics (those which will actually be varied)
  for(int iSyst = 0; iSyst < nXSecSystematics; iSyst++){
    rwXSec.Systematics().Include(XSecSystematics[iSyst]);
    rwXSec.Systematics().SetAbsTwk(XSecSystematics[iSyst]);
  }

  // Actually set them
  rwXSec.Systematics().SetTwkDial(XSecTune[0],  2.); // NXSec_VecFFCCQE;   
  rwXSec.Systematics().SetTwkDial(XSecTune[1],  1.); // NIWG2014a_SF_RFG;  
  rwXSec.Systematics().SetTwkDial(XSecTune[2],  1.); // NIWG_rpaCCQE_norm; 
  rwXSec.Systematics().SetTwkDial(XSecTune[3], -1.); // NIWG_rpaCCQE_shape;

  t2krew::T2KSyst sys;
  for(int iSyst = 0; iSyst < nXSecSystematics; iSyst++){
    // rwXSec.Systematics().SetTwkDial(XSecSystematics[iSyst],  (*BANFF_nominal_copy)(iSyst + nFluxSystematics + nFSISystematics));
    // if(iSyst + nFluxSystematics + nFSISystematics >= 23 && iSyst + nFluxSystematics + nFSISystematics <= 25)
    rwXSec.Systematics().SetTwkDial(XSecSystematics[iSyst], 0.);
    (*DialVector)[iSyst+ nFSISystematics + nFluxSystematics] = (*BANFF_nominal_copy)(iSyst + nFluxSystematics + nFSISystematics);
    std::cout << iSyst+ nFSISystematics + nFluxSystematics << " Filling XSec --> " <<  (*BANFF_nominal_copy)(iSyst + nFluxSystematics + nFSISystematics) << " name : " << sys.AsString(XSecSystematics[iSyst]) << std::endl;
    
  }
  rwXSec.Reconfigure();
  
  // Set the neut tree
  TFile* input = new TFile(infilename.c_str(), "READ");
  TTree& neuttree = *(static_cast<TTree*>(input->Get("neuttree")));
  neuttree.LoadTree(0);

  // Save in the output tree (entry 0 is the nominal)
  data.ThrowId = 0; // this is the enty of the weight array so it has to be >= 0 and smaller that data.nThrows otherwise the tree is unreadable (and that sucks)
  data.nThrows = 1; // one throw, the nominal

  
  // -----------------------
  // Save the nominal first
  // -----------------------
  SaveData(neuttree, DialVector, rwXSec, XSecSystematics, rwFSI, FSISystematics, data, out_tree_nominal);
  out_tree_nominal->Fill();

  // -----------------------
  // Now do the throws
  // -----------------------
  data.nThrows = nThrows;
  DoThrowsAndSaveData(neuttree, DialVector, rwXSec, XSecSystematics, rwFSI, FSISystematics, data, out_tree);
  
  out->cd();
  out_tree->Write();
  out_tree_nominal->Write();
  input->Close();
    
  return;
}


void AddToSystematics(t2krew::T2KReWeight& rw, const t2krew::T2KSyst_t& key) {
  rw.Systematics().Include(key);
  rw.Systematics().SetAbsTwk(key);
  return;
}


void DoThrowsAndSaveData(TTree& NEUTinputTree,
			 TVectorT<Double_t> * dial,
			 t2krew::T2KReWeight& rwXSec, std::vector<t2krew::T2KSyst_t> XSecSystematics,
			 t2krew::T2KReWeight& rwFSI,  std::vector<t2krew::T2KSyst_t> FSISystematics,
			 Data& data,  TTree* outputTree){

  // NEUT Input
  NeutVect* neutvect = new NeutVect();

  NEUTinputTree.SetBranchStatus("vectorbranch", true);
  NEUTinputTree.SetBranchAddress("vectorbranch", &neutvect);
  NEUTinputTree.SetBranchStatus("vertexbranch", false);
  
  int nEventsGenerated = NEUTinputTree.GetEntries();
 
  if(nEventsGenerated > nEvents)
    nEventsGenerated = nEvents;
  if(data.ThrowId < 0 || data.ThrowId > data.nThrows){
    std::cout << "There is a problem in the value of data.ThrowId: " << data.ThrowId << std::endl;
    std::cout << "It should be >= 0 and < " << data.nThrows << std::endl;
    exit(1);
  }
  data.flux = 1. / (double)nEventsGenerated;

  for(int iThrow = 0; iThrow < nThrows; iThrow++){      
    std::cout << "throw number " << iThrow << " / " << nThrows << std::endl;
    data.ThrowId = iThrow;
    std::vector<double> par_throw;
    data.throwParams->ThrowSet(par_throw);

    // We need lots of random numbers for the MEC parameters for every target
    // We are throwing them totally uncorrelated
    // this is a simple hack, but if I understood well it should be the equivalent to the C_MEC and O_MEC parameters for the heavy targets.
    double rand[6];
    double z[2];
    TRandom3* rd = new TRandom3(iThrow * 1947207);
    for(int iRand = 0; iRand < 3; iRand++){  //Stolen from T2KReWeight's ThrowParms.cxx
      double u = 2. * (rd->Rndm()) - 1.;
      double v = 2. * (rd->Rndm()) - 1.;
      double s = u*u+v*v;
      
      while(s==0 || s>=1.){
	u = 2. * (rd->Rndm()) - 1.;
	v = 2. * (rd->Rndm()) - 1.;
	s = u*u+v*v;
      }
      
      z[0] = u * sqrt(-2.*TMath::Log(s)/s);
      z[1] = v * sqrt(-2.*TMath::Log(s)/s);
      rand[iRand] = z[0];
      rand[iRand+3] = z[1];
    }
    
    std::cout << rand[0] << std::endl;

    for(int j = 0; j < (Int_t) par_throw.size(); j++){
      std::cout << j << "   " << par_throw.at(j) << std::endl;
      (*dial)[j] = par_throw.at(j);
    }
    // int i = 16;
    // std::cout << "(*dial_nom)[" << i << "] " << (*dial)[i] << std::endl;
    // std::cout << "(*data.BANFFMatrix)["<< i << "][" << i << "]  " << (*data.BANFFMatrix)[i][i] << std::endl;
    // std::cout << "(*data.BANFFNominal)[" << i << "] " << (*data.BANFFNominal)[i] << std::endl;

    // (*dial)[i] = (*data.BANFFMatrix)[i][i] + (*data.BANFFNominal)[i];
    // std::cout << "(*dial_thr)[" << i << "] " << (*dial)[i] << std::endl;

    // This would work if I'd used neutgeom to generate events, however I'm quite lazy and didn't bother to learn how to do it
      
    // // Set the dials
    // for(int iSyst = 0;  iSyst < nFluxSystematics; iSyst++){
    //   rwFlux.Systematics().SetTwkDial(FluxSystematics[iSyst], (*DialVector)(iSyst)); 
    // }
    // rwFlux.Reconfigure();
      
    for(int iSyst = 0;  iSyst < nFSISystematics; iSyst++){
      rwFSI.Systematics().SetTwkDial(FSISystematics[iSyst], (*dial)(iSyst+nFluxSystematics));
    }
    rwFSI.Reconfigure();
      
    for(int iSyst = 0;  iSyst < nXSecSystematics; iSyst++){
      //rwXSec.Systematics().SetTwkDial(XSecSystematics[iSyst], (*dial)(iSyst+nFluxSystematics+nFSISystematics));
      rwXSec.Systematics().SetTwkDial(XSecSystematics[iSyst], (*dial)(iSyst+nFluxSystematics+nFSISystematics) - 1.);
    }
    
    rwXSec.Reconfigure();

    for (int iEvent = 0; iEvent < nEventsGenerated; ++iEvent) {
      NEUTinputTree.GetEntry(iEvent);
      data.FluxWeight  [iEvent][data.ThrowId] = CalcFluxWeight(neutvect, data, dial);
      data.XSecWeight  [iEvent][data.ThrowId] = rwXSec.CalcWeight(neutvect) * CalcMECWeight(neutvect, rand);
      data.FSIWeight   [iEvent][data.ThrowId] = rwFSI.CalcWeight(neutvect);
      data.TotalWeight [iEvent][data.ThrowId] = data.FluxWeight[iEvent][data.ThrowId] * data.XSecWeight[iEvent][data.ThrowId] * data.FSIWeight[iEvent][data.ThrowId];
      // Save the weight, this is where the calculation is made  (takes a long time)
    }
    
  }
  
  std::cout << "Reorganizing...." << std::endl;
  for (int iEvent = 0; iEvent < nEventsGenerated; ++iEvent) {
    NEUTinputTree.GetEntry(iEvent);
    data.FillNeutData(neutvect);
    for(int iThrow2 = 0; iThrow2 < nThrows; iThrow2++){
      data.FluxWeightByEvent  [iThrow2] = data.FluxWeight  [iEvent][iThrow2];
      data.XSecWeightByEvent  [iThrow2] = data.XSecWeight  [iEvent][iThrow2];
      data.FSIWeightByEvent   [iThrow2] = data.FSIWeight   [iEvent][iThrow2];
      data.TotalWeightByEvent [iThrow2] = data.TotalWeight [iEvent][iThrow2];  
    }
    // Save the tree for the first pass
    outputTree->Fill();
  }

}


void SaveData(TTree& NEUTinputTree,
	      TVectorT<Double_t> * dial,
	      t2krew::T2KReWeight& rwXSec, std::vector<t2krew::T2KSyst_t> XSecSystematics,
	      t2krew::T2KReWeight& rwFSI,  std::vector<t2krew::T2KSyst_t> FSISystematics,
	      Data& data,  TTree* outputTree){

  // NEUT Input
  NeutVect* neutvect = new NeutVect();

  NEUTinputTree.SetBranchStatus("vectorbranch", true);
  NEUTinputTree.SetBranchAddress("vectorbranch", &neutvect);
  NEUTinputTree.SetBranchStatus("vertexbranch", false);
  
  int nEventsGenerated = NEUTinputTree.GetEntries();
 
  if(nEventsGenerated > nEvents)
    nEventsGenerated = nEvents;
  if(data.ThrowId < 0 || data.ThrowId > data.nThrows){
    std::cout << "There is a problem in the value of data.ThrowId: " << data.ThrowId << std::endl;
    std::cout << "It should be >= 0 and < " << data.nThrows << std::endl;
    exit(1);
  }
 
  data.flux = 1. / (double)nEventsGenerated;
  for (int iEvent = 0; iEvent < nEventsGenerated; ++iEvent) {
    NEUTinputTree.GetEntry(iEvent);
    
    // Save the tree
    data.FillNeutData(neutvect);
    
    data.FluxWeightByEvent  [data.ThrowId] = CalcFluxWeight(neutvect, data, dial);
    data.XSecWeightByEvent  [data.ThrowId] = rwXSec.CalcWeight(neutvect);
    data.FSIWeightByEvent   [data.ThrowId] = rwFSI.CalcWeight(neutvect);
    data.TotalWeightByEvent [data.ThrowId] = data.FluxWeightByEvent[data.ThrowId] * data.XSecWeightByEvent[data.ThrowId] * data.FSIWeightByEvent[data.ThrowId];
  
    outputTree->Fill();
  }
  
}


Double_t CalcMECWeight(const NeutVect *nvect, double rand[6]){

  if(abs(nvect->Mode) != 2)
    return 1.;

  if(nvect->TargetZ == 6 || nvect->TargetZ == 8)
    return 1.;
  else if(nvect->TargetZ == 30) //Zinc
    return rand[0] + 1;
  else if(nvect->TargetZ == 29) // Copper
    return rand[1] + 1;
  else if(nvect->TargetZ == 13) // Al
    return rand[2] + 1;
  else if(nvect->TargetZ == 26) //Fe
    return rand[3] + 1;
  else if(nvect->TargetZ == 82) // Pb
    return rand[4] + 1;
  else{
    std::cout << "Problem in the target!!" << std::endl;
    return rand[5] + 1;
  }

}

Double_t CalcFluxWeight(const NeutVect* nvect, Data& d, TVectorT<Double_t>* dial){
  
  const NeutPart& neutrino = *(const_cast<NeutVect*>(nvect)->PartInfo(0));
  const TLorentzVector& incoming_neutrino_momentum = neutrino.fP;
  d.enu = incoming_neutrino_momentum.E();
  if(neutrino.fPID == 14){
    return (*dial)[d.FluxAxis[0]->FindBin(d.enu/1000.) - 1];
  }
  else if(neutrino.fPID == -14){
    return (*dial)[d.FluxAxis[1]->FindBin(d.enu/1000.) + d.FluxAxis[0]->GetNbins() - 3];
  }
  else if(neutrino.fPID == 12){
    return (*dial)[d.FluxAxis[2]->FindBin(d.enu/1000.) + d.FluxAxis[0]->GetNbins() + d.FluxAxis[1]->GetNbins() - 5];
  }
  else if(neutrino.fPID == -12){
    return (*dial)[d.FluxAxis[3]->FindBin(d.enu/1000.) + d.FluxAxis[0]->GetNbins() + d.FluxAxis[1]->GetNbins() + d.FluxAxis[2]->GetNbins() - 7];
  }
  else{
    std::cout << "error in the flux!" << std::endl;
    exit(1);
  }
    
}

