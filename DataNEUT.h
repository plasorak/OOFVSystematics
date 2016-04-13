// #include <vector>

// #include <TMatrixDSym.h>
// #include <TLorentzVector.h>
// #include <TVectorT.h>

// #include "T2KReWeight.h"
// #include "T2KSyst.h"
// #include "ThrowParms.h"

// #include "neutvect.h"
// #include "neutpart.h"



class Data {
 public:
  // static const int NMAXPARTICLE = 100;
  // static const int NMAXTHROWS = 1000;
  // static const int NMAXEVENTS = 100000;
  int nEvents;
  int EventNo;
  int ThrowId;
  
  // Some general stuff for cross section
  double flux;
  double xsec;
  double enu;
  double pnu_x;
  double pnu_y;
  double pnu_z;
  double elep;
  double plep;
  double plep_x;
  double plep_y;
  double plep_z;
  double costh;
  double theta;
  double eccqe;
  double q0;
  double q3;

  int TargetA;
  int TargetZ;
  int TargetH;
  int Ibound; 
  int VNuclIni;
  int VNuclFin;
  int PFSurf;
  int PFMax;
  int IntMode;

  int nparticle;
  int part_id[100];
  double mom[100];
  double mom_x[100];
  double mom_y[100];
  double mom_z[100];
  
  double t2krw_RPA;
  double t2krw_RelRPA;

  double liberated_nucleon_mass;
  double target_nucleon_mass;
  double mlep;
  double mnu;
  double pnu;

  // ThrowParms *throwParams;
  int nThrows;
  // double FluxWeight [100000][1000];
  // double FSIWeight  [100000][1000];
  // double XSecWeight [100000][1000];
  // double TotalWeight[100000][1000];

  double FluxWeightByEvent[1000];
  double FSIWeightByEvent[1000];
  double XSecWeightByEvent[1000];
  double TotalWeightByEvent[1000];

  // TVectorT<Double_t> *BANFFNominal;
  // TMatrixDSym *BANFFMatrix;
  // std::vector<TAxis*> FluxAxis;

  void SetBranchesAddress(TTree* input){
    input->SetBranchAddress("ThrowId", &this->ThrowId);
    input->SetBranchAddress("EventNo", &this->EventNo);
    
    input->SetBranchAddress("nparticle",   &this->nparticle);
    input->SetBranchAddress("part_id",      this->part_id);
    input->SetBranchAddress("mom",          this->mom);
    input->SetBranchAddress("mom_x",        this->mom_x);
    input->SetBranchAddress("mom_y",        this->mom_y);
    input->SetBranchAddress("mom_z",        this->mom_z);
    input->SetBranchAddress("nThrows",     &this->nThrows);
    input->SetBranchAddress("FluxWeight",   this->FluxWeightByEvent);
    input->SetBranchAddress("FSIWeight",    this->FSIWeightByEvent);
    input->SetBranchAddress("XSecWeight",   this->XSecWeightByEvent);
    input->SetBranchAddress("TotalWeight",  this->TotalWeightByEvent);

    return;
  }


  // void SetBranches(TTree *out_tree){  
  //   out_tree->Branch("nEvents", &this->nEvents, "nEvents/I");
  //   out_tree->Branch("ThrowId", &this->ThrowId, "ThrowId/I");
  //   out_tree->Branch("EventNo", &this->EventNo, "EventNo/I");
  //   out_tree->Branch("flux",    &this->flux,    "flux/D");
  //   out_tree->Branch("xsec",    &this->xsec,    "xsec/D");
  //   out_tree->Branch("enu",     &this->enu,     "enu/D");
  //   out_tree->Branch("pnu_x",   &this->pnu_x,   "pnu_x/D");
  //   out_tree->Branch("pnu_y",   &this->pnu_y,   "pnu_y/D");
  //   out_tree->Branch("pnu_z",   &this->pnu_z,   "pnu_z/D");
  //   out_tree->Branch("elep",    &this->elep,    "elep/D");
  //   out_tree->Branch("plep",    &this->plep,    "plep/D");
  //   out_tree->Branch("plep_x",  &this->plep_x,  "plep_x/D");
  //   out_tree->Branch("plep_y",  &this->plep_y,  "plep_y/D");
  //   out_tree->Branch("plep_z",  &this->plep_z,  "plep_z/D");
  //   out_tree->Branch("costh",   &this->costh,   "costh/D");
  //   out_tree->Branch("theta",   &this->theta,   "theta/D");
  //   out_tree->Branch("eccqe",   &this->eccqe,   "eccqe/D");
  //   out_tree->Branch("q0",      &this->q0,      "q0/D");
  //   out_tree->Branch("q3",      &this->q3,      "q3/D");
  //   out_tree->Branch("IntMode", &this->IntMode, "IntMode/I");
  //   out_tree->Branch("nparticle", &this->nparticle, "nparticle/I");
  //   out_tree->Branch("part_id",    this->part_id,   "part_id[nparticle]/I");
  //   out_tree->Branch("mom",        this->mom,       "mom[nparticle]/D");
  //   out_tree->Branch("mom_x",      this->mom_x,     "mom_x[nparticle]/D");
  //   out_tree->Branch("mom_y",      this->mom_y,     "mom_y[nparticle]/D");
  //   out_tree->Branch("mom_z",      this->mom_z,     "mom_z[nparticle]/D");
    
  //   out_tree->Branch("nThrows",     &this->nThrows,     "nThrows/I");
  //   out_tree->Branch("FluxWeight",   this->FluxWeightByEvent,  "FluxWeight[nThrows]/D");
  //   out_tree->Branch("FSIWeight",    this->FSIWeightByEvent,   "FSIWeight[nThrows]/D");
  //   out_tree->Branch("XSecWeight",   this->XSecWeightByEvent,  "XSecWeight[nThrows]/D");
  //   out_tree->Branch("TotalWeight",  this->TotalWeightByEvent, "TotalWeight[nThrows]/D");
    
  //   return;
  // }
  
  // void FillNeutData(const NeutVect* nvect){
  //   // In here we save the data from the neut format in the data class
  //   this->EventNo  = nvect->EventNo;
  //   this->TargetA  = nvect->TargetA;
  //   this->TargetZ  = nvect->TargetZ;
  //   this->TargetH  = nvect->TargetH;
  //   this->Ibound   = nvect->Ibound;
  //   this->VNuclIni = nvect->VNuclIni;
  //   this->VNuclFin = nvect->VNuclFin;
  //   this->PFSurf   = nvect->PFSurf;
  //   this->PFMax    = nvect->PFMax;
  //   this->xsec     = nvect->Totcrs;
  //   this->IntMode  = nvect->Mode;

  //   const NeutPart& neutrino = *(const_cast<NeutVect*>(nvect)->PartInfo(0));
  //   const TLorentzVector& incoming_neutrino_momentum = neutrino.fP;
  //   this->enu = incoming_neutrino_momentum.E();
  //   this->mnu = neutrino.fMass;
  //   this->pnu = std::sqrt(this->enu*this->enu - this->mnu*this->mnu);
  //   this->pnu_x = incoming_neutrino_momentum.Px();
  //   this->pnu_y = incoming_neutrino_momentum.Py();
  //   this->pnu_z = incoming_neutrino_momentum.Pz();
    
  //   const NeutPart& target_nucleon = *(const_cast<NeutVect*>(nvect)->PartInfo(1));
  //   const int target_nucleon_PID = target_nucleon.fPID;
  //   this->target_nucleon_mass = target_nucleon.fMass;
    
  //   const int parts = nvect->Npart();
  //   bool outgoing_lepton_found = false;
  //   bool outgoing_nucleon_found = false;
    
  //   this->nparticle = parts - 2; // remove the neutrino and the struck nucleon
    
  //   for (int k=2; k<parts; ++k) {
  //     const NeutPart& kpart = *(const_cast<NeutVect*>(nvect)->PartInfo(k));
      
  //     bool k_is_charged_lepton = false;
  //     bool k_is_nucleon = false;
      
  //     this->part_id[k-2] = kpart.fPID;
  //     const TLorentzVector& four_momentum = kpart.fP;
  //     Double_t Energy = four_momentum.E();
  //     Double_t Mass   = kpart.fMass;
  //     this->mom[k-2]     = TMath::Sqrt(Energy * Energy - Mass * Mass);
  //     this->mom_x[k-2]   = four_momentum.Px();
  //     this->mom_y[k-2]   = four_momentum.Py();
  //     this->mom_z[k-2]   = four_momentum.Pz();
      
  //     switch (kpart.fPID) {
  //     case -15: //< tau^+
  //     case -13: //< mu^+
  //     case -11: //< e^+
  //     case  11: //< e^-
  //     case  13: //< mu^-
  //     case  15: //< tau^-
  // 	k_is_charged_lepton = true;
  // 	break;
  //     case 2112: //< neutron
  //     case 2212: //< proton
  // 	k_is_nucleon = true;
  // 	break;
  //     default:
  // 	break;
  //     }
      
  //     if (k_is_charged_lepton) {
  // 	if (outgoing_lepton_found) { continue; }
	
  // 	// We've found the outgoing lepton!
  // 	const TLorentzVector& four_momentum = kpart.fP;
  // 	this->elep = four_momentum.E();
  // 	this->mlep = kpart.fMass;
  // 	this->plep = std::sqrt(this->elep*this->elep - this->mlep*this->mlep);
  // 	this->theta = four_momentum.Angle(incoming_neutrino_momentum.Vect());
  // 	this->costh = std::cos(this->theta);
	
  // 	this->plep_x = four_momentum.Px();
  // 	this->plep_y = four_momentum.Py();
  // 	this->plep_z = four_momentum.Pz();
	
  // 	this->eccqe = calc_q0();
  // 	this->q0 = calc_q3();
  // 	this->q3 = calc_ereco();
	
  // 	// We only want the first lepton, forget about tertiary processes.
  // 	outgoing_lepton_found = true;
	
  //     }
      
  //     if(k_is_nucleon) {
  // 	if(outgoing_nucleon_found) { continue; }
  // 	if(kpart.fPID == target_nucleon_PID) { continue; }
  // 	this->liberated_nucleon_mass = kpart.fMass;
	
  // 	// We only want the first nucleon, forget about tertiary processes.
  // 	outgoing_nucleon_found = true;
  //     }
  //   }
  // }
  
  
  // This we don't care
  
  // double calc_q0() {
  //   double retval = (this->enu - this->elep);
  //   return retval;
  // }
  
  // double calc_q3() {
  //   double retval = std::sqrt(  (this->pnu_x-this->plep_x)*(this->pnu_x-this->plep_x)
  // 				+ (this->pnu_y-this->plep_y)*(this->pnu_y-this->plep_y)
  // 				+ (this->pnu_z-this->plep_z)*(this->pnu_z-this->plep_z)
  // 				);
  //   return retval;
  // }
  
  // double calc_ereco() {
  //   const double MeV_per_GeV = 1.e3;
  //   const double binding_energy = (this->VNuclFin - this->VNuclIni)*MeV_per_GeV;
    
  //   const double mp_sq = this->liberated_nucleon_mass * this->liberated_nucleon_mass;
  //   const double mn_minus_Eb = (this->target_nucleon_mass - binding_energy);
  //   const double mn_minus_Eb_sq = mn_minus_Eb*mn_minus_Eb;
  //   const double ml_sq = this->mlep*this->mlep;
    
  //   const double numerator = mp_sq - mn_minus_Eb_sq - ml_sq + 2.*mn_minus_Eb*this->elep;
  //   const double denominator = 2.*(mn_minus_Eb - this->elep + this->plep*this->costh);
    
  //   const double retval = (numerator / denominator);
  //   return retval;
  // }
  
};
