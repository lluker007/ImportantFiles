/*
Author: Aritra Bal
13/08/2021

This converts ATLAS 2016 data to CMS nanoAODplus ntuple format.

The following branches have been added and given dummy values, not dependent on ATLAS data. 
Muon_sip3d[] -> set to 0.0
Muon_isPFcand[] -> set to true
Muon_isGlobal[] -> set to true

Electron_isEE[] -> set to true
Electron_isEB[] -> set to false

Also note that relative isolations have been normalised to pt and scaled to the square of the radii. 

To be used primarily for the Higgs->4L analysis.
*/
/*
Edited by Louise Luker
17/07/2025

The edits to the original file are to update the code for use with the converted PHYSLITE files

The original file can be found:
/afs/desy.de/user/g/geiser/public/summerstudents/Aritra/public/ATLAS2016_to_CMS/atlas2016_to_cms_nanoAODplus_revived.C
*/ 



// system include files
#include <memory>
#include <iostream>

// for histogramming
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TTree.h" //
#include "TMath.h"
#include "TFile.h"
//#include "Math/VectorUtil_cint.h"

#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>

#include "TMatrixT.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
//#include "Math/VectorUtil.h"
//#include "Math/GenVector/VectorUtil.h"

#include <cmath> //apparenetly need to include this

const int pdgmuon = 13;
const int pdgelec = 11;
const int pdgtau = 15;



using namespace std;

void atlastocmsfile (const std::string& path,const std::string& outfile, int year){

gROOT->Reset();
gStyle->SetOptStat("nemruo");
gSystem->Load("libPhysics"); //had to add because of errors //edit
gSystem->Load("libMathCore"); //had to add because of errors //edit
 
//inputing the data & declaring histograms from ATLAS DATA

TFile *f = new TFile(path.c_str());
TTree *ntree = (TTree*)f->Get("analysis"); //will probably need changing as the trees in the root file do not have the same name ///changed from "mini" to "analysis" //edit
vector<Float_t> *lep_pt=0, *lep_eta=0, *lep_phi=0, *lep_E=0,*lep_ptcone30=0,*lep_etcone20=0,*lep_z0=0,*lep_trackd0pvunbiased=0,*lep_tracksigd0pvunbiased=0;
vector<UInt_t> *lep_type=0,*lep_charge=0,*lep_flag=0;
Int_t lep_n; //need to change UInt_t to Int_t because of how it is stored in the root file //edit
Float_t mass_dimuon;
//Int_t       runNumber,eventNumber,channelNumber,pvxp_n; //this line in comments as is original code but does not seem to work for our input file //edit
UInt_t runNumber,channelNumber; //added because of errors //edit
ULong64_t eventNumber; //added because of errors //edit
float       sqm1, s1, s2, s, w,mcWeight,vxp_z;
float       px_mu1, py_mu1, pz_mu1, p_mu1;
float       px_mu2, py_mu2, pz_mu2, p_mu2;
Bool_t      trigE, trigM, passGRL,hasGoodVertex;
vector<Bool_t> *lep_truthMatched=0,*lep_trigMatched=0;   

/*
ntree->SetBranchStatus("*", 0);
sqm1 = (0.105658) * (0.105658);
ntree->SetBranchStatus("lep_n", 1);  
ntree->SetBranchStatus("lep_pt", 1); 
ntree->SetBranchStatus("lep_eta", 1);
ntree->SetBranchStatus("lep_phi", 1);
ntree->SetBranchStatus("lep_E", 1);
ntree->SetBranchStatus("lep_charge", 1);
*/


ntree->SetBranchAddress("lep_n", &lep_n);  
ntree->SetBranchAddress("lep_pt", &lep_pt);
ntree->SetBranchAddress("lep_eta", &lep_eta);
ntree->SetBranchAddress("lep_phi", &lep_phi);
ntree->SetBranchAddress("lep_e", &lep_E); //change ///changed "lep_E" to "lep_e" //edit
ntree->SetBranchAddress("lep_type",&lep_type);
ntree->SetBranchAddress("lep_charge", &lep_charge);
ntree->SetBranchAddress("lep_z0", &lep_z0);
ntree->SetBranchAddress("runNumber",&runNumber);
ntree->SetBranchAddress("eventNumber",&eventNumber);
ntree->SetBranchAddress("channelNumber",&channelNumber);
ntree->SetBranchAddress("mcWeight",&mcWeight);
ntree->SetBranchAddress("lep_ptvarcone30",&lep_ptcone30); //change ///changed "lep_ptcone30" to "lep_ptvarcone30" //edit
ntree->SetBranchAddress("lep_topoetcone20",&lep_etcone20); //change ///changed "lep_etcone20" to "lep_topoetcone20" //edit
ntree->SetBranchAddress("lep_d0sig",&lep_tracksigd0pvunbiased); //change (i think) ///changed "lep_tracksigd0pvunbiased" to "lep_d0sig" //edit
ntree->SetBranchAddress("lep_d0",&lep_trackd0pvunbiased); //change (i think) ///changed "lep_trackd0pvunbiased" to "lep_d0" //edit
//ntree->SetBranchAddress("lep_truthMatched",&lep_truthMatched); //not sure, doesn't exist ///commented out because it throws an error due to not existing in the ROOT file //edit
ntree->SetBranchAddress("lep_isTrigMatched",&lep_trigMatched); //change (i think) ///changed "lep_trigMatched" to "lep_isTrigMatched" //edit
ntree->SetBranchAddress("trigE",&trigE);
ntree->SetBranchAddress("trigM",&trigM);
// branches present only in 2012 data but not in 2016 data 
///edit - will probably not need at all ///commented it out because its throwing an error >:( //edit
/*
if (year==2012)
{
    ntree->SetBranchAddress("pvxp_n",&pvxp_n);
    ntree->SetBranchAddress("vxp_z",&vxp_z);
    ntree->SetBranchAddress("passGRL",&passGRL);
    ntree->SetBranchAddress("hasGoodVertex",&hasGoodVertex);
    ntree->SetBranchAddress("lep_flag",&lep_flag);  
}
*/

cout << "Input root file complete" << endl;

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Muon Collection Start ////////////////////////////
///////////////////////////////////////////////////////////////////////////////


Float_t mu_E = -9999.;
Float_t mu_relPFIsoR03 = -9999.;
Float_t mu_relPFIsoR04 = -9999.;

Float_t mu_ip3d = -9999.;
Float_t mu_Errip3d = -9999.;
Float_t mu_sip3d = -9999.;
Float_t mu_M = 0.;


Float_t el_E = -9999.;
Float_t el_relPFIsoR03 = -9999.;
Float_t el_relPFIsoR04 = -9999.;

Float_t el_ip3d = -9999.;
Float_t el_Errip3d = -9999.;
Float_t el_sip3d = -9999.;
Float_t el_M = 0.;

  // tlorentzvector
TLorentzVector p4Mu;  
p4Mu.SetPtEtaPhiE(0., 0., 0., 0.);

// unique_ptr<TTree> t_event;
 
//declaration of output file and tree
TFile *file = new TFile(outfile.c_str(), "RECREATE");
TTree *t_event = new TTree("Events", "Events");   // renamed from Event to Events 
//t_event = unique_ptr<TTree>(new TTree("Event", "Event"));
//t_event->SetAutoSave(500000000);

//--------------- declare variables you want to put into the tree ------------//
UInt_t run;
ULong64_t event;
UInt_t luminosityBlock;
Float_t Generator_weight;
//Int_t b4_nMuon;

Int_t PV_npvs;
Float_t PV_z;

//---------------------------- for Muon variables ----------------------------//

// variables // * are filled with nonsense atm
// because some variables are not available in old CMSSW eg. softId
// and some of it I don't fill it properly yet eg. pdgId
/* 
vector<Int_t> Muon_charge;
vector<Float_t> Muon_dxy;
vector<Float_t> Muon_dxyErr;
vector<Float_t> Muon_dz;

vector<Int_t> Muon_genPartIdx; 
vector<Float_t> Muon_pfRelIso03_all;
vector<Float_t> Muon_pfRelIso03_chg;

vector<Float_t> Muon_phi;  
vector<Float_t> Muon_pt;
vector<Float_t> Muon_e;
vector<Float_t> Muon_eta;
//vector<Float_t> Muon_mass;
*/
const double mumass = 0.105658;
const double emass = 0.0005109989;
const double taumass=0.177686;
const unsigned nReserve_Muon = 32, nReserve_Electron=32, nReserve_Tau=32;
Int_t TrigObj_id;
UInt_t  nMuon, nElectron, nTau;
Float_t Muon_mass[nReserve_Muon], Electron_mass[nReserve_Electron], Tau_mass[nReserve_Tau];
// declaring the muon arrays - won't be using vectors but I keep them still in comments if needed.
Int_t Muon_charge[nReserve_Muon];
Float_t Muon_dxy[nReserve_Muon], Muon_dxyErr[nReserve_Muon], Muon_dz[nReserve_Muon], Muon_pfRelIso03_all[nReserve_Muon], Muon_pfRelIso03_chg[nReserve_Muon], Muon_pfRelIso04_all[nReserve_Muon];
Float_t Muon_phi[nReserve_Muon], Muon_pt[nReserve_Muon], Muon_e[nReserve_Muon], Muon_eta[nReserve_Muon],Muon_sip3d[nReserve_Muon],Muon_ip3d[nReserve_Muon];
Int_t Muon_genPartIdx[nReserve_Muon];   // does not exist in CMS nanoAODv8, to be removed 
Bool_t Muon_isPFcand[nReserve_Muon], Muon_isGlobal[nReserve_Muon];
cout << "Muon arrays declared" << endl;

// declaring the electron arrays - no vectors were declared earlier so nothing in comments :(
  
Int_t Electron_charge[nReserve_Electron];
Float_t Electron_dxy[nReserve_Electron], Electron_dxyErr[nReserve_Electron], Electron_dz[nReserve_Electron], Electron_pfRelIso03_all[nReserve_Electron], Electron_pfRelIso03_chg[nReserve_Electron], Electron_pfRelIso04_all[nReserve_Electron];
Float_t Electron_phi[nReserve_Electron], Electron_pt[nReserve_Electron], Electron_e[nReserve_Electron], Electron_eta[nReserve_Electron],Electron_sip3d[nReserve_Electron],Electron_ip3d[nReserve_Electron];
Int_t Electron_genPartIdx[nReserve_Electron];   // does not exist in CMS nanoAODv8, to be removed 
UChar_t Electron_lostHits[nReserve_Electron];
Bool_t Electron_isPFcand[nReserve_Electron]; 
Bool_t Electron_isEE[nReserve_Electron],Electron_isEB[nReserve_Electron];
cout << "Electron arrays declared" << endl;

Int_t Tau_charge[nReserve_Tau];
Float_t Tau_dxy[nReserve_Tau], Tau_dxyErr[nReserve_Tau], Tau_dz[nReserve_Tau], Tau_pfRelIso03_all[nReserve_Tau], Tau_pfRelIso03_chg[nReserve_Tau], Tau_pfRelIso04_all[nReserve_Tau];
Float_t Tau_phi[nReserve_Tau], Tau_pt[nReserve_Tau], Tau_e[nReserve_Tau], Tau_eta[nReserve_Tau];
Int_t Tau_genPartIdx[nReserve_Tau];   // does not exist in CMS nanoAODv8, to be removed 
cout << "Tau arrays declared" << endl;


//--------------------------------- Muon reserve -----------------------------//
/*
  Muon_charge.reserve(nReserve_Muon);
  Muon_dxy.reserve(nReserve_Muon);
  Muon_dxyErr.reserve(nReserve_Muon);
  Muon_dz.reserve(nReserve_Muon);

  Muon_genPartIdx.reserve(nReserve_Muon); 
  //Muon_mass.reserve(nReserve_Muon);
  Muon_pfRelIso03_all.reserve(nReserve_Muon);
  Muon_pfRelIso03_chg.reserve(nReserve_Muon);

  Muon_phi.reserve(nReserve_Muon);  
  Muon_pt.reserve(nReserve_Muon);
  Muon_eta.reserve(nReserve_Muon);

*/
cout << "Memory reserved" << endl;
//--------------------------------- electron reserve -----------------------------//


//output branches

t_event->Branch("run", &run, "run/i");
t_event->Branch("event", &event, "event/l");
t_event->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
t_event->Branch("Generator_weight", &Generator_weight, "Generator_weight/F");
t_event->Branch("PV_z", &PV_z, "PV_z/F");
t_event->Branch("PV_npvs", &PV_npvs, "PV_npvs/I");
t_event->Branch("TrigObj_id",&TrigObj_id,"TrigObj_id/I");

//-------------------------- Create branch of Muons's tree -------------------//

  t_event->Branch("nMuon", &nMuon, "nMuon/i");
  t_event->Branch("Muon_charge", Muon_charge, "Muon_charge[nMuon]/I");
  t_event->Branch("Muon_mass", Muon_mass, "Muon_mass[nMuon]/F");
  t_event->Branch("Muon_dxy", Muon_dxy, "Muon_dxy[nMuon]/F");
  t_event->Branch("Muon_dxyErr", Muon_dxyErr, "Muon_dxyErr[nMuon]/F");
  t_event->Branch("Muon_dz", Muon_dz, "Muon_dz[nMuon]/F");
  
  t_event->Branch("Muon_eta", Muon_eta, "Muon_eta[nMuon]/F"); 
  t_event->Branch("Muon_genPartIdx", Muon_genPartIdx, "Muon_genPartIdx[nMuon]/B");
  t_event->Branch("Muon_pfRelIso03_all", Muon_pfRelIso03_all, "Muon_pfRelIso03_all[nMuon]/F");
  t_event->Branch("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, "Muon_pfRelIso03_chg[nMuon]/F");  
  t_event->Branch("Muon_pfRelIso04_all", Muon_pfRelIso04_all, "Muon_pfRelIso04_all[nMuon]/F");
  
  t_event->Branch("Muon_phi", Muon_phi, "Muon_phi[nMuon]/F");
  t_event->Branch("Muon_pt", Muon_pt, "Muon_pt[nMuon]/F");
 // t_event->Branch("Muon_e", Muon_e, "Muon_e[nMuon]/F");
  t_event->Branch("Muon_sip3d", Muon_sip3d, "Muon_sip3d[nMuon]/F");
  t_event->Branch("Muon_ip3d", Muon_ip3d, "Muon_ip3d[nMuon]/F");
  
  t_event->Branch("Muon_isPFcand", Muon_isPFcand,"Muon_isPFcand[nMuon]/O");
  t_event->Branch("Muon_isGlobal", Muon_isGlobal,"Muon_isGlobal[nMuon]/O");


//*****************************************************//

  cout << "Muon Branches declared" << endl;
// Muon branches have been declared 

//-------------------------- Create branch of Electrons's tree -------------------//

  t_event->Branch("nElectron", &nElectron, "nElectron/i");
  t_event->Branch("Electron_charge", Electron_charge, "Electron_charge[nElectron]/I");
  t_event->Branch("Electron_mass", Electron_mass, "Electron_mass[nMuon]/F");
  t_event->Branch("Electron_dxy", Electron_dxy, "Electron_dxy[nElectron]/F");
  t_event->Branch("Electron_dxyErr", Electron_dxyErr, "Electron_dxyErr[nElectron]/F");
  t_event->Branch("Electron_dz", Electron_dz, "Electron_dz[nElectron]/F");
  
  t_event->Branch("Electron_eta", Electron_eta, "Electron_eta[nElectron]/F"); 
  t_event->Branch("Electron_genPartIdx", Electron_genPartIdx, "Electron_genPartIdx[nElectron]/B");
  t_event->Branch("Electron_pfRelIso03_all", Electron_pfRelIso03_all, "Electron_pfRelIso03_all[nElectron]/F");
  t_event->Branch("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, "Electron_pfRelIso03_chg[nElectron]/F");  
  t_event->Branch("Electron_sip3d", Electron_sip3d, "Electron_sip3d[nElectron]/F");
  t_event->Branch("Electron_ip3d", Electron_ip3d, "Electron_ip3d[nElectron]/F");
  //t_event->Branch("Electron_pfRelIso04_all", Electron_pfRelIso04_all, "Electron_pfRelIso04_all[nElectron]/F");
  

  t_event->Branch("Electron_phi", Electron_phi, "Electron_phi[nElectron]/F");
  t_event->Branch("Electron_pt", Electron_pt, "Electron_pt[nElectron]/F");
  t_event->Branch("Electron_e", Electron_e, "Electron_e[nElectron]/F");
  t_event->Branch("Electron_lostHits", Electron_lostHits, "Electron_lostHits[nElectron]/b");
  t_event->Branch("Electron_isPFcand", Electron_isPFcand,"Electron_isPFcand[nElectron]/O");
  t_event->Branch("Electron_isEE", Electron_isEE, "Electron_isEE[nElectron]/O");
  t_event->Branch("Electron_isEB", Electron_isEB, "Electron_isEB[nElectron]/O");
  cout << "Electron Branches declared" << endl;

//*****************************************************//


//-------------------------- Create branch of Taus's tree -------------------//

  t_event->Branch("nTau", &nTau, "nTau/i");
  t_event->Branch("Tau_charge", Tau_charge, "Tau_charge[nTau]/I");
  t_event->Branch("Tau_dxy", Tau_dxy, "Tau_dxy[nTau]/F");
  t_event->Branch("Tau_dxyErr", Tau_dxyErr, "Tau_dxyErr[nTau]/F");
  t_event->Branch("Tau_dz", Tau_dz, "Tau_dz[nTau]/F");
  
  t_event->Branch("Tau_eta", Tau_eta, "Tau_eta[nTau]/F"); 
  t_event->Branch("Tau_genPartIdx", Tau_genPartIdx, "Tau_genPartIdx[nTau]/B");
  t_event->Branch("Tau_pfRelIso03_all", Tau_pfRelIso03_all, "Tau_pfRelIso03_all[nTau]/F");
  t_event->Branch("Tau_pfRelIso03_chg", Tau_pfRelIso03_chg, "Tau_pfRelIso03_chg[nTau]/F");  
  //t_event->Branch("Tau_pfRelIso04_all", Tau_pfRelIso04_all, "Tau_pfRelIso04_all[nTau]/F");
  

  t_event->Branch("Tau_phi", Tau_phi, "Tau_phi[nTau]/F");
  t_event->Branch("Tau_pt", Tau_pt, "Tau_pt[nTau]/F");
  t_event->Branch("Tau_e", Tau_e, "Tau_e[nTau]/F");
  


  cout << "Tau Branches declared" << endl;

//*****************************************************//

//data transfer

Int_t nentries = (Int_t)ntree->GetEntries();
cout << "Total number of entries: " << nentries << endl; //added whilst debugging, shows the neumber of entries in the input file //edit
//ntree->Print(); //added for debugging, now commented out //edit
//ntree->Show(0); //added for debugging, now commented out //edit

/*run.clear();
event.clear();
luminosityBlock.clear();
Generator_weight.clear();
PV_z.clear();
PV_npvs.clear(); */ 


for (Int_t i = 0; i<nentries; i++) 
{
  if (i%((int)nentries/10)==0) cout << "Iteration no. " <<i <<endl;
  //cout << "A" << endl;
  TStopwatch sw; //added for debugging //edit
  sw.Start(); //added for debugging //edit
  ntree->GetEntry(i);
  cout << "Trying GetEntry(" << i << ")..." << endl; //added for debugging //edit
  ntree->GetEntry(i);
  cout << "Loaded entry " << i << endl; //added for debugging //edit
  cout << "lep_n = " << lep_n << ", lep_pt->size() = " << lep_pt->size() << endl; //added for debugging //edit


  /*  Keeping the vector part of the code intact if needed later 
  Muon_phi.clear();  
  Muon_pt.clear(); 
  Muon_eta.clear();
  Muon_charge.clear(); 
  Muon_dxy.clear();
  Muon_dxyErr.clear();
  Muon_dz.clear();
  Muon_genPartIdx.clear();  
  //Muon_mass.clear(); 
  Muon_pfRelIso03_all.clear();
  Muon_pfRelIso03_chg.clear();
  */
  
  
  int muon_num=0,ele_num=0,tau_num=0;

  run=(runNumber);
  event=(eventNumber);
  luminosityBlock=(channelNumber);
  Generator_weight=(mcWeight);
  
  //had to comment the follwing bit out because it is also throwing an error //edit
  /*
  if (year==2012)
  {
  PV_npvs=(pvxp_n);
  PV_z=(vxp_z); 
  }  
  */
  //cout << "B" << endl;
  for (UInt_t bb = 0; bb < lep_n; bb++) 
  { 
    
    // muon filtered data 
    if (lep_type->at(bb) == pdgmuon)
    {   // muons
      muon_num++;
      Muon_pt[bb]=(lep_pt->at(bb)/1000.); // MeV to GeV conversion | ATLAS to CMS
      Muon_eta[bb]=(lep_eta->at(bb));
      Muon_phi[bb]=(lep_phi->at(bb));
      Muon_mass[bb]=mumass;
      
      Muon_charge[bb]=(lep_charge->at(bb));
    //  Muon_e[bb]=(lep_E->at(bb)/1000.);  // This can be calculated from pt and eta and mass so its not actually there in nanoAODv8, will be removed later. 
     
     //dummy variables not dependent on ATLAS data
      Muon_isPFcand[bb]=true;
      Muon_isGlobal[bb]=true;   // Quality cuts already applied by ATLAS, so set to true by default. //physlite not skimmed or cut initially so will have to change eventually, just not right now //edit
      Muon_sip3d[bb]=0.0;  
      /*************************************************************/
      Muon_ip3d[bb]=sqrt(lep_trackd0pvunbiased->at(bb)*lep_trackd0pvunbiased->at(bb)+lep_z0->at(bb)*lep_z0->at(bb))/10.0;

      Muon_pfRelIso03_all[bb]=(lep_etcone20->at(bb)/lep_pt->at(bb))*(9.0/4.0);
      Muon_pfRelIso03_chg[bb]=(lep_ptcone30->at(bb)/lep_pt->at(bb));
      Muon_pfRelIso04_all[bb]=Muon_pfRelIso03_all[bb]*(16.0/9.0);
      Muon_dxy[bb]=(lep_trackd0pvunbiased->at(bb))/10.0;
      Muon_dxyErr[bb]=(lep_tracksigd0pvunbiased->at(bb))/10.0;
      Muon_dz[bb]=(lep_z0->at(bb))/10.0;
    
     /* if (bb<lep_truthMatched->size())     // this variable doesn't exist in nanoAODv8 so it's kept commented out, can be restored | same for e and tau
      Muon_genPartIdx[bb]=(lep_truthMatched->at(bb));
      else 
      cout << "WARNING: size exceeded for Muon_genPartIdx" << endl; */  
    
    }  // end of muon triggered data
   
    if (lep_type->at(bb) == pdgelec)
    {   // Electrons
      ele_num++;
      Electron_pt[bb]=(lep_pt->at(bb)/1000.); // MeV to GeV conversion | ATLAS to CMS
      Electron_eta[bb]=(lep_eta->at(bb));
      Electron_phi[bb]=(lep_phi->at(bb));
      Electron_mass[bb]=emass;
      Electron_charge[bb]=(lep_charge->at(bb));
  //    Electron_e[bb]=(lep_E->at(bb)/1000.);
      
      if (abs(lep_eta->at(bb))<1.37) {
       Electron_isEE[bb]=false;
       Electron_isEB[bb]=true;
      }
      if (abs(lep_eta->at(bb))>1.52) {
       Electron_isEE[bb]=true;
       Electron_isEB[bb]=false;
      }
      
      //dummy variables not dependent on ATLAS data
      Electron_sip3d[bb]=0.0;
      Electron_isPFcand[bb]=true;

      /*************************************************************/

      Electron_pfRelIso03_all[bb]=(lep_etcone20->at(bb)/lep_pt->at(bb))*(9.0/4.0);  // scaling by radius  
      Electron_pfRelIso03_chg[bb]=(lep_ptcone30->at(bb)/lep_pt->at(bb));
      //Electron_pfRelIso04_all[bb]=Electron_pfRelIso03_all[bb]*(16.0/9.0);
      
      // The two branches Electron_isBB and Electron_isEE have not been included here, they can be dealt with in the individual analyses. 
      Electron_ip3d[bb]=sqrt(lep_trackd0pvunbiased->at(bb)*lep_trackd0pvunbiased->at(bb)+lep_z0->at(bb)*lep_z0->at(bb))/10.0;
      
      Electron_dxy[bb]=(lep_trackd0pvunbiased->at(bb))/10.0;
      Electron_dxyErr[bb]=(lep_tracksigd0pvunbiased->at(bb))/10.0;
      Electron_dz[bb]=(lep_z0->at(bb))/10.0;
      Electron_lostHits[bb]=0; // Set this to zero so that in the Higgs decay analysis, its always satisfied.
      /*if (bb<lep_truthMatched->size())
      Electron_genPartIdx[bb]=(lep_truthMatched->at(bb));
      else 
      cout << "WARNING: size exceeded for Electron_genPartIdx" << endl; */  
    
    }  // end of Electron triggered data

/*  if (lep_type->at(bb) == pdgtau)
    {   // Taus
      tau_num++;
      Tau_pt[bb]=(lep_pt->at(bb)/1000.); // MeV to GeV conversion | ATLAS to CMS
      Tau_eta[bb]=(lep_eta->at(bb));
      Tau_phi[bb]=(lep_phi->at(bb));
      Tau_mass[bb]=taumass;
      Tau_charge[bb]=(lep_charge->at(bb));
      Tau_e[bb]=(lep_E->at(bb)/1000.);

      Tau_pfRelIso03_all[bb]=(lep_etcone20->at(bb)/lep_pt->at(bb))*(9.0/4.0);
      Tau_pfRelIso03_chg[bb]=(lep_ptcone30->at(bb)/lep_pt->at(bb));
      //Tau_pfRelIso04_all[bb]=Tau_pfRelIso03_all[bb]*(16.0/9.0);
      
      Tau_dxy[bb]=(lep_trackd0pvunbiased->at(bb));
      Tau_dxyErr[bb]=(lep_tracksigd0pvunbiased->at(bb));
      Tau_dz[bb]=(lep_z0->at(bb));
    
    if (bb<lep_truthMatched->size())
      Tau_genPartIdx[bb]=(lep_truthMatched->at(bb));
      else 
      cout << "WARNING: size exceeded for Tau_genPartIdx" << endl;   
    
    }  // end of Tau triggered data*/

  // ends the sorter loop
  }    // ends the loop over leptons

  nMuon = muon_num;   
  nElectron = ele_num;
  nTau = tau_num;
  if (trigE) TrigObj_id = pdgelec;
  if (trigM) TrigObj_id = pdgmuon;
 
 // sorting the arrays now in order of Pt
  if (lep_n == 0) continue;  // skip events with no leptons to avoid sorting empty arrays, any event with no leptons would cause problems //edit
  for (UInt_t bb = 0; bb < lep_n-1; bb++) 
  {
   for (UInt_t bb1 = bb+1; bb1 < lep_n; bb1++) 
    {
      if (Muon_pt[bb]<Muon_pt[bb1])
      {
        Float_t temp1=Muon_pt[bb];
        Muon_pt[bb]=Muon_pt[bb1];
        Muon_pt[bb1]=temp1;

        temp1=Muon_eta[bb];
        Muon_eta[bb]=Muon_eta[bb1];
        Muon_eta[bb1]=temp1;

        temp1=Muon_phi[bb];
        Muon_phi[bb]=Muon_phi[bb1];
        Muon_phi[bb1]=temp1;

        temp1=Muon_dz[bb];
        Muon_dz[bb]=Muon_dz[bb1];
        Muon_dz[bb1]=temp1;
        
        temp1=Muon_dxy[bb];
        Muon_dxy[bb]=Muon_dxy[bb1];
        Muon_dxy[bb1]=temp1;

        temp1=Muon_dxyErr[bb];
        Muon_dxyErr[bb]=Muon_dxyErr[bb1];
        Muon_dxyErr[bb1]=temp1;

        temp1=Muon_pfRelIso03_all[bb];
        Muon_pfRelIso03_all[bb]=Muon_pfRelIso03_all[bb1];
        Muon_pfRelIso03_all[bb1]=temp1;

        temp1=Muon_pfRelIso03_chg[bb];
        Muon_pfRelIso03_chg[bb]=Muon_pfRelIso03_chg[bb1];
        Muon_pfRelIso03_chg[bb1]=temp1;

        temp1=Muon_pfRelIso04_all[bb];
        Muon_pfRelIso04_all[bb]=Muon_pfRelIso04_all[bb1];
        Muon_pfRelIso04_all[bb1]=temp1;
        
        temp1=Muon_ip3d[bb];
        Muon_ip3d[bb]=Muon_ip3d[bb1];
        Muon_ip3d[bb1]=temp1;

        Int_t temp2=Muon_charge[bb];
        Muon_charge[bb]=Muon_charge[bb1];
        Muon_charge[bb1]=temp2;

      }
    } 
  }
   for (UInt_t bb = 0; bb < lep_n-1; bb++) 
  {
   for (UInt_t bb1 = bb+1; bb1 < lep_n; bb1++) 
    {
      if (Electron_pt[bb]<Electron_pt[bb1])
      {
        Float_t temp1=Electron_pt[bb];
        Electron_pt[bb]=Electron_pt[bb1];
        Electron_pt[bb1]=temp1;

        temp1=Electron_eta[bb];
        Electron_eta[bb]=Electron_eta[bb1];
        Electron_eta[bb1]=temp1;

        temp1=Electron_phi[bb];
        Electron_phi[bb]=Electron_phi[bb1];
        Electron_phi[bb1]=temp1;

        temp1=Electron_dz[bb];
        Electron_dz[bb]=Electron_dz[bb1];
        Electron_dz[bb1]=temp1;
        
        temp1=Electron_dxy[bb];
        Electron_dxy[bb]=Electron_dxy[bb1];
        Electron_dxy[bb1]=temp1;

        temp1=Electron_dxyErr[bb];
        Electron_dxyErr[bb]=Electron_dxyErr[bb1];
        Electron_dxyErr[bb1]=temp1;

        temp1=Electron_pfRelIso03_all[bb];
        Electron_pfRelIso03_all[bb]=Electron_pfRelIso03_all[bb1];
        Electron_pfRelIso03_all[bb1]=temp1;

        temp1=Electron_pfRelIso03_chg[bb];
        Electron_pfRelIso03_chg[bb]=Electron_pfRelIso03_chg[bb1];
        Electron_pfRelIso03_chg[bb1]=temp1;

        temp1=Electron_ip3d[bb];
        Electron_ip3d[bb]=Electron_ip3d[bb1];
        Electron_ip3d[bb1]=temp1;

        Int_t temp2=Electron_charge[bb];
        Electron_charge[bb]=Electron_charge[bb1];
        Electron_charge[bb1]=temp2;

        Bool_t temp3=Electron_isEE[bb];
        Electron_isEE[bb]=Electron_isEE[bb1];
        Electron_isEE[bb1]=temp3;

        temp3=Electron_isEB[bb];
        Electron_isEB[bb]=Electron_isEB[bb1];
        Electron_isEB[bb1]=temp3;
      }
    } 
  }




  t_event->Fill();   // no branch for Tau triggered events 
  sw.Stop(); //added for debugging //edit
  cout <<">>> Event " << i << " took " << sw.RealTime() << " sec" << endl; //added for debugging //edit
  
 
} // ends the loop over all entries


  //*****************************************************//

  //******************Write in OutFile*******************//

  file->cd();
  t_event->Write();

  cout << "File written" << endl;

  //*****************************************************//


  //***********************Closing***********************//

  file->Close();

  //gROOT->ProcessLine(".q");

  //*****************************************************//

}

/*
Part below adapted to deal with my file - now to hope and pray
*/

void atlas2016_to_cms_nanoAODplus_revived_adapted() {

  atlastocmsfile("/afs/desy.de/user/l/lukerlou/OpenData/data/submitDir-2025-07-15-1059-ad24/data-ANALYSIS/output.root","/afs/desy.de/user/l/lukerlou/OpenData/data/submitDir-2025-07-15-1059-ad24/data-ANALYSIS/output_cms.root",2016);
  //atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2012/Data/DataMuons.root","cms_atlas_data_new.root");
  //atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2016/Data/data_A.2lep.root","/nfs/dust/cms/user/geiser/ATLAS_CMS_MAP/cms_atlas_data_2016A_plus.root",2016);
  //atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2016/Data/data_B.2lep.root","/nfs/dust/cms/user/geiser/ATLAS_CMS_MAP/cms_atlas_data_2016B_plus.root",2016);
  //atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2016/Data/data_C.2lep.root","/nfs/dust/cms/user/geiser/ATLAS_CMS_MAP/cms_atlas_data_2016C_plus.root",2016);
  //atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2016/Data/data_D.2lep.root","/nfs/dust/cms/user/geiser/ATLAS_CMS_MAP/cms_atlas_data_2016D_plus.root",2016);

//  atlastocmsfile("/nfs/dust/cms/user/geiser/ATLAS/2012/MC/mc_147771.Zmumu.root","cms_atlas_mc_new.root",2012); //all here commented out as was old data, not my file //edit

}
