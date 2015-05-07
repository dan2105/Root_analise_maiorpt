#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TMath.h" 
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"  
#include "TLatex.h"
#include "TLorentzVector.h"
#include <math.h>
#include <algorithm>
#include <list>
#include <stdbool.h>

void histoaawwemuv3(){


  TFile *_file0 = new TFile("dpe_wwhepmc_13tev.root");  //Input File in root format
  
  TTree *T = (TTree *)_file0->Get("T");   // Acessing tree

  int const NMUMAX=100;
  int n_mu=0;
  double mu_pt[NMUMAX];
  double mu_px[NMUMAX];
  double mu_py[NMUMAX];
  double mu_pz[NMUMAX];
  double mu_energy[NMUMAX];
  double mu_eta[NMUMAX];
  double mu_phi[NMUMAX];

  int const NVMUMAX=100;
  int n_vmu=0;
  double vmu_pt[NVMUMAX];
  double vmu_px[NVMUMAX];
  double vmu_py[NVMUMAX];
  double vmu_pz[NVMUMAX];
  double vmu_energy[NVMUMAX];
  double vmu_eta[NVMUMAX];
  double vmu_phi[NVMUMAX];
  
  int const NEMAX=100;
  int n_e=0;
  double e_pt[NEMAX];
  double e_px[NEMAX];
  double e_py[NEMAX];
  double e_pz[NEMAX];
  double e_eta[NEMAX];
  double e_phi[NEMAX];
  double e_energy[NEMAX];

  int const NVEMAX=100;
  int n_ve=0;
  double ve_pt[NVEMAX];
  double ve_px[NVEMAX];
  double ve_py[NVEMAX];
  double ve_pz[NVEMAX];
  double ve_eta[NVEMAX];
  double ve_phi[NVEMAX];
  double ve_energy[NVEMAX];

  int const NCHGMAX=700;
  int n_chg=0;
  double chg_pt[NCHGMAX];
  double chg_px[NCHGMAX];
  double chg_py[NCHGMAX];
  double chg_pz[NCHGMAX];
  double chg_eta[NCHGMAX];
  double chg_phi[NCHGMAX];
  double chg_energy[NCHGMAX];

  int const GAMAMAX=50;
  int n_a=0;
  double a_pz[GAMAMAX];
  double a_energy[GAMAMAX];

  int const PROTONMAX=50;
  int n_proton=0;
  double proton_pt[PROTONMAX];
  double proton_px[PROTONMAX];
  double proton_py[PROTONMAX];
  double proton_pz[PROTONMAX];
  double proton_energy[PROTONMAX];

  // Branch of Muons
  T->SetBranchAddress("n_mu",&n_mu);
  T->SetBranchAddress("mu_pt",&mu_pt);
  T->SetBranchAddress("mu_px",&mu_px);
  T->SetBranchAddress("mu_py",&mu_py);
  T->SetBranchAddress("mu_pz",&mu_pz);
  T->SetBranchAddress("mu_energy",&mu_energy);
  T->SetBranchAddress("mu_eta",&mu_eta);
  T->SetBranchAddress("mu_phi",&mu_phi);
  // Branch of Muons neutrins
  T->SetBranchAddress("n_vmu",&n_vmu);
  T->SetBranchAddress("vmu_pt",&vmu_pt);
  T->SetBranchAddress("vmu_px",&vmu_px);
  T->SetBranchAddress("vmu_py",&vmu_py);
  T->SetBranchAddress("vmu_pz",&vmu_pz);
  T->SetBranchAddress("vmu_energy",&vmu_energy);
  T->SetBranchAddress("vmu_eta",&vmu_eta);
  T->SetBranchAddress("vmu_phi",&vmu_phi);
  // Branch of Eletrons
  T->SetBranchAddress("n_e",&n_e);
  T->SetBranchAddress("e_pt",&e_pt);
  T->SetBranchAddress("e_px",&e_px);
  T->SetBranchAddress("e_py",&e_py);
  T->SetBranchAddress("e_pz",&e_pz);
  T->SetBranchAddress("e_eta",&e_eta);
  T->SetBranchAddress("e_phi",&e_phi);
  T->SetBranchAddress("e_energy",&e_energy);
  // Branch of Eletrons neutrins
  T->SetBranchAddress("n_ve",&n_ve);
  T->SetBranchAddress("ve_pt",&ve_pt);
  T->SetBranchAddress("ve_px",&ve_px);
  T->SetBranchAddress("ve_py",&ve_py);
  T->SetBranchAddress("ve_pz",&ve_pz);
  T->SetBranchAddress("ve_eta",&ve_eta);
  T->SetBranchAddress("ve_phi",&ve_phi);
  T->SetBranchAddress("ve_energy",&ve_energy);
  // Branch of charged particles
  T->SetBranchAddress("n_chg",&n_chg);
  T->SetBranchAddress("chg_pt",&chg_pt);
  T->SetBranchAddress("chg_px",&chg_px);
  T->SetBranchAddress("chg_py",&chg_py);
  T->SetBranchAddress("chg_pz",&chg_pz);
  T->SetBranchAddress("chg_eta",&chg_eta);
  T->SetBranchAddress("chg_phi",&chg_phi);
  T->SetBranchAddress("chg_energy",&chg_energy);
  // Branch of photons
  T->SetBranchAddress("n_a",&n_a);
  T->SetBranchAddress("a_pz",&a_pz);
  T->SetBranchAddress("a_energy",&a_energy);
  //Branch of Protons
  T->SetBranchAddress("n_proton",&n_proton);
  T->SetBranchAddress("proton_px",&proton_px);
  T->SetBranchAddress("proton_py",&proton_py);
  T->SetBranchAddress("proton_pz",&proton_pz);
  T->SetBranchAddress("proton_pt",&proton_pt);
  T->SetBranchAddress("proton_energy",&proton_energy);

  //===============================================================
  //Histograms setup
  //===============================================================

  //===============================================================
  //Muon's variables
  //===============================================================
  TH1F* h_nmu                = new TH1F("nmu","nmu",100,0.,20.);
  TH1F* h_mu_pt              = new TH1F("mu_pt","mu_pt",1000,0.,1250.);
  TH1F* h_mu_pt_proton       = new TH1F("mu_pt_proton","mu_pt_proton",1000,0.,1250.);
  //================================================================
  TH1F* h_mu_px          = new TH1F("mu_px","mu_px",1000,0.,1250.);
  TH1F* h_mu_py          = new TH1F("mu_py","mu_py",1000,0.,1250.);
  TH1F* h_mu_pz          = new TH1F("mu_pz","mu_pz",1000,0.,1250.);
  TH1F* h_mu_eta         = new TH1F("mu_eta","mu_eta",1000,-2.4,2.4);
  TH1F* h_mu_phi         = new TH1F("mu_phi","mu_phi",1000,-M_PI,M_PI);
  TH1F* h_mu_energy      = new TH1F("mu_energy","mu_energy",1000,0.,1000.);
  TH1F* h_mu_Tenergy     = new TH1F("mu_Tenergy","mu_Tenergy",1000,0.,1000.);
  //===============================================================
  //Eletron's variables
  //===============================================================
  TH1F* h_ne        = new TH1F("ne","ne",100,0.,20.);
  TH1F* h_e_pt      = new TH1F("e_pt","e_pt",1000,0.,1250.);
  TH1F* h_e_pt_proton = new TH1F("e_pt_proton","e_pt_proton",1000,0.,1250.);
  TH1F* h_e_px      = new TH1F("e_px","e_px",1000,0.,1250.);
  TH1F* h_e_py      = new TH1F("e_py","e_py",1000,0.,1250.);
  TH1F* h_e_pz      = new TH1F("e_pz","e_pz",1000,0.,1250.);
  TH1F* h_e_eta     = new TH1F("e_eta","e_eta",1000,-2.5,2.5);
  TH1F* h_e_phi     = new TH1F("e_phi","e_phi",1000,-M_PI,M_PI);
  TH1F* h_e_energy  = new TH1F("e_energy","e_energy",100,0.,1250.);
  TH1F* h_e_Tenergy = new TH1F("e_Tenergy","e_Tenergy",100,0.,1250.);
  //========================================================================
  TH1F* h_nvmu       = new TH1F("nvmu","nvmu",100,0.,20.);
  TH1F* h_vmu_pt     = new TH1F("vmu_pt","vmu_pt",1000,0.,1250.);
  TH1F* h_vmu_px     = new TH1F("vmu_px","vmu_px",1000,0.,1250.);
  TH1F* h_vmu_py     = new TH1F("vmu_py","vmu_py",1000,0.,1250.);
  TH1F* h_vmu_pz     = new TH1F("vmu_pz","vmu_pz",1000,0.,1250.);
  TH1F* h_vmu_eta    = new TH1F("vmu_eta","vmu_eta",1000,-2.5,2.5);
  TH1F* h_vmu_phi    = new TH1F("vmu_phi","vmu_phi",1000,-M_PI,M_PI);
  TH1F* h_vmu_energy = new TH1F("vmu_energy","vmu_energy",100,0.,1250.);
  //===============================================================
  //Eletron's neutrins variables
  //===============================================================
  TH1F* h_nve       = new TH1F("nve","nve",1000,0.,20.);
  TH1F* h_ve_pt     = new TH1F("ve_pt","ve_pt",1000,0.,1250.);
  TH1F* h_ve_px     = new TH1F("ve_px","ve_px",1000,0.,1250.);
  TH1F* h_ve_py     = new TH1F("ve_py","ve_py",1000,0.,1250.);
  TH1F* h_ve_pz     = new TH1F("ve_pz","ve_pz",1000,0.,1250.);
  TH1F* h_ve_eta    = new TH1F("ve_eta","ve_eta",1000,-2.5,2.5);
  TH1F* h_ve_phi    = new TH1F("ve_phi","ve_phi",1000,-M_PI,M_PI);
  TH1F* h_ve_energy = new TH1F("ve_energy","ve_energy",1000,0.,1250.);
  //===============================================================
  //Charged particle's variables
  //===============================================================
  TH1F* h_nchg        = new TH1F("nchg","nchg",100,0.,100.);
  TH1F* h_chg_px      = new TH1F("chg_px","chg_px",1000,0.,1250.);
  TH1F* h_chg_py      = new TH1F("chg_py","chg_py",1000,0.,1250.);
  TH1F* h_chg_pz      = new TH1F("chg_pz","chg_pz",1000,0.,1250.);
  TH1F* h_chg_pt      = new TH1F("chg_pt","chg_pt",1000,0.,1250.);
  TH1F* h_chg_eta     = new TH1F("chg_eta","chg_eta",1000,-5.0,5.0);
  TH1F* h_chg_phi     = new TH1F("chg_phi","chg_phi",1000,-M_PI,M_PI);
  TH1F* h_chg_energy  = new TH1F("chg_energy","chg_energy",1000,0,1250);
  TH1F* h_chg_Tenergy = new TH1F("chg_Tenergy","chg_Tenergy",1000,0,1250);
  //===============================================================
  // Selected charged particle's variables
  //===============================================================
  TH1F* h_nchg_sel   = new TH1F("nchg_sel","nchg_sel",100,0.,100.);
  TH1F* h_chg_ptdiff = new TH1F("chg_ptdiff","chg_ptdiff",1000,0.,1500.);
  TH1F* h_chg_ptsel  = new TH1F("chg_ptsel","chg_ptsel",2000,0.,1000.);
  //======================================================================
  TH1F* h_somaemu_pt = new TH1F("somaemu_pt","somaemu_pt", 100,0.,1000.);
  TH1F* h_somaemu_m  = new TH1F("somaemu_m","somaemu_m", 100,0.,1000.);
  //===========================================================
  TH1F* h_dphiemu = new TH1F("dphiemu","dphiemu",1000,0,M_PI);
  TH1F* h_emu_aco = new TH1F("emu_aco","emu_aco",1000,0.,1);
  //======================================================================
  //Final state protons
  //===================================================================
  TH1F* h_nproton     = new TH1F("nproton","nproton",10,0.,10.);
  TH1F* h_nproton_sel     = new TH1F("nproton_sel","nproton_sel",10,0.,10.);
  TH1F* h_proton_pt   = new TH1F("proton_pt","proton_pt",100,0.,5.);
  TH1F* h_proton1_pt  = new TH1F("proton1_pt","proton1_pt",100,0.,10.);
  TH1F* h_proton2_pt  = new TH1F("proton2_pt","proton2_pt",100,0.,10.);
  TH1F* h_proton1_eta  = new TH1F("proton1_eta","proton1_eta",100,-20.,20.);
  TH1F* h_proton2_eta  = new TH1F("proton2_eta","proton2_eta",100,-20.,20.);
  TH1F* h_proton1_phi  = new TH1F("proton1_phi","proton1_phi",100,-3.14,3.14);
  TH1F* h_proton2_phi  = new TH1F("proton2_phi","proton2_phi",100,-3.14,3.14);
  TH1F* h_dphiproton = new TH1F("dphi_proton","dphi_proton",1000,0,M_PI);
  TH1F* h_protons_aco = new TH1F("proton_aco","proton_aco",1000,0.,1);
  TH1F* h_proton1_pt2 = new TH1F("proton1_pt2","proton1_pt2",100,0.,10.);
  TH1F* h_proton2_pt2 = new TH1F("proton2_pt2","proton2_pt2",100,0.,10.);
  TH1F* h_proton_px   = new TH1F("proton_px","proton_px",100,0.,5.);
  TH1F* h_proton_py   = new TH1F("proton_py","proton_py",100,0.,5.);
  TH1F* h_proton_pz   = new TH1F("proton_pz","proton_pz",100,0.,5.);
  TH1F* h_proton_energy = new TH1F("proton_energy","proton_energy",100,-4000.,4000.);
  TH1F* h_proton_pt2  = new TH1F("proton_pt2","proton_pt2",100,0.,5.);
  TH1F* h_proton1_xi  = new TH1F("proton1_xi","proton1_xi",1000,0.,1.);
  TH1F* h_proton2_xi  = new TH1F("proton2_xi","proton2_xi",1000,0.,1.);
  TH1F* h_ecm         = new TH1F("e_cm","e_cm",300,0.,3000.);
  TH1F* h_ecm_160_500 = new TH1F("e_cm_160_500","e_cm_160_500",300,0.,3000.);
  TH1F* h_ecm_500_1000 = new TH1F("e_cm_500_1000","e_cm_500_1000",300,0.,3000.);
  TH1F* h_ecm_1000 = new TH1F("e_cm_1000","e_cm_1000",300,0.,3000.);
  //======================================================================
  TH1F* h_W1mass   = new TH1F("W1mass","W1mass",100,0,100);
  TH1F* h_W2mass   = new TH1F("W2mass","W2mass",100,0,100);
  TH1F* h_WWmass   = new TH1F("WWmass","WWmass",500,0,1500);
  //============================================================================
  TH1F* h_vlpt         = new TH1F("vlpt","vlpt",1000,0,1500);
  TH1F* h_vlEt         = new TH1F("vlEt","vlEt",1000,0,1500);
  //======================================================================
  //tree loop 
  //======================================================================
  Int_t nentries = (Int_t)T->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    T->GetEntry(i);
    //======================================================================
    TLorentzVector chgp4;
    for(int chg = 0; chg < n_chg; chg++){
  
      double chg1_pt = chg_pt[chg];
      double chg1_eta =chg_eta[chg];
      double chg1_phi = chg_phi[chg];
      double chg1_energy = chg_energy[chg];
      double chg1_px = chg_px[chg];
      double chg1_py = chg_py[chg];
      double chg1_pz = chg_pz[chg];
    
      chgp4.SetPxPyPzE(chg_px[chg],chg_py[chg],chg_pz[chg],chg_energy[chg]);
    }  

    h_nchg->Fill(n_chg);          
    h_chg_px->Fill(chgp4.Px());
    h_chg_py->Fill(chgp4.Py());
    h_chg_pz->Fill(chgp4.Pz());
    // h_chg_eta->Fill(chgp4.Eta());
    h_chg_phi->Fill(chgp4.Phi());    
    h_chg_energy->Fill(chgp4.E());
    h_chg_Tenergy->Fill(chgp4.Et());
    h_chg_pt->Fill(chgp4.Pt());
    

    //======================================================================        
    //proton (Declare 4vector of each particle and make it acessible to the whole macro)
    //======================================================================
    TLorentzVector proton1p4,proton2p4, protonp4;
    for(int iproton = 0; iproton < n_proton; iproton++) {

      double proton1_px = proton_px[iproton];
      double proton1_py = proton_py[iproton];
      double proton1_pz = proton_pz[iproton];
      double proton1_pt = proton_pt[iproton];
      double proton1_energy = proton_energy[iproton];
      protonp4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);

      if(protonp4.Pz()>0)
	{
	  proton1p4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);
	}// 4 vector-Pomeron  Pz>0
      else{
	proton2p4.SetPxPyPzE(proton1_px,proton1_py,proton1_pz,proton1_energy);
      }// 4 vector-Reggeon Pz<0
    }

    h_proton_pt->Fill(protonp4.Pt());
    h_proton_pt2->Fill(protonp4.Pt()*protonp4.Pt());
    h_proton_px->Fill(protonp4.Px());
    h_proton_py->Fill(protonp4.Py());
    h_proton_pz->Fill(protonp4.Pz());
    h_proton_energy->Fill(protonp4.E());    
    h_nproton->Fill(n_proton);       
    
    double proton1_xi = 1-fabs(proton1p4.Pz())/6500;
    if(proton1_xi > 0.0015 && proton1_xi < 0.15){
      h_proton1_xi->Fill(proton1_xi);
      h_proton1_pt->Fill(proton1p4.Pt());
      h_proton1_pt2->Fill(proton1p4.Pt()*proton1p4.Pt());
      h_proton1_eta->Fill(proton1p4.Eta());
      h_proton1_phi->Fill(proton1p4.Phi());
    }
    double proton2_xi= 1-fabs(proton2p4.Pz())/6500;
    if(proton2_xi > 0.0015 && proton2_xi < 0.15){
      h_proton2_xi->Fill(proton2_xi);  
      h_proton2_pt->Fill(proton2p4.Pt());
      h_proton2_pt2->Fill(proton2p4.Pt()*proton2p4.Pt());
      h_proton2_eta->Fill(proton2p4.Eta()); 
      h_proton2_phi->Fill(proton2p4.Phi());
    }
    double ROOTS=13000;
    if(proton1_xi > 0.0015 && proton1_xi < 0.15 && proton2_xi > 0.0015 && proton2_xi < 0.15 ){
      double Ecm= sqrt((1-abs(proton1p4.Pz())/6500)*(1-abs(proton2p4.Pz())/6500)) *ROOTS ;
      h_ecm->Fill(Ecm);
      double dphi_proton= fabs(proton1p4.Phi()-proton2p4.Phi());
      if(dphi_proton>3.1415927) dphi_proton= 6.2831853-dphi_proton;
      h_dphiproton->Fill(dphi_proton);
      h_protons_aco-> Fill(1- dphi_proton/M_PI);  	
      h_nproton_sel->Fill(n_proton);       


    if (Ecm >160 && Ecm < 500){
      h_ecm_160_500->Fill(Ecm);
    }
    if (Ecm >500 && Ecm < 1000){
      h_ecm_500_1000->Fill(Ecm);
    }
    if (Ecm >1000 ){
      h_ecm_1000->Fill(Ecm);
    }

    }//proton loop

    //==========================================================================
    //Proton analysis in progress
    //===========================================================================

    //======================================================================
    //Eletron neutrins
    //======================================================================
    TLorentzVector ve1p4;
    for(int ive = 0; ive < n_ve; ive++) {
	    
      double ve1pt = ve_pt[ive];
      double ve1px = ve_px[ive];
      double ve1py = ve_py[ive];
      double ve1pz = ve_pz[ive];
      double ve1eta =ve_eta[ive];
      double ve1phi = ve_phi[ive];
      double ve1energy = ve_energy[ive];
      
      ve1p4.SetPxPyPzE(ve1px,ve1py,ve1pz,ve1energy);	
    }//end neutrin loop

    if(n_ve>0){
      h_ve_pt->Fill(ve1p4.Pt());
      h_ve_px->Fill(ve1p4.Px());
      h_ve_py->Fill(ve1p4.Py());
      h_ve_pz->Fill(ve1p4.Pz());
      h_ve_energy->Fill(ve1p4.E());    
      // h_ve_eta->Fill(ve1p4.Eta());    
      h_ve_phi->Fill(ve1p4.Phi());
      h_nve->Fill(n_ve);         	
    }

    //======================================================================
    //Eletron neutrins
    //======================================================================
    TLorentzVector vmu1p4;
    for(int ivmu = 0; ivmu < n_vmu; ivmu++) {
	    
      double vmu1pt = vmu_pt[ivmu];
      double vmu1px = vmu_px[ivmu];
      double vmu1py = vmu_py[ivmu];
      double vmu1pz = vmu_pz[ivmu];
      double vmu1eta =vmu_eta[ivmu];
      double vmu1phi = vmu_phi[ivmu];
      double vmu1energy = vmu_energy[ivmu];
	
      vmu1p4.SetPxPyPzE(vmu1px,vmu1py,vmu1pz,vmu1energy);
    }
    if(n_vmu>0){	
      h_vmu_pt->Fill(vmu1p4.Pt());
      h_vmu_px->Fill(vmu1p4.Px());
      h_vmu_py->Fill(vmu1p4.Py());
      h_vmu_pz->Fill(vmu1p4.Pz());
      h_vmu_energy->Fill(vmu1p4.E());    
      //h_vmu_eta->Fill(vmu1p4.Eta());    
      h_vmu_phi->Fill(vmu1p4.Phi());
      h_nvmu->Fill(n_vmu);         	
    }

    //======================================================================
    //Eletron
    //======================================================================
    TLorentzVector e1p4;
    int ie_maxpt = -1;     
    double e_pt_max=0.;
    for(int ie = 0; ie < n_e; ie++) {
      if( e_pt[ie] > e_pt_max ){
	ie_maxpt = ie;
	e_pt_max = e_pt[ie];
      }
    }
    if(ie_maxpt >0 || n_e > 0){
      double e1pt = e_pt[ie_maxpt];
      double e1px = e_px[ie_maxpt];
      double e1py = e_py[ie_maxpt];
      double e1pz = e_pz[ie_maxpt];
      double e1eta =e_eta[ie_maxpt];
      double e1phi = e_phi[ie_maxpt];
      double e1energy = e_energy[ie_maxpt];
      
      e1p4.SetPxPyPzE(e1px,e1py,e1pz,e1energy);
      
      if(e1p4.Pt()>10 && fabs(e1p4.Eta())<2.5){
	h_e_pt->Fill(e1p4.Pt());
	h_e_px->Fill(e1p4.Px());
	h_e_py->Fill(e1p4.Py());
	h_e_pz->Fill(e1p4.Pz());
	h_e_energy->Fill(e1p4.E());    
	h_e_Tenergy->Fill(e1p4.Et());    
	h_e_eta->Fill(e1p4.Eta());    
	h_e_phi->Fill(e1p4.Phi());
	h_ne->Fill(n_e);         
      }

      if(e1p4.Pt()>10 && fabs(e1p4.Eta())<2.5  && proton1_xi > 0.0015 && proton1_xi < 0.15 && proton2_xi > 0.0015 && proton2_xi < 0.15){
	h_e_pt_proton->Fill(e1p4.Pt());
      }

    }
    //======================================================================
    //muon transverse momentum//Selection of the highest muon pt
    //======================================================================
    TLorentzVector mu1p4;
    int imu_maxpt = -1;     
    double mu_pt_max=0.;
    for(int imu = 0; imu < n_mu; imu++) {
      if( mu_pt[imu] > mu_pt_max ){
	imu_maxpt = imu ;
	mu_pt_max = mu_pt[imu];
      }
    }
    //colocar o mesmo no eletron
    if(imu_maxpt >0 || n_mu > 0){
      double mu1pt = mu_pt[imu_maxpt];
      double mu1px = mu_px[imu_maxpt];
      double mu1py = mu_py[imu_maxpt];
      double mu1pz = mu_pz[imu_maxpt];
      double mu1eta = mu_eta[imu_maxpt];
      double mu1phi = mu_phi[imu_maxpt];
      double mu1energy = mu_energy[imu_maxpt];   
    
      mu1p4.SetPxPyPzE(mu1px,mu1py,mu1pz,mu1energy);
      

      //colocar a seleção nos léptons centrais 
      if(mu1p4.Pt()>10 && fabs(mu1p4.Eta())<2.5){        
	h_mu_pt->Fill(mu1p4.Pt());
	h_mu_px->Fill(mu1p4.Px());
	h_mu_py->Fill(mu1p4.Py());
	h_mu_pz->Fill(mu1p4.Pz());
	h_mu_eta->Fill(mu1p4.Eta());    
	h_mu_phi->Fill(mu1p4.Phi());    
	h_mu_energy->Fill(mu1p4.E());    
	h_mu_Tenergy->Fill(mu1p4.Et());    
	h_nmu->Fill(n_mu);         
      }

      if(mu1p4.Pt()>10 && fabs(mu1p4.Eta())<2.5  && proton1_xi > 0.0015 && proton1_xi < 0.15 && proton2_xi > 0.0015 && proton2_xi < 0.15){
	h_mu_pt_proton->Fill(mu1p4.Pt());
      }
    }

    //=====================================================
    //W variables
    //========================================================
    // A partir daqui apenas se eletron e muon foram encontrados, only if electrons and muons were found!
    
    double dphi= fabs(e1p4.Phi()-mu1p4.Phi());
    if(dphi>3.1415927) dphi= 6.2831853-dphi;
    h_dphiemu->Fill(dphi);
    h_emu_aco-> Fill(1- dphi/M_PI);  	
    
    
    if(ie_maxpt == -1 || imu_maxpt == -1 ) continue;       
    if( fabs(mu1p4.Eta()) > 2.4 || fabs(e1p4.Eta()) > 2.5) continue;
    if( mu1p4.Pt() < 10 || e1p4.Pt() < 10) continue;
    
    double Wref= 80.43;
    double W1mass = (mu1p4+vmu1p4).M();
    h_W1mass->Fill(W1mass);
    double W2mass = (e1p4+ve1p4).M();
    h_W2mass->Fill(W2mass);
    double WWmass =(mu1p4+vmu1p4+e1p4+ve1p4).M();
    h_WWmass->Fill(WWmass); //proton missing mass
    
    //==============================================================
    // Varíaveis de 4-momento com um lépton e um neutrino
    //==============================================================
    
    //  if(fabs(Wref-W1mass)<1 || fabs(Wref-W2mass)<1){ //leptons coming from Ws
    double WWmassEt =(mu1p4+vmu1p4+e1p4+ve1p4).Et();
    double emusumEt =(mu1p4+e1p4).Et();
    double emumassT =(mu1p4+e1p4).Mt();
    double missingpt = (vmu1p4+ve1p4).Pt();
    h_vlpt->Fill(missingpt);
    double missingEt = (vmu1p4+ve1p4).Et();
    h_vlEt->Fill(missingEt);
    double Minv = (mu1p4+e1p4).M();
    h_somaemu_m->Fill(Minv);          
    double emupt = (e1p4+mu1p4).Pt();
    h_somaemu_pt->Fill(emupt);	
  

    //=============================================================	  	
    //Track selection of the charged particles
    //=============================================================
    double dR_e = sqrt(pow(chgp4.Phi()-e1p4.Phi(),2)+ pow(chgp4.Eta()-e1p4.Eta(),2));
    double dR_mu = sqrt(  pow(chgp4.Phi()-mu1p4.Phi(),2)+ pow(chgp4.Eta()-mu1p4.Eta(),2));
    double dpt_e = fabs(chgp4.Pt()-e1p4.Pt());
    double dpt_mu = fabs(chgp4.Pt()-mu1p4.Pt());
    if(dR_e > 0.0001 && dR_mu >0.0001 && dpt_e>0.1 && dpt_mu>0.1)
      {
	h_chg_ptsel->Fill(chgp4.Pt());
	h_nchg_sel->Fill(n_chg);        
      }
  
  } //Fim do loop na tree
  
    
  //===================================================================
  //Normalization
  //====================================================================
   
  Double_t  L_int = 1;
  Double_t sigma_wwinc=2.9;  
  Double_t N_ww = 50000;    
  Double_t scale1 = ((sigma_wwinc * L_int)/(N_ww));  
  
  
  h_nmu->Scale(scale1);
  h_mu_px->Scale(scale1);
  h_mu_py->Scale(scale1);
  h_mu_pz->Scale(scale1);
  h_mu_pt->Scale(scale1);
  h_mu_pt_proton->Scale(scale1);
  h_mu_eta->Scale(scale1);
  h_mu_phi->Scale(scale1);
  h_mu_energy->Scale(scale1);
  h_mu_Tenergy->Scale(scale1);
  //===================================
  h_somaemu_pt->Scale(scale1);
  //================================================
  h_vlpt->Scale(scale1);
  h_vlEt->Scale(scale1);
  //=================================================    
  //====================================
  h_ne->Scale(scale1);
  h_e_px->Scale(scale1);
  h_e_py->Scale(scale1);
  h_e_pz->Scale(scale1);
  h_e_pt->Scale(scale1);
  h_e_pt_proton->Scale(scale1);
  h_e_eta->Scale(scale1);
  h_e_phi->Scale(scale1);
  h_e_energy->Scale(scale1);
  //====================================     
  h_nvmu->Scale(scale1);
  h_vmu_px->Scale(scale1);
  h_vmu_py->Scale(scale1);
  h_vmu_pz->Scale(scale1);
  h_vmu_pt->Scale(scale1);
  h_vmu_eta->Scale(scale1);
  h_vmu_phi->Scale(scale1);
  h_vmu_energy->Scale(scale1);
  //====================================
  h_nve->Scale(scale1);
  h_ve_px->Scale(scale1);
  h_ve_py->Scale(scale1);
  h_ve_pz->Scale(scale1);
  h_ve_pt->Scale(scale1);
  h_ve_eta->Scale(scale1);
  h_ve_phi->Scale(scale1);
  h_ve_energy->Scale(scale1);
  //==========================================================
  h_nchg->Scale(scale1);
  h_chg_px->Scale(scale1);
  h_chg_py->Scale(scale1);
  h_chg_pz->Scale(scale1);
  h_chg_pt->Scale(scale1);
  h_chg_eta->Scale(scale1);
  h_chg_energy->Scale(scale1);
  h_chg_Tenergy->Scale(scale1);
  h_chg_phi->Scale(scale1);
  h_nchg_sel->Scale(scale1);
  h_chg_ptsel->Scale(scale1);
  h_chg_ptdiff->Scale(scale1);
  //===================================================================
  h_emu_aco->Scale(scale1);
  h_somaemu_m->Scale(scale1);
  h_somaemu_pt->Scale(scale1);
  h_dphiemu->Scale(scale1);
  //===================================================================
  h_nproton->Scale(scale1);
  h_nproton_sel->Scale(scale1);
  h_proton_px->Scale(scale1);
  h_proton_py->Scale(scale1);
  h_proton_pz->Scale(scale1);
  h_proton_pt->Scale(scale1);
  h_proton1_pt->Scale(scale1);
  h_proton2_pt->Scale(scale1);
  h_proton1_eta->Scale(scale1);
  h_proton2_eta->Scale(scale1);
  h_proton1_phi->Scale(scale1);
  h_proton2_phi->Scale(scale1);
  h_dphiproton->Scale(scale1);
  h_protons_aco->Scale(scale1);  	
  h_proton1_pt2->Scale(scale1);
  h_proton2_pt2->Scale(scale1);
  h_proton_pt2->Scale(scale1);
  h_proton1_xi->Scale(scale1);
  h_proton2_xi->Scale(scale1);
  h_ecm->Scale(scale1);
  h_ecm_160_500->Scale(scale1);
  h_ecm_500_1000->Scale(scale1);
  h_ecm_1000->Scale(scale1);
  h_proton_energy->Scale(scale1);
  //===================================================================
  h_W1mass->Scale(scale1);
  h_W2mass->Scale(scale1);
  h_WWmass->Scale(scale1);
  //===================================================================

  //========================================================================
  //Output File
  TFile* output = new TFile("dpe_wwhepmclumi1_13tev.root","RECREATE");
  //========================================================================
  h_nmu->Write();     
  h_mu_pt->Write();
  h_mu_pt_proton->Write();
  h_mu_px->Write();
  h_mu_py->Write();
  h_mu_pz->Write();
  h_mu_energy->Write();
  h_mu_Tenergy->Write();
  h_mu_eta->Write();
  h_mu_phi->Write();
  //=================================================================
  h_ne->Write();     
  h_e_pt->Write();
  h_e_pt_proton->Write();
  h_e_px->Write();
  h_e_py->Write();
  h_e_pz->Write();
  h_e_energy->Write();
  h_e_Tenergy->Write();
  h_e_eta->Write();
  h_e_phi->Write();
  //===================================================================
  h_nvmu->Write();     
  h_vmu_pt->Write();
  h_vmu_px->Write();
  h_vmu_py->Write();
  h_vmu_pz->Write();
  h_vmu_energy->Write();
  h_vmu_eta->Write();
  h_vmu_phi->Write();
  //=================================================================
  h_nve->Write();     
  h_ve_pt->Write();
  h_ve_px->Write();
  h_ve_py->Write();
  h_ve_pz->Write();
  h_ve_energy->Write();
  h_ve_eta->Write();
  h_ve_phi->Write();
  //========================================================================
  h_nchg->Write();     
  h_nchg_sel->Write();     
  h_chg_pt->Write();
  h_chg_ptsel->Write();
  h_chg_px->Write();
  h_chg_py->Write();
  h_chg_pz->Write();
  h_chg_eta->Write();
  h_chg_energy->Write();
  h_chg_Tenergy->Write();
  h_chg_phi->Write();
  //========================================================================
  h_somaemu_pt->Write();
  h_somaemu_m->Write();
  //=======================================================================
  h_vlpt->Write();
  h_vlEt->Write();
  //========================================================================
  h_emu_aco->Write();
  h_dphiemu->Write();
  //========================================================================
  h_nproton->Write();
  h_nproton_sel->Write();
  h_proton_px->Write();
  h_proton_py->Write();
  h_proton_pz->Write();
  h_proton_pt->Write();
  h_proton1_pt->Write();
  h_proton2_pt->Write();
  h_proton1_eta->Write();
  h_proton2_eta->Write();
  h_proton1_phi->Write();
  h_proton2_phi->Write();
  h_dphiproton->Write();
  h_protons_aco->Write();  	
  h_proton1_pt2->Write();
  h_proton2_pt2->Write();
  h_proton_pt2->Write();
  h_proton1_xi->Write();
  h_proton2_xi->Write();
  h_ecm->Write();
  h_ecm_160_500->Write();
  h_ecm_500_1000->Write();
  h_ecm_1000->Write();
  h_proton_energy->Write();
  //=======================================================================
  h_W1mass->Write();
  h_W2mass->Write();
  h_WWmass->Write();

  output->Close();
} 

