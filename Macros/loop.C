// Includes //
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
//
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TEntryList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector2.h"

// Namespaces //
using namespace std;

// Typedefs //

//   lumi  events
typedef map< UInt_t , vector<ULong64_t> > MapLumiEvent;


// Declare functions //

Int_t processChain(TChain* fChain, Bool_t doCutOnMetFlags, Bool_t doPrintBeforeCuts);


// MAIN //
Int_t loop(Bool_t doCutOnMetFlags=kTRUE, Bool_t doPrintBeforeCuts=kTRUE)
{

  TChain* chain = new TChain("tree");
  chain->Add("/eos/cms/store/group/phys_exotica/xtracks/6Mar2019-Hadded/QCD_HT500to700/treeProducerXtracks/tree.root");

  processChain(chain, doCutOnMetFlags, doPrintBeforeCuts);

  return 0;
}

Int_t processChain(TChain* fChain, Bool_t doCutOnMetFlags, Bool_t doPrintBeforeCuts)
{

  // Log file
  ofstream outlog("log.txt", ios::out); //fixme
  ofstream outmap("map.txt", ios::out); //fixme

  // Map selected events
  MapLumiEvent mapEvents; // < lumi , <events> >
  MapLumiEvent::iterator itMapEvents;
  vector<ULong64_t> vEvents;

  // Branches //
  //
  // Event
  UInt_t          run;
  UInt_t          lumi;
  ULong64_t       evt;
  Int_t           isData;
  //
  // MC
  Float_t         xsec;
  Float_t         puWeight;
  Float_t         nTrueInt;
  Float_t         genWeight;
  Int_t           nVert;
  //
  // Trigger
  Int_t           HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Int_t           HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
  Int_t           HLT_MET;
  Int_t           HLT_BIT_HLT_IsoMu24_v;
  Int_t           HLT_BIT_HLT_IsoMu24_eta2p1_v;
  Int_t           HLT_BIT_HLT_IsoMu27_v;
  Int_t           HLT_SingleMu;
  Int_t           HLT_BIT_HLT_Ele32_WPTight_Gsf_v;
  Int_t           HLT_BIT_HLT_Ele35_WPTight_Gsf_v;
  Int_t           HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v;
  Int_t           HLT_SingleEl;
  //
  // MET Noise Filters
  Int_t           Flag_goodVertices;
  Int_t           Flag_BadPFMuonFilter;
  Int_t           Flag_HBHENoiseIsoFilter;
  Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
  Int_t           Flag_eeBadScFilter;
  Int_t           Flag_BadChargedCandidateFilter;
  Int_t           Flag_ecalBadCalibFilter;
  Int_t           Flag_HBHENoiseFilter;
  Int_t           Flag_globalTightHalo2016Filter;
  //
  // MET
  Float_t         met_sumEt;
  Float_t         met_rawPt;
  Float_t         met_rawPhi;
  Float_t         met_rawSumEt;
  Float_t         met_genPt;
  Float_t         met_genPhi;
  Float_t         met_genEta;
  Float_t         met_pt;
  Float_t         met_eta;
  Float_t         met_phi;
  Float_t         met_mass;
  Float_t         metNoMu_pt;
  Float_t         metNoMu_eta;
  Float_t         metNoMu_phi;
  Float_t         metNoMu_mass;
  //
  // Jets
  const UInt_t    NJET=15;
  Int_t           nJet;
  Int_t           nJet30;
  Int_t           nJet30a;
  Float_t         Jet_pt[NJET];
  Float_t         Jet_eta[NJET];
  Float_t         Jet_phi[NJET];
  Float_t         Jet_chHEF[NJET];
  Float_t         Jet_neHEF[NJET];
  //
  // Jets forward
  const UInt_t    NJETF=6;
  Int_t           nJetFwd;
  Float_t         JetFwd_pt[NJETF];
  Float_t         JetFwd_eta[NJETF];
  Float_t         JetFwd_phi[NJETF];
  Float_t         JetFwd_chHEF[NJETF];
  Float_t         JetFwd_neHEF[NJETF];
  //
  // Leptons
  Int_t           nTauGood;
  Int_t           nLepGood;
  //
  // Tracks
  const UInt_t    NTRACK=4;
  Int_t           nIsoTrack;
  //
  Float_t         IsoTrack_pt[NTRACK];   //[nIsoTrack]             
  Float_t         IsoTrack_eta[NTRACK];   //[nIsoTrack]            
  Float_t         IsoTrack_phi[NTRACK];   //[nIsoTrack]            
  Float_t         IsoTrack_mass[NTRACK];   //[nIsoTrack]           
  //
  Int_t           IsoTrack_charge[NTRACK];   //[nIsoTrack]         
  Float_t         IsoTrack_dxy[NTRACK];   //[nIsoTrack]            
  Float_t         IsoTrack_dz[NTRACK];   //[nIsoTrack]             
  Float_t         IsoTrack_edxy[NTRACK];   //[nIsoTrack]           
  Float_t         IsoTrack_edz[NTRACK];   //[nIsoTrack]            
  //
  Float_t         IsoTrack_miniRelIso[NTRACK];
  Float_t         IsoTrack_relIso03[NTRACK];
  //
  Int_t           IsoTrack_trackerLayers[NTRACK];   //[nIsoTrack]  
  Int_t           IsoTrack_pixelLayers[NTRACK];   //[nIsoTrack]    
  Int_t           IsoTrack_trackerHits[NTRACK];   //[nIsoTrack]    
  Int_t           IsoTrack_pixelHits[NTRACK];   //[nIsoTrack]      
  Int_t           IsoTrack_missingInnerPixelHits[NTRACK];   //[nIsoTrack]                                      
  Int_t           IsoTrack_missingOuterPixelHits[NTRACK];   //[nIsoTrack]                                      
  Int_t           IsoTrack_missingInnerStripHits[NTRACK];   //[nIsoTrack]                                      
  Int_t           IsoTrack_missingOuterStripHits[NTRACK];   //[nIsoTrack]                                      
  Int_t           IsoTrack_missingInnerTrackerHits[NTRACK];   //[nIsoTrack]                                    
  Int_t           IsoTrack_missingOuterTrackerHits[NTRACK];   //[nIsoTrack]                                    
  Int_t           IsoTrack_missingMiddleTrackerHits[NTRACK];   //[nIsoTrack]                                   
  //
  Float_t         IsoTrack_dedxByHit0[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit1[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit2[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit3[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit4[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit5[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit6[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit7[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit8[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit9[NTRACK];   //[nIsoTrack]     
  Float_t         IsoTrack_dedxByHit10[NTRACK];   //[nIsoTrack]    
  Float_t         IsoTrack_dedxByHit11[NTRACK];   //[nIsoTrack]    
  Float_t         IsoTrack_dedxByHit12[NTRACK];   //[nIsoTrack]    
  Float_t         IsoTrack_dedxByHit13[NTRACK];   //[nIsoTrack]    
  Float_t         IsoTrack_dedxUnSmearedByHit0[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit1[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit2[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit3[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit4[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit5[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit6[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit7[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit8[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit9[NTRACK];   //[nIsoTrack]                                        
  Float_t         IsoTrack_dedxUnSmearedByHit10[NTRACK];   //[nIsoTrack]                                       
  Float_t         IsoTrack_dedxUnSmearedByHit11[NTRACK];   //[nIsoTrack]                                       
  Float_t         IsoTrack_dedxUnSmearedByHit12[NTRACK];   //[nIsoTrack]                                       
  Float_t         IsoTrack_dedxUnSmearedByHit13[NTRACK];   //[nIsoTrack]                                       
  //

  // Set branch addresses //
  fChain->SetBranchAddress("run", &run); // , &b_run);
  fChain->SetBranchAddress("lumi", &lumi); // , &b_lumi);
  fChain->SetBranchAddress("evt", &evt); // , &b_evt);
  fChain->SetBranchAddress("isData", &isData); // , &b_isData);
  //
  fChain->SetBranchAddress("xsec", &xsec); // , &b_xsec);
  fChain->SetBranchAddress("puWeight", &puWeight); // , &b_puWeight);
  fChain->SetBranchAddress("nTrueInt", &nTrueInt); // , &b_nTrueInt);
  fChain->SetBranchAddress("genWeight", &genWeight); // , &b_genWeight);
  fChain->SetBranchAddress("nVert", &nVert); // , &b_nVert);
  //
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60); // , &b_HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60); // , &b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight); // , &b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight); // , &b_HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
  fChain->SetBranchAddress("HLT_MET", &HLT_MET); // , &b_HLT_MET);
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v", &HLT_BIT_HLT_IsoMu24_v); // , &b_HLT_BIT_HLT_IsoMu24_v);
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_v", &HLT_BIT_HLT_IsoMu24_eta2p1_v); // , &b_HLT_BIT_HLT_IsoMu24_eta2p1_v);
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v", &HLT_BIT_HLT_IsoMu27_v); // , &b_HLT_BIT_HLT_IsoMu27_v);
  fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu); // , &b_HLT_SingleMu);
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_WPTight_Gsf_v", &HLT_BIT_HLT_Ele32_WPTight_Gsf_v); // , &b_HLT_BIT_HLT_Ele32_WPTight_Gsf_v);
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele35_WPTight_Gsf_v", &HLT_BIT_HLT_Ele35_WPTight_Gsf_v); // , &b_HLT_BIT_HLT_Ele35_WPTight_Gsf_v);
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v", &HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v); // , &b_HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
  fChain->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl); // , &b_HLT_SingleEl);
  //
  fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices); // , &b_Flag_goodVertices);
  fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter); // , &b_Flag_BadPFMuonFilter);
  fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter); // , &b_Flag_HBHENoiseIsoFilter);
  fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // , &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
  fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // , &b_Flag_eeBadScFilter);
  fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter); // , &b_Flag_BadChargedCandidateFilter);
  fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter); // , &b_Flag_ecalBadCalibFilter);
  fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // , &b_Flag_HBHENoiseFilter);
  fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter); // , &b_Flag_globalTightHalo2016Filter);
  //
  fChain->SetBranchAddress("met_sumEt", &met_sumEt); // , &b_met_sumEt);
  fChain->SetBranchAddress("met_rawPt", &met_rawPt); // , &b_met_rawPt);
  fChain->SetBranchAddress("met_rawPhi", &met_rawPhi); // , &b_met_rawPhi);
  fChain->SetBranchAddress("met_rawSumEt", &met_rawSumEt); // , &b_met_rawSumEt);
  fChain->SetBranchAddress("met_genPt", &met_genPt); // , &b_met_genPt);
  fChain->SetBranchAddress("met_genPhi", &met_genPhi); // , &b_met_genPhi);
  fChain->SetBranchAddress("met_genEta", &met_genEta); // , &b_met_genEta);
  fChain->SetBranchAddress("met_pt", &met_pt); // , &b_met_pt);
  fChain->SetBranchAddress("met_eta", &met_eta); // , &b_met_eta);
  fChain->SetBranchAddress("met_phi", &met_phi); // , &b_met_phi);
  fChain->SetBranchAddress("met_mass", &met_mass); // , &b_met_mass);
  fChain->SetBranchAddress("metNoMu_pt", &metNoMu_pt); // , &b_metNoMu_pt);
  fChain->SetBranchAddress("metNoMu_eta", &metNoMu_eta); // , &b_metNoMu_eta);
  fChain->SetBranchAddress("metNoMu_phi", &metNoMu_phi); // , &b_metNoMu_phi);
  fChain->SetBranchAddress("metNoMu_mass", &metNoMu_mass); // , &b_metNoMu_mass);
  //
  fChain->SetBranchAddress("nJet", &nJet); // , &b_nJet);
  fChain->SetBranchAddress("nJet30", &nJet30); // , &b_nJet30);
  fChain->SetBranchAddress("nJet30a", &nJet30a); // , &b_nJet30a);
  fChain->SetBranchAddress("Jet_pt", Jet_pt); // , &b_Jet_pt);
  fChain->SetBranchAddress("Jet_eta", Jet_eta); // , &b_Jet_eta);
  fChain->SetBranchAddress("Jet_phi", Jet_phi); // , &b_Jet_phi);
  fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF); // , &b_Jet_chHEF);
  fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF); // , &b_Jet_neHEF);
  //
  fChain->SetBranchAddress("nJetFwd", &nJetFwd); // , &b_nJetFwd);
  fChain->SetBranchAddress("JetFwd_pt", JetFwd_pt); // , &b_JetFwd_pt);
  fChain->SetBranchAddress("JetFwd_eta", JetFwd_eta); // , &b_JetFwd_eta);
  fChain->SetBranchAddress("JetFwd_phi", JetFwd_phi); // , &b_JetFwd_phi);
  fChain->SetBranchAddress("JetFwd_chHEF", JetFwd_chHEF); // , &b_JetFwd_chHEF);
  fChain->SetBranchAddress("JetFwd_neHEF", JetFwd_neHEF); // , &b_JetFwd_neHEF);
  //
  fChain->SetBranchAddress("nTauGood", &nTauGood); // , &b_nTauGood);
  fChain->SetBranchAddress("nLepGood", &nLepGood); // , &b_nLepGood);
  //
  fChain->SetBranchAddress("nIsoTrack", &nIsoTrack); // , &b_nIsoTrack);
  //
  fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt); // , &b_IsoTrack_pt);
  fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta); // , &b_IsoTrack_eta);
  fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi); // , &b_IsoTrack_phi);
  fChain->SetBranchAddress("IsoTrack_mass", IsoTrack_mass); // , &b_IsoTrack_mass);
  //
  fChain->SetBranchAddress("IsoTrack_charge", IsoTrack_charge); // , &b_IsoTrack_charge);
  fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy); // , &b_IsoTrack_dxy);
  fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz); // , &b_IsoTrack_dz);
  fChain->SetBranchAddress("IsoTrack_edxy", IsoTrack_edxy); // , &b_IsoTrack_edxy);
  fChain->SetBranchAddress("IsoTrack_edz", IsoTrack_edz); // , &b_IsoTrack_edz);
  //
  fChain->SetBranchAddress("IsoTrack_miniRelIso", IsoTrack_miniRelIso); // , &b_IsoTrack_miniRelIso);
  fChain->SetBranchAddress("IsoTrack_relIso03", IsoTrack_relIso03); // , &b_IsoTrack_relIso03);
  //
  fChain->SetBranchAddress("IsoTrack_trackerLayers", IsoTrack_trackerLayers); // , &b_IsoTrack_trackerLayers);
  fChain->SetBranchAddress("IsoTrack_pixelLayers", IsoTrack_pixelLayers); // , &b_IsoTrack_pixelLayers);
  fChain->SetBranchAddress("IsoTrack_trackerHits", IsoTrack_trackerHits); // , &b_IsoTrack_trackerHits);
  fChain->SetBranchAddress("IsoTrack_pixelHits", IsoTrack_pixelHits); // , &b_IsoTrack_pixelHits);
  fChain->SetBranchAddress("IsoTrack_missingInnerPixelHits", IsoTrack_missingInnerPixelHits); // , &b_IsoTrack_missingInnerPixelHits);
  fChain->SetBranchAddress("IsoTrack_missingOuterPixelHits", IsoTrack_missingOuterPixelHits); // , &b_IsoTrack_missingOuterPixelHits);
  fChain->SetBranchAddress("IsoTrack_missingInnerStripHits", IsoTrack_missingInnerStripHits); // , &b_IsoTrack_missingInnerStripHits);
  fChain->SetBranchAddress("IsoTrack_missingOuterStripHits", IsoTrack_missingOuterStripHits); // , &b_IsoTrack_missingOuterStripHits);
  fChain->SetBranchAddress("IsoTrack_missingInnerTrackerHits", IsoTrack_missingInnerTrackerHits); // , &b_IsoTrack_missingInnerTrackerHits);
  fChain->SetBranchAddress("IsoTrack_missingOuterTrackerHits", IsoTrack_missingOuterTrackerHits); // , &b_IsoTrack_missingOuterTrackerHits);
  fChain->SetBranchAddress("IsoTrack_missingMiddleTrackerHits", IsoTrack_missingMiddleTrackerHits); // , &b_IsoTrack_missingMiddleTrackerHits);
  fChain->SetBranchAddress("IsoTrack_dedxByHit0", IsoTrack_dedxByHit0); // , &b_IsoTrack_dedxByHit0);
  fChain->SetBranchAddress("IsoTrack_dedxByHit1", IsoTrack_dedxByHit1); // , &b_IsoTrack_dedxByHit1);
  fChain->SetBranchAddress("IsoTrack_dedxByHit2", IsoTrack_dedxByHit2); // , &b_IsoTrack_dedxByHit2);
  fChain->SetBranchAddress("IsoTrack_dedxByHit3", IsoTrack_dedxByHit3); // , &b_IsoTrack_dedxByHit3);
  fChain->SetBranchAddress("IsoTrack_dedxByHit4", IsoTrack_dedxByHit4); // , &b_IsoTrack_dedxByHit4);
  fChain->SetBranchAddress("IsoTrack_dedxByHit5", IsoTrack_dedxByHit5); // , &b_IsoTrack_dedxByHit5);
  fChain->SetBranchAddress("IsoTrack_dedxByHit6", IsoTrack_dedxByHit6); // , &b_IsoTrack_dedxByHit6);
  fChain->SetBranchAddress("IsoTrack_dedxByHit7", IsoTrack_dedxByHit7); // , &b_IsoTrack_dedxByHit7);
  fChain->SetBranchAddress("IsoTrack_dedxByHit8", IsoTrack_dedxByHit8); // , &b_IsoTrack_dedxByHit8);
  fChain->SetBranchAddress("IsoTrack_dedxByHit9", IsoTrack_dedxByHit9); // , &b_IsoTrack_dedxByHit9);
  fChain->SetBranchAddress("IsoTrack_dedxByHit10", IsoTrack_dedxByHit10); // , &b_IsoTrack_dedxByHit10);
  fChain->SetBranchAddress("IsoTrack_dedxByHit11", IsoTrack_dedxByHit11); // , &b_IsoTrack_dedxByHit11);
  fChain->SetBranchAddress("IsoTrack_dedxByHit12", IsoTrack_dedxByHit12); // , &b_IsoTrack_dedxByHit12);
  fChain->SetBranchAddress("IsoTrack_dedxByHit13", IsoTrack_dedxByHit13); // , &b_IsoTrack_dedxByHit13);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit0", IsoTrack_dedxUnSmearedByHit0); // , &b_IsoTrack_dedxUnSmearedByHit0);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit1", IsoTrack_dedxUnSmearedByHit1); // , &b_IsoTrack_dedxUnSmearedByHit1);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit2", IsoTrack_dedxUnSmearedByHit2); // , &b_IsoTrack_dedxUnSmearedByHit2);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit3", IsoTrack_dedxUnSmearedByHit3); // , &b_IsoTrack_dedxUnSmearedByHit3);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit4", IsoTrack_dedxUnSmearedByHit4); // , &b_IsoTrack_dedxUnSmearedByHit4);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit5", IsoTrack_dedxUnSmearedByHit5); // , &b_IsoTrack_dedxUnSmearedByHit5);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit6", IsoTrack_dedxUnSmearedByHit6); // , &b_IsoTrack_dedxUnSmearedByHit6);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit7", IsoTrack_dedxUnSmearedByHit7); // , &b_IsoTrack_dedxUnSmearedByHit7);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit8", IsoTrack_dedxUnSmearedByHit8); // , &b_IsoTrack_dedxUnSmearedByHit8);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit9", IsoTrack_dedxUnSmearedByHit9); // , &b_IsoTrack_dedxUnSmearedByHit9);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit10", IsoTrack_dedxUnSmearedByHit10); // , &b_IsoTrack_dedxUnSmearedByHit10);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit11", IsoTrack_dedxUnSmearedByHit11); // , &b_IsoTrack_dedxUnSmearedByHit11);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit12", IsoTrack_dedxUnSmearedByHit12); // , &b_IsoTrack_dedxUnSmearedByHit12);
  fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit13", IsoTrack_dedxUnSmearedByHit13); // , &b_IsoTrack_dedxUnSmearedByHit13);
  //

  // Utility variables //
  UInt_t  iMax=0;
  Float_t ptMax=0;
  Float_t dphi_jme=-999;

  // Counters
  const UInt_t nEVC = 9;
  UInt_t  nEventCount[nEVC];
  for(UInt_t i=0; i<nEVC; i++) nEventCount[i]=0;
  TString nameStep[nEVC] = { "MET Trigger" , "MET Flags" , "Lepton Veto" , "Tau Veto" , 
			     "MET Cut" , "#Jets" , "#Tracks" , "Leading Jet" , "DeltaPhi(Jet,MET)" }; 

  Int_t   foundPathologicalJet = -1;
  Bool_t leadJetIsFwd = kFALSE;
  Bool_t passDeltaPhi = kTRUE;

  // < <index in array, isForward>
  vector< pair<UInt_t,Bool_t> > idxJetID;
  vector<UInt_t> idxTrackID;

  // Loop over the chain //
  UInt_t nEntries = fChain->GetEntries();
  cout << "- About to process: " << nEntries << " entries." << endl;
  for(UInt_t iE=0; iE<nEntries; iE++) {
    
    fChain->GetEntry(iE);


    // 1- Loop over jets to find the leading one and the ID'ed ones

    iMax=0;
    ptMax=0;
    leadJetIsFwd=kFALSE;
    dphi_jme=-999;
    foundPathologicalJet = -1;
    idxJetID.clear();

    for(UInt_t iJ=0; iJ<NJET && (Int_t)iJ<nJet ; iJ++) {

      // Slim down jet collection
      if( Jet_pt[iJ]>30 
	  && Jet_chHEF[iJ]>0.01 && Jet_chHEF[iJ]<0.99
	  && Jet_neHEF[iJ]>0.01 && Jet_neHEF[iJ]<0.99 ) {
	
	idxJetID.push_back(make_pair(iJ,kFALSE));

	// Search for leading jets only among the ID'ed ones
	if(Jet_pt[iJ] > ptMax) {
	  ptMax = Jet_pt[iJ];
	  iMax  = iJ;
	}
      }

      // Check for pathological jets (array issues)
      if( Jet_pt[iJ]>1e+04 ) {
	foundPathologicalJet = iJ;
	cout << "-- Just found a pathological jet: Jet_pt[" << iJ << "] = " << Jet_pt[iJ] << endl;
      }

    } // end loop: jets

    /// Loop over forward jets
    for(UInt_t iJ=0; iJ<NJETF && (Int_t)iJ<nJetFwd ; iJ++) {

      // Slim down jet collection
      if( JetFwd_pt[iJ]>30 
	  && JetFwd_chHEF[iJ]>0.01 && JetFwd_chHEF[iJ]<0.99
	  && JetFwd_neHEF[iJ]>0.01 && JetFwd_neHEF[iJ]<0.99 ) {
	
	idxJetID.push_back(make_pair(iJ,kTRUE));

	// Search for leading jets only among the ID'ed ones
	if(JetFwd_pt[iJ] > ptMax) {
	  ptMax = JetFwd_pt[iJ];
	  iMax  = iJ;
	  leadJetIsFwd = kTRUE;
	}
      }

      // Check for pathological jets (array issues)
      if( JetFwd_pt[iJ]>1e+04 ) {
	foundPathologicalJet = iJ;
	cout << "-- Just found a pathological jet: JetFwd_pt[" << iJ << "] = " << JetFwd_pt[iJ] << endl;
      }

    } // end loop: forward jets

    
    // 2- Print-out leading-jet search
    if( doPrintBeforeCuts && (iE%10000==0 || iMax>=NJET || foundPathologicalJet>=0) ) {

      cout << "----------------------------------------" 
	   << endl
	   << "ptMax = " << ptMax << " | iMax = " << iMax ;

      if(iMax<NJET) cout << " | Jet_pt[iMax] = " << Jet_pt[iMax];
      else          cout << " WRONG iMax=" << iMax ;

      cout << endl
	   << "nJet = " << nJet 
	   << endl << endl;

      if(Jet_pt[1]>Jet_pt[0])
	cout << Jet_pt[0] << " " << Jet_pt[1] << endl << endl;
      
      if(foundPathologicalJet>=0) {
	cout << "Found pathological jet: Jet_pt[" << foundPathologicalJet  << "] = ";
	if(foundPathologicalJet<(Int_t)NJET) cout << Jet_pt[foundPathologicalJet] << endl;
	else                                 cout << "foundPathologicalJet>NJET"  << endl;
      } 
    } // end print-out leading jet search


    // 3- Track selection //
    idxTrackID.clear();
    for(UInt_t iT=0; iT<NTRACK && (Int_t)iT<nIsoTrack ; iT++) {
      if( TMath::Abs(IsoTrack_eta[iT]) < 2.1
	  && IsoTrack_relIso03[iT]     < 0.5 // confirmed
	  //&& IsoTrack_miniRelIso[iT]   < 0.5
	  && IsoTrack_missingInnerPixelHits[iT]    == 0
	  && IsoTrack_missingMiddleTrackerHits[iT] == 0 ) {
	  //&& IsoTrack_trackerLayers[iT]            >= 2 ) { // fixme: ask Jeremi // does not change event count

	idxTrackID.push_back(iT);
      }
    }

    // 4- Level-1 Cuts //

    //// MET trigger and flags
    //if(!HLT_MET)          continue; // fixme: ask Jeremi
    if(!HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight) continue;
    nEventCount[0]++ ;

    if(doCutOnMetFlags) {
      if(!Flag_goodVertices) continue;
      if(!Flag_BadPFMuonFilter)    continue;
      if(!Flag_HBHENoiseIsoFilter) continue;
      if(!Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
      if(!Flag_eeBadScFilter)             continue; // not suggested for MC, see twiki MissingETOptionalFiltersRun2
      //if(!Flag_BadChargedCandidateFilter) continue; // not in L1 cuts, not recommended (under review), see twiki
      if(!Flag_ecalBadCalibFilter)        continue;
      if(!Flag_HBHENoiseFilter)           continue;
      if(!Flag_globalTightHalo2016Filter) continue;
    }
    nEventCount[1]++ ;
    //// 

    //// Lepton veto
    if(nLepGood>0) continue; 
    nEventCount[2]++ ;
    if(nTauGood>0) continue;
    nEventCount[3]++ ;

    //// MET Cut
    if(metNoMu_pt<=200)   continue;
    nEventCount[4]++ ;

    //// Number of jets passing ID
    if(idxJetID.size()==0) continue;
    nEventCount[5]++ ;

    //// Number of tracks passing ID
    if(idxTrackID.size()==0) continue;
    nEventCount[6]++ ;

    //// Intermediate printout
    outlog << "Run: " << run << " Lumi: " << lumi << " Event: " << evt << endl << endl;

    ///// Loop over preselected central jets
    for(UInt_t idxJ=0; idxJ<idxJetID.size(); idxJ++) {

      // break this loop if we reached forward jets
      if(idxJetID[idxJ].second) break; 

      UInt_t iJ = idxJetID[idxJ].first;

      outlog << "pt: "   << Jet_pt[iJ]    << "\teta: "      << Jet_eta[iJ]   << endl
	     << "phi: "  << Jet_phi[iJ]   << "\tf_CH: "     << Jet_chHEF[iJ] << endl
	     << "f_NE: " << Jet_neHEF[iJ] << "\tforward: 0" << endl;

      if(iJ==iMax) outlog << "=> LEADING JET" << endl;

      if(Jet_pt[iJ]<=100)                outlog << "=> Fails cut: pT>100"    << endl;
      if(TMath::Abs(Jet_eta[iJ]) >= 2.4) outlog << "=> Fails cut: |eta|<2.4" << endl;

      if(Jet_neHEF[iJ]>=0.8) outlog << "=> Fails cut: NHEF<0.8" << endl;
      if(Jet_chHEF[iJ]<=0.1) outlog << "=> Fails cut: CHEF>0.1" << endl;

      outlog << endl;
    }

    ///// Loop over preselected forward jets
    for(UInt_t idxJ=0; idxJ<idxJetID.size(); idxJ++) {

      // go beyond central jets to reach forward jets
      if(!idxJetID[idxJ].second) continue; 

      UInt_t iJ = idxJetID[idxJ].first;

      outlog << "pt: "   << JetFwd_pt[iJ]    << "\teta: "      << JetFwd_eta[iJ]   << endl
	     << "phi: "  << JetFwd_phi[iJ]   << "\tf_CH: "     << JetFwd_chHEF[iJ] << endl
	     << "f_NE: " << JetFwd_neHEF[iJ] << "\tforward: 1" << endl;

      if(iJ==iMax) outlog << "=> LEADING JET" << endl;

      if(JetFwd_pt[iJ]<=100)                outlog << "=> Fails cut: pT>100"    << endl;
      if(TMath::Abs(JetFwd_eta[iJ]) >= 2.4) outlog << "=> Fails cut: |eta|<2.4" << endl;

      if(JetFwd_neHEF[iJ]>=0.8) outlog << "=> Fails cut: NHEF<0.8" << endl;
      if(JetFwd_chHEF[iJ]<=0.1) outlog << "=> Fails cut: CHEF>0.1" << endl;

      outlog << endl;
    }

    ///// MET Flags
    if(!Flag_goodVertices) outlog << "Fails MET Flag_goodVertices" << endl;
    if(!Flag_BadPFMuonFilter)    outlog << "Fails MET Flag_BadPFMuonFilter" << endl;
    if(!Flag_HBHENoiseIsoFilter) outlog << "Fails MET Flag_HBHENoiseIsoFilter" << endl;
    if(!Flag_EcalDeadCellTriggerPrimitiveFilter) outlog << "Fails MET Flag_EcalDeadCellTriggerPrimitiveFilter" << endl;
    if(!Flag_eeBadScFilter)             outlog << "Fails MET Flag_eeBadScFilter" << endl; 
    if(!Flag_ecalBadCalibFilter)        outlog << "Fails MET Flag_ecalBadCalibFilter" << endl;
    if(!Flag_HBHENoiseFilter)           outlog << "Fails MET Flag_HBHENoiseFilter" << endl;
    if(!Flag_globalTightHalo2016Filter) outlog << "Fails MET Flag_globalTightHalo2016Filter" << endl;
    //
    outlog << endl << endl;

    //// Leading jet
    if(leadJetIsFwd) continue;
    else {
      if(Jet_pt[iMax]<=100)    continue;
      if(TMath::Abs(Jet_eta[iMax]) >= 2.4) continue;
      //
      if(Jet_neHEF[iMax]>=0.8) continue;
      if(Jet_chHEF[iMax]<=0.1) continue;
      //
      if(Jet_neHEF[iMax]<=0.01) continue; //sanity check
      if(Jet_chHEF[iMax]>=0.99) continue; //sanity check
    }
    //
    nEventCount[7]++ ;
    ////

    //// QCD suppression // fixme: use MET or METNoMu ?
    passDeltaPhi = kTRUE;
    for(UInt_t idxJ=0; idxJ<idxJetID.size(); idxJ++) { 
      UInt_t iJ = idxJetID[idxJ].first;

      if(idxJetID[idxJ].second) // forward jet
	dphi_jme = TVector2::Phi_mpi_pi( JetFwd_phi[iJ] - metNoMu_phi);
      else // central jet
	dphi_jme = TVector2::Phi_mpi_pi( Jet_phi[iJ] - metNoMu_phi);

      if(TMath::Abs(dphi_jme) <= 0.5) passDeltaPhi = kFALSE;
    }
    if(!passDeltaPhi) continue;
    nEventCount[8]++ ;
    ////

    //// Map selected events up to the previous step
    mapEvents[lumi].push_back(evt);

    // Printouts after cuts //
    for(UInt_t iJ=0; iJ<NJET; iJ++) {
      if(Jet_pt[iJ]>2000) {
	cout << "Entry #" << iE << " Jet_pt[" << iJ << "] = " << Jet_pt[iJ] << endl;
      }
    }

  } // end loop: entries

  // Printout # accepted events
  cout << endl << "----------------" << endl << "ACCEPTED EVENTS" << endl;
  for(UInt_t i=0; i<nEVC; i++) 
    cout << "#Events passing step #" << i << " = " << nEventCount[i] 
	 << "\t(" << nameStep[i] << ")" << endl;

  // Printout # rejected events
  cout << endl << "----------------" << endl << "REJECTED EVENTS" << endl;
  for(UInt_t i=0; i<nEVC; i++) 
    cout << "#Events rejected by step #" << i << " = " << nEntries - nEventCount[i] 
	 << "\t(" << nameStep[i] << ")" << endl;

  // Printout map of selected events
  lumi = 0;
  vEvents.clear();
  //
  for(itMapEvents = mapEvents.begin(); 
      itMapEvents != mapEvents.end() ; 
      itMapEvents++) {

    lumi    = itMapEvents->first;
    vEvents = itMapEvents->second;
    sort(vEvents.begin(), vEvents.end()); 

    for(UInt_t iE=0; iE<vEvents.size(); iE++) {
      outmap << lumi << ":" << vEvents[iE] << endl;
    }
  }


  return 0;
}
