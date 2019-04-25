//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr  4 17:24:06 2019 by ROOT version 6.10/09
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef XTrackNtuple_h
#define XTrackNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class XTrackNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       evt;
   Int_t           isData;
   Float_t         xsec;
   Float_t         puWeight;
   Float_t         nTrueInt;
   Float_t         genWeight;
   Float_t         rho;
   Int_t           nVert;
   Float_t         vertex_x;
   Float_t         vertex_y;
   Float_t         vertex_z;
   Int_t           nJet30;
   Int_t           nJet30a;
   Float_t         lheHT;
   Float_t         lheHTIncoming;
   Float_t         lheVpt;
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
   Int_t           Flag_goodVertices;
   Int_t           Flag_BadPFMuonFilter;
   Int_t           Flag_HBHENoiseIsoFilter;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_BadChargedCandidateFilter;
   Int_t           Flag_ecalBadCalibFilter;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_globalTightHalo2016Filter;
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
   Float_t         met_jecUp_sumEt;
   Float_t         met_jecUp_rawPt;
   Float_t         met_jecUp_rawPhi;
   Float_t         met_jecUp_rawSumEt;
   Float_t         met_jecUp_genPt;
   Float_t         met_jecUp_genPhi;
   Float_t         met_jecUp_genEta;
   Float_t         met_jecUp_pt;
   Float_t         met_jecUp_eta;
   Float_t         met_jecUp_phi;
   Float_t         met_jecUp_mass;
   Float_t         metNoMu_jecUp_pt;
   Float_t         metNoMu_jecUp_eta;
   Float_t         metNoMu_jecUp_phi;
   Float_t         metNoMu_jecUp_mass;
   Float_t         metNoMu_jecUp_pt;
   Float_t         metNoMu_jecUp_eta;
   Float_t         metNoMu_jecUp_phi;
   Float_t         metNoMu_jecUp_mass;
   Float_t         met_jecDown_sumEt;
   Float_t         met_jecDown_rawPt;
   Float_t         met_jecDown_rawPhi;
   Float_t         met_jecDown_rawSumEt;
   Float_t         met_jecDown_genPt;
   Float_t         met_jecDown_genPhi;
   Float_t         met_jecDown_genEta;
   Float_t         met_jecDown_pt;
   Float_t         met_jecDown_eta;
   Float_t         met_jecDown_phi;
   Float_t         met_jecDown_mass;
   Int_t           nGenChargino;
   Int_t           GenChargino_motherId[1];   //[nGenChargino]
   Int_t           GenChargino_grandmotherId[1];   //[nGenChargino]
   Float_t         GenChargino_charge[1];   //[nGenChargino]
   Int_t           GenChargino_status[1];   //[nGenChargino]
   Int_t           GenChargino_isPromptHard[1];   //[nGenChargino]
   Int_t           GenChargino_pdgId[1];   //[nGenChargino]
   Float_t         GenChargino_pt[1];   //[nGenChargino]
   Float_t         GenChargino_eta[1];   //[nGenChargino]
   Float_t         GenChargino_phi[1];   //[nGenChargino]
   Float_t         GenChargino_mass[1];   //[nGenChargino]
   Float_t         GenChargino_beta[1];   //[nGenChargino]
   Float_t         GenChargino_decayR[1];   //[nGenChargino]
   Float_t         GenChargino_decayZ[1];   //[nGenChargino]
   Int_t           nIsoTrack;
   Int_t           IsoTrack_pdgId[4];   //[nIsoTrack]
   Float_t         IsoTrack_pt[4];   //[nIsoTrack]
   Float_t         IsoTrack_eta[4];   //[nIsoTrack]
   Float_t         IsoTrack_phi[4];   //[nIsoTrack]
   Float_t         IsoTrack_mass[4];   //[nIsoTrack]
   Int_t           IsoTrack_trackerLayers[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixelLayers[4];   //[nIsoTrack]
   Int_t           IsoTrack_trackerHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixelHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingInnerPixelHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingOuterPixelHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingInnerStripHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingOuterStripHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingInnerTrackerHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingOuterTrackerHits[4];   //[nIsoTrack]
   Int_t           IsoTrack_missingMiddleTrackerHits[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit0[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit1[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit2[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit3[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit4[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit5[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit6[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit7[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit8[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit9[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit10[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit11[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit12[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxByHit13[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit0[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit1[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit2[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit3[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit4[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit5[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit6[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit7[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit8[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit9[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit10[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit11[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit12[4];   //[nIsoTrack]
   Float_t         IsoTrack_dedxUnSmearedByHit13[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit0[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit1[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit2[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit3[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit4[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit5[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit6[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit7[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit8[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit9[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit10[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit11[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit12[4];   //[nIsoTrack]
   Float_t         IsoTrack_deByHit13[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit0[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit1[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit2[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit3[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit4[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit5[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit6[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit7[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit8[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit9[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit10[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit11[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit12[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_subDetIdByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeXbyHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_sizeYbyHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_layerOrSideByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_ladderOrBladeByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_sideByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit0[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit1[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit2[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit3[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit4[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit5[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit6[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit7[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit8[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit9[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit10[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit11[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit12[4];   //[nIsoTrack]
   Int_t           IsoTrack_pixByHit13[4];   //[nIsoTrack]
   Int_t           IsoTrack_charge[4];   //[nIsoTrack]
   Float_t         IsoTrack_dxy[4];   //[nIsoTrack]
   Float_t         IsoTrack_dz[4];   //[nIsoTrack]
   Float_t         IsoTrack_edxy[4];   //[nIsoTrack]
   Float_t         IsoTrack_edz[4];   //[nIsoTrack]
   Int_t           IsoTrack_highPurity[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniIsoCH[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniIsoNH[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniIsoPH[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniIsoPU[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniRelIso[4];   //[nIsoTrack]
   Float_t         IsoTrack_dR03IsoCH[4];   //[nIsoTrack]
   Float_t         IsoTrack_dR03IsoNH[4];   //[nIsoTrack]
   Float_t         IsoTrack_dR03IsoPH[4];   //[nIsoTrack]
   Float_t         IsoTrack_dR03IsoPU[4];   //[nIsoTrack]
   Float_t         IsoTrack_relIso03[4];   //[nIsoTrack]
   Float_t         IsoTrack_caloEmEnergy[4];   //[nIsoTrack]
   Float_t         IsoTrack_caloHadEnergy[4];   //[nIsoTrack]
   Int_t           IsoTrack_channelsGoodECAL[4];   //[nIsoTrack]
   Int_t           IsoTrack_channelsGoodHCAL[4];   //[nIsoTrack]
   Int_t           IsoTrack_awayJet_idx[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayJet_pt[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayJet_dr[4];   //[nIsoTrack]
   Int_t           IsoTrack_awayNJet[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayHTJet[4];   //[nIsoTrack]
   Int_t           IsoTrack_awayMu_idx[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayMu_dr[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayMu_mll[4];   //[nIsoTrack]
   Int_t           IsoTrack_awayEle_idx[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayEle_dr[4];   //[nIsoTrack]
   Float_t         IsoTrack_awayEle_mll[4];   //[nIsoTrack]
   Int_t           IsoTrack_closestMu_idx[4];   //[nIsoTrack]
   Int_t           IsoTrack_closestEle_idx[4];   //[nIsoTrack]
   Int_t           IsoTrack_closestTau_idx[4];   //[nIsoTrack]
   Int_t           IsoTrack_trigLepton_idx[4];   //[nIsoTrack]
   Float_t         IsoTrack_myDeDx[4];   //[nIsoTrack]
   Int_t           IsoTrack_mcMatch[4];   //[nIsoTrack]
   Int_t           IsoTrack_mcMatchAnyId[4];   //[nIsoTrack]
   Int_t           IsoTrack_mcMatchAnyPt[4];   //[nIsoTrack]
   Int_t           nJetFwd;
   Int_t           JetFwd_id[6];   //[nJetFwd]
   Int_t           JetFwd_puId[6];   //[nJetFwd]
   Float_t         JetFwd_btagCSV[6];   //[nJetFwd]
   Float_t         JetFwd_btagDeepCSV[6];   //[nJetFwd]
   Float_t         JetFwd_rawPt[6];   //[nJetFwd]
   Float_t         JetFwd_mcPt[6];   //[nJetFwd]
   Int_t           JetFwd_mcFlavour[6];   //[nJetFwd]
   Int_t           JetFwd_hadronFlavour[6];   //[nJetFwd]
   Int_t           JetFwd_mcMatchId[6];   //[nJetFwd]
   Float_t         JetFwd_corr_JECUp[6];   //[nJetFwd]
   Float_t         JetFwd_corr_JECDown[6];   //[nJetFwd]
   Float_t         JetFwd_corr[6];   //[nJetFwd]
   Float_t         JetFwd_corr_JERUp[6];   //[nJetFwd]
   Float_t         JetFwd_corr_JERDown[6];   //[nJetFwd]
   Float_t         JetFwd_corr_JER[6];   //[nJetFwd]
   Float_t         JetFwd_pt[6];   //[nJetFwd]
   Float_t         JetFwd_eta[6];   //[nJetFwd]
   Float_t         JetFwd_phi[6];   //[nJetFwd]
   Float_t         JetFwd_mass[6];   //[nJetFwd]
   Float_t         JetFwd_chHEF[6];   //[nJetFwd]
   Float_t         JetFwd_neHEF[6];   //[nJetFwd]
   Int_t           nJet;
   Float_t         Jet_area[15];   //[nJet]
   Float_t         Jet_qgl[15];   //[nJet]
   Float_t         Jet_ptd[15];   //[nJet]
   Float_t         Jet_axis2[15];   //[nJet]
   Float_t         Jet_axis1[15];   //[nJet]
   Int_t           Jet_mult[15];   //[nJet]
   Int_t           Jet_partonId[15];   //[nJet]
   Int_t           Jet_partonMotherId[15];   //[nJet]
   Float_t         Jet_nLeptons[15];   //[nJet]
   Int_t           Jet_id[15];   //[nJet]
   Int_t           Jet_puId[15];   //[nJet]
   Float_t         Jet_btagCSV[15];   //[nJet]
   Float_t         Jet_btagDeepCSV[15];   //[nJet]
   Float_t         Jet_rawPt[15];   //[nJet]
   Float_t         Jet_mcPt[15];   //[nJet]
   Int_t           Jet_mcFlavour[15];   //[nJet]
   Int_t           Jet_hadronFlavour[15];   //[nJet]
   Int_t           Jet_mcMatchId[15];   //[nJet]
   Float_t         Jet_corr_JECUp[15];   //[nJet]
   Float_t         Jet_corr_JECDown[15];   //[nJet]
   Float_t         Jet_corr[15];   //[nJet]
   Float_t         Jet_corr_JERUp[15];   //[nJet]
   Float_t         Jet_corr_JERDown[15];   //[nJet]
   Float_t         Jet_corr_JER[15];   //[nJet]
   Float_t         Jet_pt[15];   //[nJet]
   Float_t         Jet_eta[15];   //[nJet]
   Float_t         Jet_phi[15];   //[nJet]
   Float_t         Jet_mass[15];   //[nJet]
   Float_t         Jet_chHEF[15];   //[nJet]
   Float_t         Jet_neHEF[15];   //[nJet]
   Int_t           nTauGood;
   Int_t           TauGood_charge[2];   //[nTauGood]
   Int_t           TauGood_decayMode[2];   //[nTauGood]
   Int_t           TauGood_idDecayMode[2];   //[nTauGood]
   Int_t           TauGood_idDecayModeNewDMs[2];   //[nTauGood]
   Float_t         TauGood_dxy[2];   //[nTauGood]
   Float_t         TauGood_dz[2];   //[nTauGood]
   Int_t           TauGood_idMVA[2];   //[nTauGood]
   Int_t           TauGood_idMVANewDM[2];   //[nTauGood]
   Int_t           TauGood_idCI3hit[2];   //[nTauGood]
   Int_t           TauGood_idAntiMu[2];   //[nTauGood]
   Int_t           TauGood_idAntiE[2];   //[nTauGood]
   Float_t         TauGood_isoCI3hit[2];   //[nTauGood]
   Int_t           TauGood_mcMatchId[2];   //[nTauGood]
   Int_t           TauGood_pdgId[2];   //[nTauGood]
   Float_t         TauGood_pt[2];   //[nTauGood]
   Float_t         TauGood_eta[2];   //[nTauGood]
   Float_t         TauGood_phi[2];   //[nTauGood]
   Float_t         TauGood_mass[2];   //[nTauGood]
   Int_t           TauGood_idMVAdR03[2];   //[nTauGood]
   Int_t           nLepGood;
   Int_t           LepGood_charge[3];   //[nLepGood]
   Int_t           LepGood_tightId[3];   //[nLepGood]
   Int_t           LepGood_hltId[3];   //[nLepGood]
   Int_t           LepGood_eleCutIdSpring15_25ns_v1[3];   //[nLepGood]
   Int_t           LepGood_hasGainSwitchFlag[3];   //[nLepGood]
   Float_t         LepGood_dxy[3];   //[nLepGood]
   Float_t         LepGood_dz[3];   //[nLepGood]
   Float_t         LepGood_edxy[3];   //[nLepGood]
   Float_t         LepGood_edz[3];   //[nLepGood]
   Float_t         LepGood_ip3d[3];   //[nLepGood]
   Float_t         LepGood_sip3d[3];   //[nLepGood]
   Int_t           LepGood_convVeto[3];   //[nLepGood]
   Int_t           LepGood_lostHits[3];   //[nLepGood]
   Float_t         LepGood_relIso03[3];   //[nLepGood]
   Float_t         LepGood_relIso04[3];   //[nLepGood]
   Float_t         LepGood_miniRelIso[3];   //[nLepGood]
   Float_t         LepGood_relIsoAn04[3];   //[nLepGood]
   Int_t           LepGood_tightCharge[3];   //[nLepGood]
   Int_t           LepGood_mcMatchId[3];   //[nLepGood]
   Int_t           LepGood_mcMatchAny[3];   //[nLepGood]
   Int_t           LepGood_mcMatchTau[3];   //[nLepGood]
   Float_t         LepGood_mcPt[3];   //[nLepGood]
   Int_t           LepGood_mediumMuonId[3];   //[nLepGood]
   Int_t           LepGood_ICHEPsoftMuonId[3];   //[nLepGood]
   Int_t           LepGood_ICHEPmediumMuonId[3];   //[nLepGood]
   Int_t           LepGood_pdgId[3];   //[nLepGood]
   Float_t         LepGood_pt[3];   //[nLepGood]
   Float_t         LepGood_eta[3];   //[nLepGood]
   Float_t         LepGood_phi[3];   //[nLepGood]
   Float_t         LepGood_mass[3];   //[nLepGood]
   Float_t         LepGood_etaSc[3];   //[nLepGood]
   Int_t           LepGood_mcMatchPdgId[3];   //[nLepGood]
   Int_t           LepGood_mcPromptGamma[3];   //[nLepGood]
   Int_t           LepGood_triggerMatched[3];   //[nLepGood]
   Float_t         wgtsum;
   Float_t         kfact;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_puWeight;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_nJet30;   //!
   TBranch        *b_nJet30a;   //!
   TBranch        *b_lheHT;   //!
   TBranch        *b_lheHTIncoming;   //!
   TBranch        *b_lheVpt;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_MET;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu24_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu24_eta2p1_v;   //!
   TBranch        *b_HLT_BIT_HLT_IsoMu27_v;   //!
   TBranch        *b_HLT_SingleMu;   //!
   TBranch        *b_HLT_BIT_HLT_Ele32_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele35_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v;   //!
   TBranch        *b_HLT_SingleEl;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_met_sumEt;   //!
   TBranch        *b_met_rawPt;   //!
   TBranch        *b_met_rawPhi;   //!
   TBranch        *b_met_rawSumEt;   //!
   TBranch        *b_met_genPt;   //!
   TBranch        *b_met_genPhi;   //!
   TBranch        *b_met_genEta;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_mass;   //!
   TBranch        *b_metNoMu_pt;   //!
   TBranch        *b_metNoMu_eta;   //!
   TBranch        *b_metNoMu_phi;   //!
   TBranch        *b_metNoMu_mass;   //!
   TBranch        *b_met_jecUp_sumEt;   //!
   TBranch        *b_met_jecUp_rawPt;   //!
   TBranch        *b_met_jecUp_rawPhi;   //!
   TBranch        *b_met_jecUp_rawSumEt;   //!
   TBranch        *b_met_jecUp_genPt;   //!
   TBranch        *b_met_jecUp_genPhi;   //!
   TBranch        *b_met_jecUp_genEta;   //!
   TBranch        *b_met_jecUp_pt;   //!
   TBranch        *b_met_jecUp_eta;   //!
   TBranch        *b_met_jecUp_phi;   //!
   TBranch        *b_met_jecUp_mass;   //!
   TBranch        *b_metNoMu_jecUp_pt;   //!
   TBranch        *b_metNoMu_jecUp_eta;   //!
   TBranch        *b_metNoMu_jecUp_phi;   //!
   TBranch        *b_metNoMu_jecUp_mass;   //!
   TBranch        *b_metNoMu_jecUp_pt;   //!
   TBranch        *b_metNoMu_jecUp_eta;   //!
   TBranch        *b_metNoMu_jecUp_phi;   //!
   TBranch        *b_metNoMu_jecUp_mass;   //!
   TBranch        *b_met_jecDown_sumEt;   //!
   TBranch        *b_met_jecDown_rawPt;   //!
   TBranch        *b_met_jecDown_rawPhi;   //!
   TBranch        *b_met_jecDown_rawSumEt;   //!
   TBranch        *b_met_jecDown_genPt;   //!
   TBranch        *b_met_jecDown_genPhi;   //!
   TBranch        *b_met_jecDown_genEta;   //!
   TBranch        *b_met_jecDown_pt;   //!
   TBranch        *b_met_jecDown_eta;   //!
   TBranch        *b_met_jecDown_phi;   //!
   TBranch        *b_met_jecDown_mass;   //!
   TBranch        *b_nGenChargino;   //!
   TBranch        *b_GenChargino_motherId;   //!
   TBranch        *b_GenChargino_grandmotherId;   //!
   TBranch        *b_GenChargino_charge;   //!
   TBranch        *b_GenChargino_status;   //!
   TBranch        *b_GenChargino_isPromptHard;   //!
   TBranch        *b_GenChargino_pdgId;   //!
   TBranch        *b_GenChargino_pt;   //!
   TBranch        *b_GenChargino_eta;   //!
   TBranch        *b_GenChargino_phi;   //!
   TBranch        *b_GenChargino_mass;   //!
   TBranch        *b_GenChargino_beta;   //!
   TBranch        *b_GenChargino_decayR;   //!
   TBranch        *b_GenChargino_decayZ;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_mass;   //!
   TBranch        *b_IsoTrack_trackerLayers;   //!
   TBranch        *b_IsoTrack_pixelLayers;   //!
   TBranch        *b_IsoTrack_trackerHits;   //!
   TBranch        *b_IsoTrack_pixelHits;   //!
   TBranch        *b_IsoTrack_missingInnerPixelHits;   //!
   TBranch        *b_IsoTrack_missingOuterPixelHits;   //!
   TBranch        *b_IsoTrack_missingInnerStripHits;   //!
   TBranch        *b_IsoTrack_missingOuterStripHits;   //!
   TBranch        *b_IsoTrack_missingInnerTrackerHits;   //!
   TBranch        *b_IsoTrack_missingOuterTrackerHits;   //!
   TBranch        *b_IsoTrack_missingMiddleTrackerHits;   //!
   TBranch        *b_IsoTrack_dedxByHit0;   //!
   TBranch        *b_IsoTrack_dedxByHit1;   //!
   TBranch        *b_IsoTrack_dedxByHit2;   //!
   TBranch        *b_IsoTrack_dedxByHit3;   //!
   TBranch        *b_IsoTrack_dedxByHit4;   //!
   TBranch        *b_IsoTrack_dedxByHit5;   //!
   TBranch        *b_IsoTrack_dedxByHit6;   //!
   TBranch        *b_IsoTrack_dedxByHit7;   //!
   TBranch        *b_IsoTrack_dedxByHit8;   //!
   TBranch        *b_IsoTrack_dedxByHit9;   //!
   TBranch        *b_IsoTrack_dedxByHit10;   //!
   TBranch        *b_IsoTrack_dedxByHit11;   //!
   TBranch        *b_IsoTrack_dedxByHit12;   //!
   TBranch        *b_IsoTrack_dedxByHit13;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit0;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit1;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit2;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit3;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit4;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit5;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit6;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit7;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit8;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit9;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit10;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit11;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit12;   //!
   TBranch        *b_IsoTrack_dedxUnSmearedByHit13;   //!
   TBranch        *b_IsoTrack_deByHit0;   //!
   TBranch        *b_IsoTrack_deByHit1;   //!
   TBranch        *b_IsoTrack_deByHit2;   //!
   TBranch        *b_IsoTrack_deByHit3;   //!
   TBranch        *b_IsoTrack_deByHit4;   //!
   TBranch        *b_IsoTrack_deByHit5;   //!
   TBranch        *b_IsoTrack_deByHit6;   //!
   TBranch        *b_IsoTrack_deByHit7;   //!
   TBranch        *b_IsoTrack_deByHit8;   //!
   TBranch        *b_IsoTrack_deByHit9;   //!
   TBranch        *b_IsoTrack_deByHit10;   //!
   TBranch        *b_IsoTrack_deByHit11;   //!
   TBranch        *b_IsoTrack_deByHit12;   //!
   TBranch        *b_IsoTrack_deByHit13;   //!
   TBranch        *b_IsoTrack_dxByHit0;   //!
   TBranch        *b_IsoTrack_dxByHit1;   //!
   TBranch        *b_IsoTrack_dxByHit2;   //!
   TBranch        *b_IsoTrack_dxByHit3;   //!
   TBranch        *b_IsoTrack_dxByHit4;   //!
   TBranch        *b_IsoTrack_dxByHit5;   //!
   TBranch        *b_IsoTrack_dxByHit6;   //!
   TBranch        *b_IsoTrack_dxByHit7;   //!
   TBranch        *b_IsoTrack_dxByHit8;   //!
   TBranch        *b_IsoTrack_dxByHit9;   //!
   TBranch        *b_IsoTrack_dxByHit10;   //!
   TBranch        *b_IsoTrack_dxByHit11;   //!
   TBranch        *b_IsoTrack_dxByHit12;   //!
   TBranch        *b_IsoTrack_dxByHit13;   //!
   TBranch        *b_IsoTrack_subDetIdByHit0;   //!
   TBranch        *b_IsoTrack_subDetIdByHit1;   //!
   TBranch        *b_IsoTrack_subDetIdByHit2;   //!
   TBranch        *b_IsoTrack_subDetIdByHit3;   //!
   TBranch        *b_IsoTrack_subDetIdByHit4;   //!
   TBranch        *b_IsoTrack_subDetIdByHit5;   //!
   TBranch        *b_IsoTrack_subDetIdByHit6;   //!
   TBranch        *b_IsoTrack_subDetIdByHit7;   //!
   TBranch        *b_IsoTrack_subDetIdByHit8;   //!
   TBranch        *b_IsoTrack_subDetIdByHit9;   //!
   TBranch        *b_IsoTrack_subDetIdByHit10;   //!
   TBranch        *b_IsoTrack_subDetIdByHit11;   //!
   TBranch        *b_IsoTrack_subDetIdByHit12;   //!
   TBranch        *b_IsoTrack_subDetIdByHit13;   //!
   TBranch        *b_IsoTrack_sizeXbyHit0;   //!
   TBranch        *b_IsoTrack_sizeXbyHit1;   //!
   TBranch        *b_IsoTrack_sizeXbyHit2;   //!
   TBranch        *b_IsoTrack_sizeXbyHit3;   //!
   TBranch        *b_IsoTrack_sizeXbyHit4;   //!
   TBranch        *b_IsoTrack_sizeXbyHit5;   //!
   TBranch        *b_IsoTrack_sizeXbyHit6;   //!
   TBranch        *b_IsoTrack_sizeXbyHit7;   //!
   TBranch        *b_IsoTrack_sizeXbyHit8;   //!
   TBranch        *b_IsoTrack_sizeXbyHit9;   //!
   TBranch        *b_IsoTrack_sizeXbyHit10;   //!
   TBranch        *b_IsoTrack_sizeXbyHit11;   //!
   TBranch        *b_IsoTrack_sizeXbyHit12;   //!
   TBranch        *b_IsoTrack_sizeXbyHit13;   //!
   TBranch        *b_IsoTrack_sizeYbyHit0;   //!
   TBranch        *b_IsoTrack_sizeYbyHit1;   //!
   TBranch        *b_IsoTrack_sizeYbyHit2;   //!
   TBranch        *b_IsoTrack_sizeYbyHit3;   //!
   TBranch        *b_IsoTrack_sizeYbyHit4;   //!
   TBranch        *b_IsoTrack_sizeYbyHit5;   //!
   TBranch        *b_IsoTrack_sizeYbyHit6;   //!
   TBranch        *b_IsoTrack_sizeYbyHit7;   //!
   TBranch        *b_IsoTrack_sizeYbyHit8;   //!
   TBranch        *b_IsoTrack_sizeYbyHit9;   //!
   TBranch        *b_IsoTrack_sizeYbyHit10;   //!
   TBranch        *b_IsoTrack_sizeYbyHit11;   //!
   TBranch        *b_IsoTrack_sizeYbyHit12;   //!
   TBranch        *b_IsoTrack_sizeYbyHit13;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit0;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit1;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit2;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit3;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit4;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit5;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit6;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit7;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit8;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit9;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit10;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit11;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit12;   //!
   TBranch        *b_IsoTrack_layerOrSideByHit13;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit0;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit1;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit2;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit3;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit4;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit5;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit6;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit7;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit8;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit9;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit10;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit11;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit12;   //!
   TBranch        *b_IsoTrack_ladderOrBladeByHit13;   //!
   TBranch        *b_IsoTrack_sideByHit0;   //!
   TBranch        *b_IsoTrack_sideByHit1;   //!
   TBranch        *b_IsoTrack_sideByHit2;   //!
   TBranch        *b_IsoTrack_sideByHit3;   //!
   TBranch        *b_IsoTrack_sideByHit4;   //!
   TBranch        *b_IsoTrack_sideByHit5;   //!
   TBranch        *b_IsoTrack_sideByHit6;   //!
   TBranch        *b_IsoTrack_sideByHit7;   //!
   TBranch        *b_IsoTrack_sideByHit8;   //!
   TBranch        *b_IsoTrack_sideByHit9;   //!
   TBranch        *b_IsoTrack_sideByHit10;   //!
   TBranch        *b_IsoTrack_sideByHit11;   //!
   TBranch        *b_IsoTrack_sideByHit12;   //!
   TBranch        *b_IsoTrack_sideByHit13;   //!
   TBranch        *b_IsoTrack_pixByHit0;   //!
   TBranch        *b_IsoTrack_pixByHit1;   //!
   TBranch        *b_IsoTrack_pixByHit2;   //!
   TBranch        *b_IsoTrack_pixByHit3;   //!
   TBranch        *b_IsoTrack_pixByHit4;   //!
   TBranch        *b_IsoTrack_pixByHit5;   //!
   TBranch        *b_IsoTrack_pixByHit6;   //!
   TBranch        *b_IsoTrack_pixByHit7;   //!
   TBranch        *b_IsoTrack_pixByHit8;   //!
   TBranch        *b_IsoTrack_pixByHit9;   //!
   TBranch        *b_IsoTrack_pixByHit10;   //!
   TBranch        *b_IsoTrack_pixByHit11;   //!
   TBranch        *b_IsoTrack_pixByHit12;   //!
   TBranch        *b_IsoTrack_pixByHit13;   //!
   TBranch        *b_IsoTrack_charge;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_edxy;   //!
   TBranch        *b_IsoTrack_edz;   //!
   TBranch        *b_IsoTrack_highPurity;   //!
   TBranch        *b_IsoTrack_miniIsoCH;   //!
   TBranch        *b_IsoTrack_miniIsoNH;   //!
   TBranch        *b_IsoTrack_miniIsoPH;   //!
   TBranch        *b_IsoTrack_miniIsoPU;   //!
   TBranch        *b_IsoTrack_miniRelIso;   //!
   TBranch        *b_IsoTrack_dR03IsoCH;   //!
   TBranch        *b_IsoTrack_dR03IsoNH;   //!
   TBranch        *b_IsoTrack_dR03IsoPH;   //!
   TBranch        *b_IsoTrack_dR03IsoPU;   //!
   TBranch        *b_IsoTrack_relIso03;   //!
   TBranch        *b_IsoTrack_caloEmEnergy;   //!
   TBranch        *b_IsoTrack_caloHadEnergy;   //!
   TBranch        *b_IsoTrack_channelsGoodECAL;   //!
   TBranch        *b_IsoTrack_channelsGoodHCAL;   //!
   TBranch        *b_IsoTrack_awayJet_idx;   //!
   TBranch        *b_IsoTrack_awayJet_pt;   //!
   TBranch        *b_IsoTrack_awayJet_dr;   //!
   TBranch        *b_IsoTrack_awayNJet;   //!
   TBranch        *b_IsoTrack_awayHTJet;   //!
   TBranch        *b_IsoTrack_awayMu_idx;   //!
   TBranch        *b_IsoTrack_awayMu_dr;   //!
   TBranch        *b_IsoTrack_awayMu_mll;   //!
   TBranch        *b_IsoTrack_awayEle_idx;   //!
   TBranch        *b_IsoTrack_awayEle_dr;   //!
   TBranch        *b_IsoTrack_awayEle_mll;   //!
   TBranch        *b_IsoTrack_closestMu_idx;   //!
   TBranch        *b_IsoTrack_closestEle_idx;   //!
   TBranch        *b_IsoTrack_closestTau_idx;   //!
   TBranch        *b_IsoTrack_trigLepton_idx;   //!
   TBranch        *b_IsoTrack_myDeDx;   //!
   TBranch        *b_IsoTrack_mcMatch;   //!
   TBranch        *b_IsoTrack_mcMatchAnyId;   //!
   TBranch        *b_IsoTrack_mcMatchAnyPt;   //!
   TBranch        *b_nJetFwd;   //!
   TBranch        *b_JetFwd_id;   //!
   TBranch        *b_JetFwd_puId;   //!
   TBranch        *b_JetFwd_btagCSV;   //!
   TBranch        *b_JetFwd_btagDeepCSV;   //!
   TBranch        *b_JetFwd_rawPt;   //!
   TBranch        *b_JetFwd_mcPt;   //!
   TBranch        *b_JetFwd_mcFlavour;   //!
   TBranch        *b_JetFwd_hadronFlavour;   //!
   TBranch        *b_JetFwd_mcMatchId;   //!
   TBranch        *b_JetFwd_corr_JECUp;   //!
   TBranch        *b_JetFwd_corr_JECDown;   //!
   TBranch        *b_JetFwd_corr;   //!
   TBranch        *b_JetFwd_corr_JERUp;   //!
   TBranch        *b_JetFwd_corr_JERDown;   //!
   TBranch        *b_JetFwd_corr_JER;   //!
   TBranch        *b_JetFwd_pt;   //!
   TBranch        *b_JetFwd_eta;   //!
   TBranch        *b_JetFwd_phi;   //!
   TBranch        *b_JetFwd_mass;   //!
   TBranch        *b_JetFwd_chHEF;   //!
   TBranch        *b_JetFwd_neHEF;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_ptd;   //!
   TBranch        *b_Jet_axis2;   //!
   TBranch        *b_Jet_axis1;   //!
   TBranch        *b_Jet_mult;   //!
   TBranch        *b_Jet_partonId;   //!
   TBranch        *b_Jet_partonMotherId;   //!
   TBranch        *b_Jet_nLeptons;   //!
   TBranch        *b_Jet_id;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_Jet_btagCSV;   //!
   TBranch        *b_Jet_btagDeepCSV;   //!
   TBranch        *b_Jet_rawPt;   //!
   TBranch        *b_Jet_mcPt;   //!
   TBranch        *b_Jet_mcFlavour;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_mcMatchId;   //!
   TBranch        *b_Jet_corr_JECUp;   //!
   TBranch        *b_Jet_corr_JECDown;   //!
   TBranch        *b_Jet_corr;   //!
   TBranch        *b_Jet_corr_JERUp;   //!
   TBranch        *b_Jet_corr_JERDown;   //!
   TBranch        *b_Jet_corr_JER;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_nTauGood;   //!
   TBranch        *b_TauGood_charge;   //!
   TBranch        *b_TauGood_decayMode;   //!
   TBranch        *b_TauGood_idDecayMode;   //!
   TBranch        *b_TauGood_idDecayModeNewDMs;   //!
   TBranch        *b_TauGood_dxy;   //!
   TBranch        *b_TauGood_dz;   //!
   TBranch        *b_TauGood_idMVA;   //!
   TBranch        *b_TauGood_idMVANewDM;   //!
   TBranch        *b_TauGood_idCI3hit;   //!
   TBranch        *b_TauGood_idAntiMu;   //!
   TBranch        *b_TauGood_idAntiE;   //!
   TBranch        *b_TauGood_isoCI3hit;   //!
   TBranch        *b_TauGood_mcMatchId;   //!
   TBranch        *b_TauGood_pdgId;   //!
   TBranch        *b_TauGood_pt;   //!
   TBranch        *b_TauGood_eta;   //!
   TBranch        *b_TauGood_phi;   //!
   TBranch        *b_TauGood_mass;   //!
   TBranch        *b_TauGood_idMVAdR03;   //!
   TBranch        *b_nLepGood;   //!
   TBranch        *b_LepGood_charge;   //!
   TBranch        *b_LepGood_tightId;   //!
   TBranch        *b_LepGood_hltId;   //!
   TBranch        *b_LepGood_eleCutIdSpring15_25ns_v1;   //!
   TBranch        *b_LepGood_hasGainSwitchFlag;   //!
   TBranch        *b_LepGood_dxy;   //!
   TBranch        *b_LepGood_dz;   //!
   TBranch        *b_LepGood_edxy;   //!
   TBranch        *b_LepGood_edz;   //!
   TBranch        *b_LepGood_ip3d;   //!
   TBranch        *b_LepGood_sip3d;   //!
   TBranch        *b_LepGood_convVeto;   //!
   TBranch        *b_LepGood_lostHits;   //!
   TBranch        *b_LepGood_relIso03;   //!
   TBranch        *b_LepGood_relIso04;   //!
   TBranch        *b_LepGood_miniRelIso;   //!
   TBranch        *b_LepGood_relIsoAn04;   //!
   TBranch        *b_LepGood_tightCharge;   //!
   TBranch        *b_LepGood_mcMatchId;   //!
   TBranch        *b_LepGood_mcMatchAny;   //!
   TBranch        *b_LepGood_mcMatchTau;   //!
   TBranch        *b_LepGood_mcPt;   //!
   TBranch        *b_LepGood_mediumMuonId;   //!
   TBranch        *b_LepGood_ICHEPsoftMuonId;   //!
   TBranch        *b_LepGood_ICHEPmediumMuonId;   //!
   TBranch        *b_LepGood_pdgId;   //!
   TBranch        *b_LepGood_pt;   //!
   TBranch        *b_LepGood_eta;   //!
   TBranch        *b_LepGood_phi;   //!
   TBranch        *b_LepGood_mass;   //!
   TBranch        *b_LepGood_etaSc;   //!
   TBranch        *b_LepGood_mcMatchPdgId;   //!
   TBranch        *b_LepGood_mcPromptGamma;   //!
   TBranch        *b_LepGood_triggerMatched;   //!
   TBranch        *b_wgtsum;   //!
   TBranch        *b_kfact;   //!

   XTrackNtuple(TTree *tree=0);
   virtual ~XTrackNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef XTrackNtuple_cxx
XTrackNtuple::XTrackNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("/eos/cms/store/group/phys_exotica/xtracks/6Mar2019-Hadded/QCD_HT500to700/treeProducerXtracks/tree.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

XTrackNtuple::~XTrackNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t XTrackNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t XTrackNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void XTrackNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("nJet30", &nJet30, &b_nJet30);
   fChain->SetBranchAddress("nJet30a", &nJet30a, &b_nJet30a);
   fChain->SetBranchAddress("lheHT", &lheHT, &b_lheHT);
   fChain->SetBranchAddress("lheHTIncoming", &lheHTIncoming, &b_lheHTIncoming);
   fChain->SetBranchAddress("lheVpt", &lheVpt, &b_lheVpt);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60, &b_HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_MET", &HLT_MET, &b_HLT_MET);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v", &HLT_BIT_HLT_IsoMu24_v, &b_HLT_BIT_HLT_IsoMu24_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_v", &HLT_BIT_HLT_IsoMu24_eta2p1_v, &b_HLT_BIT_HLT_IsoMu24_eta2p1_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v", &HLT_BIT_HLT_IsoMu27_v, &b_HLT_BIT_HLT_IsoMu27_v);
   fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu, &b_HLT_SingleMu);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_WPTight_Gsf_v", &HLT_BIT_HLT_Ele32_WPTight_Gsf_v, &b_HLT_BIT_HLT_Ele32_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele35_WPTight_Gsf_v", &HLT_BIT_HLT_Ele35_WPTight_Gsf_v, &b_HLT_BIT_HLT_Ele35_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v", &HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v, &b_HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
   fChain->SetBranchAddress("HLT_SingleEl", &HLT_SingleEl, &b_HLT_SingleEl);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("met_sumEt", &met_sumEt, &b_met_sumEt);
   fChain->SetBranchAddress("met_rawPt", &met_rawPt, &b_met_rawPt);
   fChain->SetBranchAddress("met_rawPhi", &met_rawPhi, &b_met_rawPhi);
   fChain->SetBranchAddress("met_rawSumEt", &met_rawSumEt, &b_met_rawSumEt);
   fChain->SetBranchAddress("met_genPt", &met_genPt, &b_met_genPt);
   fChain->SetBranchAddress("met_genPhi", &met_genPhi, &b_met_genPhi);
   fChain->SetBranchAddress("met_genEta", &met_genEta, &b_met_genEta);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_mass", &met_mass, &b_met_mass);
   fChain->SetBranchAddress("metNoMu_pt", &metNoMu_pt, &b_metNoMu_pt);
   fChain->SetBranchAddress("metNoMu_eta", &metNoMu_eta, &b_metNoMu_eta);
   fChain->SetBranchAddress("metNoMu_phi", &metNoMu_phi, &b_metNoMu_phi);
   fChain->SetBranchAddress("metNoMu_mass", &metNoMu_mass, &b_metNoMu_mass);
   fChain->SetBranchAddress("met_jecUp_sumEt", &met_jecUp_sumEt, &b_met_jecUp_sumEt);
   fChain->SetBranchAddress("met_jecUp_rawPt", &met_jecUp_rawPt, &b_met_jecUp_rawPt);
   fChain->SetBranchAddress("met_jecUp_rawPhi", &met_jecUp_rawPhi, &b_met_jecUp_rawPhi);
   fChain->SetBranchAddress("met_jecUp_rawSumEt", &met_jecUp_rawSumEt, &b_met_jecUp_rawSumEt);
   fChain->SetBranchAddress("met_jecUp_genPt", &met_jecUp_genPt, &b_met_jecUp_genPt);
   fChain->SetBranchAddress("met_jecUp_genPhi", &met_jecUp_genPhi, &b_met_jecUp_genPhi);
   fChain->SetBranchAddress("met_jecUp_genEta", &met_jecUp_genEta, &b_met_jecUp_genEta);
   fChain->SetBranchAddress("met_jecUp_pt", &met_jecUp_pt, &b_met_jecUp_pt);
   fChain->SetBranchAddress("met_jecUp_eta", &met_jecUp_eta, &b_met_jecUp_eta);
   fChain->SetBranchAddress("met_jecUp_phi", &met_jecUp_phi, &b_met_jecUp_phi);
   fChain->SetBranchAddress("met_jecUp_mass", &met_jecUp_mass, &b_met_jecUp_mass);
   fChain->SetBranchAddress("metNoMu_jecUp_pt", &metNoMu_jecUp_pt, &b_metNoMu_jecUp_pt);
   fChain->SetBranchAddress("metNoMu_jecUp_eta", &metNoMu_jecUp_eta, &b_metNoMu_jecUp_eta);
   fChain->SetBranchAddress("metNoMu_jecUp_phi", &metNoMu_jecUp_phi, &b_metNoMu_jecUp_phi);
   fChain->SetBranchAddress("metNoMu_jecUp_mass", &metNoMu_jecUp_mass, &b_metNoMu_jecUp_mass);
//    fChain->SetBranchAddress("metNoMu_jecUp_pt", &metNoMu_jecUp_pt, &b_metNoMu_jecUp_pt);
//    fChain->SetBranchAddress("metNoMu_jecUp_eta", &metNoMu_jecUp_eta, &b_metNoMu_jecUp_eta);
//    fChain->SetBranchAddress("metNoMu_jecUp_phi", &metNoMu_jecUp_phi, &b_metNoMu_jecUp_phi);
//    fChain->SetBranchAddress("metNoMu_jecUp_mass", &metNoMu_jecUp_mass, &b_metNoMu_jecUp_mass);
   fChain->SetBranchAddress("met_jecDown_sumEt", &met_jecDown_sumEt, &b_met_jecDown_sumEt);
   fChain->SetBranchAddress("met_jecDown_rawPt", &met_jecDown_rawPt, &b_met_jecDown_rawPt);
   fChain->SetBranchAddress("met_jecDown_rawPhi", &met_jecDown_rawPhi, &b_met_jecDown_rawPhi);
   fChain->SetBranchAddress("met_jecDown_rawSumEt", &met_jecDown_rawSumEt, &b_met_jecDown_rawSumEt);
   fChain->SetBranchAddress("met_jecDown_genPt", &met_jecDown_genPt, &b_met_jecDown_genPt);
   fChain->SetBranchAddress("met_jecDown_genPhi", &met_jecDown_genPhi, &b_met_jecDown_genPhi);
   fChain->SetBranchAddress("met_jecDown_genEta", &met_jecDown_genEta, &b_met_jecDown_genEta);
   fChain->SetBranchAddress("met_jecDown_pt", &met_jecDown_pt, &b_met_jecDown_pt);
   fChain->SetBranchAddress("met_jecDown_eta", &met_jecDown_eta, &b_met_jecDown_eta);
   fChain->SetBranchAddress("met_jecDown_phi", &met_jecDown_phi, &b_met_jecDown_phi);
   fChain->SetBranchAddress("met_jecDown_mass", &met_jecDown_mass, &b_met_jecDown_mass);
   fChain->SetBranchAddress("nGenChargino", &nGenChargino, &b_nGenChargino);
   fChain->SetBranchAddress("GenChargino_motherId", &GenChargino_motherId, &b_GenChargino_motherId);
   fChain->SetBranchAddress("GenChargino_grandmotherId", &GenChargino_grandmotherId, &b_GenChargino_grandmotherId);
   fChain->SetBranchAddress("GenChargino_charge", &GenChargino_charge, &b_GenChargino_charge);
   fChain->SetBranchAddress("GenChargino_status", &GenChargino_status, &b_GenChargino_status);
   fChain->SetBranchAddress("GenChargino_isPromptHard", &GenChargino_isPromptHard, &b_GenChargino_isPromptHard);
   fChain->SetBranchAddress("GenChargino_pdgId", &GenChargino_pdgId, &b_GenChargino_pdgId);
   fChain->SetBranchAddress("GenChargino_pt", &GenChargino_pt, &b_GenChargino_pt);
   fChain->SetBranchAddress("GenChargino_eta", &GenChargino_eta, &b_GenChargino_eta);
   fChain->SetBranchAddress("GenChargino_phi", &GenChargino_phi, &b_GenChargino_phi);
   fChain->SetBranchAddress("GenChargino_mass", &GenChargino_mass, &b_GenChargino_mass);
   fChain->SetBranchAddress("GenChargino_beta", &GenChargino_beta, &b_GenChargino_beta);
   fChain->SetBranchAddress("GenChargino_decayR", &GenChargino_decayR, &b_GenChargino_decayR);
   fChain->SetBranchAddress("GenChargino_decayZ", &GenChargino_decayZ, &b_GenChargino_decayZ);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_mass", IsoTrack_mass, &b_IsoTrack_mass);
   fChain->SetBranchAddress("IsoTrack_trackerLayers", IsoTrack_trackerLayers, &b_IsoTrack_trackerLayers);
   fChain->SetBranchAddress("IsoTrack_pixelLayers", IsoTrack_pixelLayers, &b_IsoTrack_pixelLayers);
   fChain->SetBranchAddress("IsoTrack_trackerHits", IsoTrack_trackerHits, &b_IsoTrack_trackerHits);
   fChain->SetBranchAddress("IsoTrack_pixelHits", IsoTrack_pixelHits, &b_IsoTrack_pixelHits);
   fChain->SetBranchAddress("IsoTrack_missingInnerPixelHits", IsoTrack_missingInnerPixelHits, &b_IsoTrack_missingInnerPixelHits);
   fChain->SetBranchAddress("IsoTrack_missingOuterPixelHits", IsoTrack_missingOuterPixelHits, &b_IsoTrack_missingOuterPixelHits);
   fChain->SetBranchAddress("IsoTrack_missingInnerStripHits", IsoTrack_missingInnerStripHits, &b_IsoTrack_missingInnerStripHits);
   fChain->SetBranchAddress("IsoTrack_missingOuterStripHits", IsoTrack_missingOuterStripHits, &b_IsoTrack_missingOuterStripHits);
   fChain->SetBranchAddress("IsoTrack_missingInnerTrackerHits", IsoTrack_missingInnerTrackerHits, &b_IsoTrack_missingInnerTrackerHits);
   fChain->SetBranchAddress("IsoTrack_missingOuterTrackerHits", IsoTrack_missingOuterTrackerHits, &b_IsoTrack_missingOuterTrackerHits);
   fChain->SetBranchAddress("IsoTrack_missingMiddleTrackerHits", IsoTrack_missingMiddleTrackerHits, &b_IsoTrack_missingMiddleTrackerHits);
   fChain->SetBranchAddress("IsoTrack_dedxByHit0", IsoTrack_dedxByHit0, &b_IsoTrack_dedxByHit0);
   fChain->SetBranchAddress("IsoTrack_dedxByHit1", IsoTrack_dedxByHit1, &b_IsoTrack_dedxByHit1);
   fChain->SetBranchAddress("IsoTrack_dedxByHit2", IsoTrack_dedxByHit2, &b_IsoTrack_dedxByHit2);
   fChain->SetBranchAddress("IsoTrack_dedxByHit3", IsoTrack_dedxByHit3, &b_IsoTrack_dedxByHit3);
   fChain->SetBranchAddress("IsoTrack_dedxByHit4", IsoTrack_dedxByHit4, &b_IsoTrack_dedxByHit4);
   fChain->SetBranchAddress("IsoTrack_dedxByHit5", IsoTrack_dedxByHit5, &b_IsoTrack_dedxByHit5);
   fChain->SetBranchAddress("IsoTrack_dedxByHit6", IsoTrack_dedxByHit6, &b_IsoTrack_dedxByHit6);
   fChain->SetBranchAddress("IsoTrack_dedxByHit7", IsoTrack_dedxByHit7, &b_IsoTrack_dedxByHit7);
   fChain->SetBranchAddress("IsoTrack_dedxByHit8", IsoTrack_dedxByHit8, &b_IsoTrack_dedxByHit8);
   fChain->SetBranchAddress("IsoTrack_dedxByHit9", IsoTrack_dedxByHit9, &b_IsoTrack_dedxByHit9);
   fChain->SetBranchAddress("IsoTrack_dedxByHit10", IsoTrack_dedxByHit10, &b_IsoTrack_dedxByHit10);
   fChain->SetBranchAddress("IsoTrack_dedxByHit11", IsoTrack_dedxByHit11, &b_IsoTrack_dedxByHit11);
   fChain->SetBranchAddress("IsoTrack_dedxByHit12", IsoTrack_dedxByHit12, &b_IsoTrack_dedxByHit12);
   fChain->SetBranchAddress("IsoTrack_dedxByHit13", IsoTrack_dedxByHit13, &b_IsoTrack_dedxByHit13);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit0", IsoTrack_dedxUnSmearedByHit0, &b_IsoTrack_dedxUnSmearedByHit0);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit1", IsoTrack_dedxUnSmearedByHit1, &b_IsoTrack_dedxUnSmearedByHit1);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit2", IsoTrack_dedxUnSmearedByHit2, &b_IsoTrack_dedxUnSmearedByHit2);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit3", IsoTrack_dedxUnSmearedByHit3, &b_IsoTrack_dedxUnSmearedByHit3);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit4", IsoTrack_dedxUnSmearedByHit4, &b_IsoTrack_dedxUnSmearedByHit4);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit5", IsoTrack_dedxUnSmearedByHit5, &b_IsoTrack_dedxUnSmearedByHit5);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit6", IsoTrack_dedxUnSmearedByHit6, &b_IsoTrack_dedxUnSmearedByHit6);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit7", IsoTrack_dedxUnSmearedByHit7, &b_IsoTrack_dedxUnSmearedByHit7);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit8", IsoTrack_dedxUnSmearedByHit8, &b_IsoTrack_dedxUnSmearedByHit8);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit9", IsoTrack_dedxUnSmearedByHit9, &b_IsoTrack_dedxUnSmearedByHit9);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit10", IsoTrack_dedxUnSmearedByHit10, &b_IsoTrack_dedxUnSmearedByHit10);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit11", IsoTrack_dedxUnSmearedByHit11, &b_IsoTrack_dedxUnSmearedByHit11);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit12", IsoTrack_dedxUnSmearedByHit12, &b_IsoTrack_dedxUnSmearedByHit12);
   fChain->SetBranchAddress("IsoTrack_dedxUnSmearedByHit13", IsoTrack_dedxUnSmearedByHit13, &b_IsoTrack_dedxUnSmearedByHit13);
   fChain->SetBranchAddress("IsoTrack_deByHit0", IsoTrack_deByHit0, &b_IsoTrack_deByHit0);
   fChain->SetBranchAddress("IsoTrack_deByHit1", IsoTrack_deByHit1, &b_IsoTrack_deByHit1);
   fChain->SetBranchAddress("IsoTrack_deByHit2", IsoTrack_deByHit2, &b_IsoTrack_deByHit2);
   fChain->SetBranchAddress("IsoTrack_deByHit3", IsoTrack_deByHit3, &b_IsoTrack_deByHit3);
   fChain->SetBranchAddress("IsoTrack_deByHit4", IsoTrack_deByHit4, &b_IsoTrack_deByHit4);
   fChain->SetBranchAddress("IsoTrack_deByHit5", IsoTrack_deByHit5, &b_IsoTrack_deByHit5);
   fChain->SetBranchAddress("IsoTrack_deByHit6", IsoTrack_deByHit6, &b_IsoTrack_deByHit6);
   fChain->SetBranchAddress("IsoTrack_deByHit7", IsoTrack_deByHit7, &b_IsoTrack_deByHit7);
   fChain->SetBranchAddress("IsoTrack_deByHit8", IsoTrack_deByHit8, &b_IsoTrack_deByHit8);
   fChain->SetBranchAddress("IsoTrack_deByHit9", IsoTrack_deByHit9, &b_IsoTrack_deByHit9);
   fChain->SetBranchAddress("IsoTrack_deByHit10", IsoTrack_deByHit10, &b_IsoTrack_deByHit10);
   fChain->SetBranchAddress("IsoTrack_deByHit11", IsoTrack_deByHit11, &b_IsoTrack_deByHit11);
   fChain->SetBranchAddress("IsoTrack_deByHit12", IsoTrack_deByHit12, &b_IsoTrack_deByHit12);
   fChain->SetBranchAddress("IsoTrack_deByHit13", IsoTrack_deByHit13, &b_IsoTrack_deByHit13);
   fChain->SetBranchAddress("IsoTrack_dxByHit0", IsoTrack_dxByHit0, &b_IsoTrack_dxByHit0);
   fChain->SetBranchAddress("IsoTrack_dxByHit1", IsoTrack_dxByHit1, &b_IsoTrack_dxByHit1);
   fChain->SetBranchAddress("IsoTrack_dxByHit2", IsoTrack_dxByHit2, &b_IsoTrack_dxByHit2);
   fChain->SetBranchAddress("IsoTrack_dxByHit3", IsoTrack_dxByHit3, &b_IsoTrack_dxByHit3);
   fChain->SetBranchAddress("IsoTrack_dxByHit4", IsoTrack_dxByHit4, &b_IsoTrack_dxByHit4);
   fChain->SetBranchAddress("IsoTrack_dxByHit5", IsoTrack_dxByHit5, &b_IsoTrack_dxByHit5);
   fChain->SetBranchAddress("IsoTrack_dxByHit6", IsoTrack_dxByHit6, &b_IsoTrack_dxByHit6);
   fChain->SetBranchAddress("IsoTrack_dxByHit7", IsoTrack_dxByHit7, &b_IsoTrack_dxByHit7);
   fChain->SetBranchAddress("IsoTrack_dxByHit8", IsoTrack_dxByHit8, &b_IsoTrack_dxByHit8);
   fChain->SetBranchAddress("IsoTrack_dxByHit9", IsoTrack_dxByHit9, &b_IsoTrack_dxByHit9);
   fChain->SetBranchAddress("IsoTrack_dxByHit10", IsoTrack_dxByHit10, &b_IsoTrack_dxByHit10);
   fChain->SetBranchAddress("IsoTrack_dxByHit11", IsoTrack_dxByHit11, &b_IsoTrack_dxByHit11);
   fChain->SetBranchAddress("IsoTrack_dxByHit12", IsoTrack_dxByHit12, &b_IsoTrack_dxByHit12);
   fChain->SetBranchAddress("IsoTrack_dxByHit13", IsoTrack_dxByHit13, &b_IsoTrack_dxByHit13);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit0", IsoTrack_subDetIdByHit0, &b_IsoTrack_subDetIdByHit0);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit1", IsoTrack_subDetIdByHit1, &b_IsoTrack_subDetIdByHit1);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit2", IsoTrack_subDetIdByHit2, &b_IsoTrack_subDetIdByHit2);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit3", IsoTrack_subDetIdByHit3, &b_IsoTrack_subDetIdByHit3);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit4", IsoTrack_subDetIdByHit4, &b_IsoTrack_subDetIdByHit4);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit5", IsoTrack_subDetIdByHit5, &b_IsoTrack_subDetIdByHit5);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit6", IsoTrack_subDetIdByHit6, &b_IsoTrack_subDetIdByHit6);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit7", IsoTrack_subDetIdByHit7, &b_IsoTrack_subDetIdByHit7);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit8", IsoTrack_subDetIdByHit8, &b_IsoTrack_subDetIdByHit8);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit9", IsoTrack_subDetIdByHit9, &b_IsoTrack_subDetIdByHit9);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit10", IsoTrack_subDetIdByHit10, &b_IsoTrack_subDetIdByHit10);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit11", IsoTrack_subDetIdByHit11, &b_IsoTrack_subDetIdByHit11);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit12", IsoTrack_subDetIdByHit12, &b_IsoTrack_subDetIdByHit12);
   fChain->SetBranchAddress("IsoTrack_subDetIdByHit13", IsoTrack_subDetIdByHit13, &b_IsoTrack_subDetIdByHit13);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit0", IsoTrack_sizeXbyHit0, &b_IsoTrack_sizeXbyHit0);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit1", IsoTrack_sizeXbyHit1, &b_IsoTrack_sizeXbyHit1);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit2", IsoTrack_sizeXbyHit2, &b_IsoTrack_sizeXbyHit2);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit3", IsoTrack_sizeXbyHit3, &b_IsoTrack_sizeXbyHit3);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit4", IsoTrack_sizeXbyHit4, &b_IsoTrack_sizeXbyHit4);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit5", IsoTrack_sizeXbyHit5, &b_IsoTrack_sizeXbyHit5);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit6", IsoTrack_sizeXbyHit6, &b_IsoTrack_sizeXbyHit6);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit7", IsoTrack_sizeXbyHit7, &b_IsoTrack_sizeXbyHit7);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit8", IsoTrack_sizeXbyHit8, &b_IsoTrack_sizeXbyHit8);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit9", IsoTrack_sizeXbyHit9, &b_IsoTrack_sizeXbyHit9);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit10", IsoTrack_sizeXbyHit10, &b_IsoTrack_sizeXbyHit10);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit11", IsoTrack_sizeXbyHit11, &b_IsoTrack_sizeXbyHit11);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit12", IsoTrack_sizeXbyHit12, &b_IsoTrack_sizeXbyHit12);
   fChain->SetBranchAddress("IsoTrack_sizeXbyHit13", IsoTrack_sizeXbyHit13, &b_IsoTrack_sizeXbyHit13);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit0", IsoTrack_sizeYbyHit0, &b_IsoTrack_sizeYbyHit0);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit1", IsoTrack_sizeYbyHit1, &b_IsoTrack_sizeYbyHit1);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit2", IsoTrack_sizeYbyHit2, &b_IsoTrack_sizeYbyHit2);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit3", IsoTrack_sizeYbyHit3, &b_IsoTrack_sizeYbyHit3);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit4", IsoTrack_sizeYbyHit4, &b_IsoTrack_sizeYbyHit4);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit5", IsoTrack_sizeYbyHit5, &b_IsoTrack_sizeYbyHit5);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit6", IsoTrack_sizeYbyHit6, &b_IsoTrack_sizeYbyHit6);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit7", IsoTrack_sizeYbyHit7, &b_IsoTrack_sizeYbyHit7);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit8", IsoTrack_sizeYbyHit8, &b_IsoTrack_sizeYbyHit8);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit9", IsoTrack_sizeYbyHit9, &b_IsoTrack_sizeYbyHit9);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit10", IsoTrack_sizeYbyHit10, &b_IsoTrack_sizeYbyHit10);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit11", IsoTrack_sizeYbyHit11, &b_IsoTrack_sizeYbyHit11);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit12", IsoTrack_sizeYbyHit12, &b_IsoTrack_sizeYbyHit12);
   fChain->SetBranchAddress("IsoTrack_sizeYbyHit13", IsoTrack_sizeYbyHit13, &b_IsoTrack_sizeYbyHit13);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit0", IsoTrack_layerOrSideByHit0, &b_IsoTrack_layerOrSideByHit0);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit1", IsoTrack_layerOrSideByHit1, &b_IsoTrack_layerOrSideByHit1);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit2", IsoTrack_layerOrSideByHit2, &b_IsoTrack_layerOrSideByHit2);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit3", IsoTrack_layerOrSideByHit3, &b_IsoTrack_layerOrSideByHit3);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit4", IsoTrack_layerOrSideByHit4, &b_IsoTrack_layerOrSideByHit4);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit5", IsoTrack_layerOrSideByHit5, &b_IsoTrack_layerOrSideByHit5);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit6", IsoTrack_layerOrSideByHit6, &b_IsoTrack_layerOrSideByHit6);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit7", IsoTrack_layerOrSideByHit7, &b_IsoTrack_layerOrSideByHit7);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit8", IsoTrack_layerOrSideByHit8, &b_IsoTrack_layerOrSideByHit8);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit9", IsoTrack_layerOrSideByHit9, &b_IsoTrack_layerOrSideByHit9);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit10", IsoTrack_layerOrSideByHit10, &b_IsoTrack_layerOrSideByHit10);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit11", IsoTrack_layerOrSideByHit11, &b_IsoTrack_layerOrSideByHit11);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit12", IsoTrack_layerOrSideByHit12, &b_IsoTrack_layerOrSideByHit12);
   fChain->SetBranchAddress("IsoTrack_layerOrSideByHit13", IsoTrack_layerOrSideByHit13, &b_IsoTrack_layerOrSideByHit13);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit0", IsoTrack_ladderOrBladeByHit0, &b_IsoTrack_ladderOrBladeByHit0);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit1", IsoTrack_ladderOrBladeByHit1, &b_IsoTrack_ladderOrBladeByHit1);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit2", IsoTrack_ladderOrBladeByHit2, &b_IsoTrack_ladderOrBladeByHit2);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit3", IsoTrack_ladderOrBladeByHit3, &b_IsoTrack_ladderOrBladeByHit3);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit4", IsoTrack_ladderOrBladeByHit4, &b_IsoTrack_ladderOrBladeByHit4);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit5", IsoTrack_ladderOrBladeByHit5, &b_IsoTrack_ladderOrBladeByHit5);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit6", IsoTrack_ladderOrBladeByHit6, &b_IsoTrack_ladderOrBladeByHit6);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit7", IsoTrack_ladderOrBladeByHit7, &b_IsoTrack_ladderOrBladeByHit7);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit8", IsoTrack_ladderOrBladeByHit8, &b_IsoTrack_ladderOrBladeByHit8);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit9", IsoTrack_ladderOrBladeByHit9, &b_IsoTrack_ladderOrBladeByHit9);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit10", IsoTrack_ladderOrBladeByHit10, &b_IsoTrack_ladderOrBladeByHit10);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit11", IsoTrack_ladderOrBladeByHit11, &b_IsoTrack_ladderOrBladeByHit11);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit12", IsoTrack_ladderOrBladeByHit12, &b_IsoTrack_ladderOrBladeByHit12);
   fChain->SetBranchAddress("IsoTrack_ladderOrBladeByHit13", IsoTrack_ladderOrBladeByHit13, &b_IsoTrack_ladderOrBladeByHit13);
   fChain->SetBranchAddress("IsoTrack_sideByHit0", IsoTrack_sideByHit0, &b_IsoTrack_sideByHit0);
   fChain->SetBranchAddress("IsoTrack_sideByHit1", IsoTrack_sideByHit1, &b_IsoTrack_sideByHit1);
   fChain->SetBranchAddress("IsoTrack_sideByHit2", IsoTrack_sideByHit2, &b_IsoTrack_sideByHit2);
   fChain->SetBranchAddress("IsoTrack_sideByHit3", IsoTrack_sideByHit3, &b_IsoTrack_sideByHit3);
   fChain->SetBranchAddress("IsoTrack_sideByHit4", IsoTrack_sideByHit4, &b_IsoTrack_sideByHit4);
   fChain->SetBranchAddress("IsoTrack_sideByHit5", IsoTrack_sideByHit5, &b_IsoTrack_sideByHit5);
   fChain->SetBranchAddress("IsoTrack_sideByHit6", IsoTrack_sideByHit6, &b_IsoTrack_sideByHit6);
   fChain->SetBranchAddress("IsoTrack_sideByHit7", IsoTrack_sideByHit7, &b_IsoTrack_sideByHit7);
   fChain->SetBranchAddress("IsoTrack_sideByHit8", IsoTrack_sideByHit8, &b_IsoTrack_sideByHit8);
   fChain->SetBranchAddress("IsoTrack_sideByHit9", IsoTrack_sideByHit9, &b_IsoTrack_sideByHit9);
   fChain->SetBranchAddress("IsoTrack_sideByHit10", IsoTrack_sideByHit10, &b_IsoTrack_sideByHit10);
   fChain->SetBranchAddress("IsoTrack_sideByHit11", IsoTrack_sideByHit11, &b_IsoTrack_sideByHit11);
   fChain->SetBranchAddress("IsoTrack_sideByHit12", IsoTrack_sideByHit12, &b_IsoTrack_sideByHit12);
   fChain->SetBranchAddress("IsoTrack_sideByHit13", IsoTrack_sideByHit13, &b_IsoTrack_sideByHit13);
   fChain->SetBranchAddress("IsoTrack_pixByHit0", IsoTrack_pixByHit0, &b_IsoTrack_pixByHit0);
   fChain->SetBranchAddress("IsoTrack_pixByHit1", IsoTrack_pixByHit1, &b_IsoTrack_pixByHit1);
   fChain->SetBranchAddress("IsoTrack_pixByHit2", IsoTrack_pixByHit2, &b_IsoTrack_pixByHit2);
   fChain->SetBranchAddress("IsoTrack_pixByHit3", IsoTrack_pixByHit3, &b_IsoTrack_pixByHit3);
   fChain->SetBranchAddress("IsoTrack_pixByHit4", IsoTrack_pixByHit4, &b_IsoTrack_pixByHit4);
   fChain->SetBranchAddress("IsoTrack_pixByHit5", IsoTrack_pixByHit5, &b_IsoTrack_pixByHit5);
   fChain->SetBranchAddress("IsoTrack_pixByHit6", IsoTrack_pixByHit6, &b_IsoTrack_pixByHit6);
   fChain->SetBranchAddress("IsoTrack_pixByHit7", IsoTrack_pixByHit7, &b_IsoTrack_pixByHit7);
   fChain->SetBranchAddress("IsoTrack_pixByHit8", IsoTrack_pixByHit8, &b_IsoTrack_pixByHit8);
   fChain->SetBranchAddress("IsoTrack_pixByHit9", IsoTrack_pixByHit9, &b_IsoTrack_pixByHit9);
   fChain->SetBranchAddress("IsoTrack_pixByHit10", IsoTrack_pixByHit10, &b_IsoTrack_pixByHit10);
   fChain->SetBranchAddress("IsoTrack_pixByHit11", IsoTrack_pixByHit11, &b_IsoTrack_pixByHit11);
   fChain->SetBranchAddress("IsoTrack_pixByHit12", IsoTrack_pixByHit12, &b_IsoTrack_pixByHit12);
   fChain->SetBranchAddress("IsoTrack_pixByHit13", IsoTrack_pixByHit13, &b_IsoTrack_pixByHit13);
   fChain->SetBranchAddress("IsoTrack_charge", IsoTrack_charge, &b_IsoTrack_charge);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_edxy", IsoTrack_edxy, &b_IsoTrack_edxy);
   fChain->SetBranchAddress("IsoTrack_edz", IsoTrack_edz, &b_IsoTrack_edz);
   fChain->SetBranchAddress("IsoTrack_highPurity", IsoTrack_highPurity, &b_IsoTrack_highPurity);
   fChain->SetBranchAddress("IsoTrack_miniIsoCH", IsoTrack_miniIsoCH, &b_IsoTrack_miniIsoCH);
   fChain->SetBranchAddress("IsoTrack_miniIsoNH", IsoTrack_miniIsoNH, &b_IsoTrack_miniIsoNH);
   fChain->SetBranchAddress("IsoTrack_miniIsoPH", IsoTrack_miniIsoPH, &b_IsoTrack_miniIsoPH);
   fChain->SetBranchAddress("IsoTrack_miniIsoPU", IsoTrack_miniIsoPU, &b_IsoTrack_miniIsoPU);
   fChain->SetBranchAddress("IsoTrack_miniRelIso", IsoTrack_miniRelIso, &b_IsoTrack_miniRelIso);
   fChain->SetBranchAddress("IsoTrack_dR03IsoCH", IsoTrack_dR03IsoCH, &b_IsoTrack_dR03IsoCH);
   fChain->SetBranchAddress("IsoTrack_dR03IsoNH", IsoTrack_dR03IsoNH, &b_IsoTrack_dR03IsoNH);
   fChain->SetBranchAddress("IsoTrack_dR03IsoPH", IsoTrack_dR03IsoPH, &b_IsoTrack_dR03IsoPH);
   fChain->SetBranchAddress("IsoTrack_dR03IsoPU", IsoTrack_dR03IsoPU, &b_IsoTrack_dR03IsoPU);
   fChain->SetBranchAddress("IsoTrack_relIso03", IsoTrack_relIso03, &b_IsoTrack_relIso03);
   fChain->SetBranchAddress("IsoTrack_caloEmEnergy", IsoTrack_caloEmEnergy, &b_IsoTrack_caloEmEnergy);
   fChain->SetBranchAddress("IsoTrack_caloHadEnergy", IsoTrack_caloHadEnergy, &b_IsoTrack_caloHadEnergy);
   fChain->SetBranchAddress("IsoTrack_channelsGoodECAL", IsoTrack_channelsGoodECAL, &b_IsoTrack_channelsGoodECAL);
   fChain->SetBranchAddress("IsoTrack_channelsGoodHCAL", IsoTrack_channelsGoodHCAL, &b_IsoTrack_channelsGoodHCAL);
   fChain->SetBranchAddress("IsoTrack_awayJet_idx", IsoTrack_awayJet_idx, &b_IsoTrack_awayJet_idx);
   fChain->SetBranchAddress("IsoTrack_awayJet_pt", IsoTrack_awayJet_pt, &b_IsoTrack_awayJet_pt);
   fChain->SetBranchAddress("IsoTrack_awayJet_dr", IsoTrack_awayJet_dr, &b_IsoTrack_awayJet_dr);
   fChain->SetBranchAddress("IsoTrack_awayNJet", IsoTrack_awayNJet, &b_IsoTrack_awayNJet);
   fChain->SetBranchAddress("IsoTrack_awayHTJet", IsoTrack_awayHTJet, &b_IsoTrack_awayHTJet);
   fChain->SetBranchAddress("IsoTrack_awayMu_idx", IsoTrack_awayMu_idx, &b_IsoTrack_awayMu_idx);
   fChain->SetBranchAddress("IsoTrack_awayMu_dr", IsoTrack_awayMu_dr, &b_IsoTrack_awayMu_dr);
   fChain->SetBranchAddress("IsoTrack_awayMu_mll", IsoTrack_awayMu_mll, &b_IsoTrack_awayMu_mll);
   fChain->SetBranchAddress("IsoTrack_awayEle_idx", IsoTrack_awayEle_idx, &b_IsoTrack_awayEle_idx);
   fChain->SetBranchAddress("IsoTrack_awayEle_dr", IsoTrack_awayEle_dr, &b_IsoTrack_awayEle_dr);
   fChain->SetBranchAddress("IsoTrack_awayEle_mll", IsoTrack_awayEle_mll, &b_IsoTrack_awayEle_mll);
   fChain->SetBranchAddress("IsoTrack_closestMu_idx", IsoTrack_closestMu_idx, &b_IsoTrack_closestMu_idx);
   fChain->SetBranchAddress("IsoTrack_closestEle_idx", IsoTrack_closestEle_idx, &b_IsoTrack_closestEle_idx);
   fChain->SetBranchAddress("IsoTrack_closestTau_idx", IsoTrack_closestTau_idx, &b_IsoTrack_closestTau_idx);
   fChain->SetBranchAddress("IsoTrack_trigLepton_idx", IsoTrack_trigLepton_idx, &b_IsoTrack_trigLepton_idx);
   fChain->SetBranchAddress("IsoTrack_myDeDx", IsoTrack_myDeDx, &b_IsoTrack_myDeDx);
   fChain->SetBranchAddress("IsoTrack_mcMatch", IsoTrack_mcMatch, &b_IsoTrack_mcMatch);
   fChain->SetBranchAddress("IsoTrack_mcMatchAnyId", IsoTrack_mcMatchAnyId, &b_IsoTrack_mcMatchAnyId);
   fChain->SetBranchAddress("IsoTrack_mcMatchAnyPt", IsoTrack_mcMatchAnyPt, &b_IsoTrack_mcMatchAnyPt);
   fChain->SetBranchAddress("nJetFwd", &nJetFwd, &b_nJetFwd);
   fChain->SetBranchAddress("JetFwd_id", JetFwd_id, &b_JetFwd_id);
   fChain->SetBranchAddress("JetFwd_puId", JetFwd_puId, &b_JetFwd_puId);
   fChain->SetBranchAddress("JetFwd_btagCSV", JetFwd_btagCSV, &b_JetFwd_btagCSV);
   fChain->SetBranchAddress("JetFwd_btagDeepCSV", JetFwd_btagDeepCSV, &b_JetFwd_btagDeepCSV);
   fChain->SetBranchAddress("JetFwd_rawPt", JetFwd_rawPt, &b_JetFwd_rawPt);
   fChain->SetBranchAddress("JetFwd_mcPt", JetFwd_mcPt, &b_JetFwd_mcPt);
   fChain->SetBranchAddress("JetFwd_mcFlavour", JetFwd_mcFlavour, &b_JetFwd_mcFlavour);
   fChain->SetBranchAddress("JetFwd_hadronFlavour", JetFwd_hadronFlavour, &b_JetFwd_hadronFlavour);
   fChain->SetBranchAddress("JetFwd_mcMatchId", JetFwd_mcMatchId, &b_JetFwd_mcMatchId);
   fChain->SetBranchAddress("JetFwd_corr_JECUp", JetFwd_corr_JECUp, &b_JetFwd_corr_JECUp);
   fChain->SetBranchAddress("JetFwd_corr_JECDown", JetFwd_corr_JECDown, &b_JetFwd_corr_JECDown);
   fChain->SetBranchAddress("JetFwd_corr", JetFwd_corr, &b_JetFwd_corr);
   fChain->SetBranchAddress("JetFwd_corr_JERUp", JetFwd_corr_JERUp, &b_JetFwd_corr_JERUp);
   fChain->SetBranchAddress("JetFwd_corr_JERDown", JetFwd_corr_JERDown, &b_JetFwd_corr_JERDown);
   fChain->SetBranchAddress("JetFwd_corr_JER", JetFwd_corr_JER, &b_JetFwd_corr_JER);
   fChain->SetBranchAddress("JetFwd_pt", JetFwd_pt, &b_JetFwd_pt);
   fChain->SetBranchAddress("JetFwd_eta", JetFwd_eta, &b_JetFwd_eta);
   fChain->SetBranchAddress("JetFwd_phi", JetFwd_phi, &b_JetFwd_phi);
   fChain->SetBranchAddress("JetFwd_mass", JetFwd_mass, &b_JetFwd_mass);
   fChain->SetBranchAddress("JetFwd_chHEF", JetFwd_chHEF, &b_JetFwd_chHEF);
   fChain->SetBranchAddress("JetFwd_neHEF", JetFwd_neHEF, &b_JetFwd_neHEF);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_ptd", Jet_ptd, &b_Jet_ptd);
   fChain->SetBranchAddress("Jet_axis2", Jet_axis2, &b_Jet_axis2);
   fChain->SetBranchAddress("Jet_axis1", Jet_axis1, &b_Jet_axis1);
   fChain->SetBranchAddress("Jet_mult", Jet_mult, &b_Jet_mult);
   fChain->SetBranchAddress("Jet_partonId", Jet_partonId, &b_Jet_partonId);
   fChain->SetBranchAddress("Jet_partonMotherId", Jet_partonMotherId, &b_Jet_partonMotherId);
   fChain->SetBranchAddress("Jet_nLeptons", Jet_nLeptons, &b_Jet_nLeptons);
   fChain->SetBranchAddress("Jet_id", Jet_id, &b_Jet_id);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV, &b_Jet_btagCSV);
   fChain->SetBranchAddress("Jet_btagDeepCSV", Jet_btagDeepCSV, &b_Jet_btagDeepCSV);
   fChain->SetBranchAddress("Jet_rawPt", Jet_rawPt, &b_Jet_rawPt);
   fChain->SetBranchAddress("Jet_mcPt", Jet_mcPt, &b_Jet_mcPt);
   fChain->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_mcMatchId", Jet_mcMatchId, &b_Jet_mcMatchId);
   fChain->SetBranchAddress("Jet_corr_JECUp", Jet_corr_JECUp, &b_Jet_corr_JECUp);
   fChain->SetBranchAddress("Jet_corr_JECDown", Jet_corr_JECDown, &b_Jet_corr_JECDown);
   fChain->SetBranchAddress("Jet_corr", Jet_corr, &b_Jet_corr);
   fChain->SetBranchAddress("Jet_corr_JERUp", Jet_corr_JERUp, &b_Jet_corr_JERUp);
   fChain->SetBranchAddress("Jet_corr_JERDown", Jet_corr_JERDown, &b_Jet_corr_JERDown);
   fChain->SetBranchAddress("Jet_corr_JER", Jet_corr_JER, &b_Jet_corr_JER);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("nTauGood", &nTauGood, &b_nTauGood);
   fChain->SetBranchAddress("TauGood_charge", TauGood_charge, &b_TauGood_charge);
   fChain->SetBranchAddress("TauGood_decayMode", TauGood_decayMode, &b_TauGood_decayMode);
   fChain->SetBranchAddress("TauGood_idDecayMode", TauGood_idDecayMode, &b_TauGood_idDecayMode);
   fChain->SetBranchAddress("TauGood_idDecayModeNewDMs", TauGood_idDecayModeNewDMs, &b_TauGood_idDecayModeNewDMs);
   fChain->SetBranchAddress("TauGood_dxy", TauGood_dxy, &b_TauGood_dxy);
   fChain->SetBranchAddress("TauGood_dz", TauGood_dz, &b_TauGood_dz);
   fChain->SetBranchAddress("TauGood_idMVA", TauGood_idMVA, &b_TauGood_idMVA);
   fChain->SetBranchAddress("TauGood_idMVANewDM", TauGood_idMVANewDM, &b_TauGood_idMVANewDM);
   fChain->SetBranchAddress("TauGood_idCI3hit", TauGood_idCI3hit, &b_TauGood_idCI3hit);
   fChain->SetBranchAddress("TauGood_idAntiMu", TauGood_idAntiMu, &b_TauGood_idAntiMu);
   fChain->SetBranchAddress("TauGood_idAntiE", TauGood_idAntiE, &b_TauGood_idAntiE);
   fChain->SetBranchAddress("TauGood_isoCI3hit", TauGood_isoCI3hit, &b_TauGood_isoCI3hit);
   fChain->SetBranchAddress("TauGood_mcMatchId", TauGood_mcMatchId, &b_TauGood_mcMatchId);
   fChain->SetBranchAddress("TauGood_pdgId", TauGood_pdgId, &b_TauGood_pdgId);
   fChain->SetBranchAddress("TauGood_pt", TauGood_pt, &b_TauGood_pt);
   fChain->SetBranchAddress("TauGood_eta", TauGood_eta, &b_TauGood_eta);
   fChain->SetBranchAddress("TauGood_phi", TauGood_phi, &b_TauGood_phi);
   fChain->SetBranchAddress("TauGood_mass", TauGood_mass, &b_TauGood_mass);
   fChain->SetBranchAddress("TauGood_idMVAdR03", TauGood_idMVAdR03, &b_TauGood_idMVAdR03);
   fChain->SetBranchAddress("nLepGood", &nLepGood, &b_nLepGood);
   fChain->SetBranchAddress("LepGood_charge", LepGood_charge, &b_LepGood_charge);
   fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId, &b_LepGood_tightId);
   fChain->SetBranchAddress("LepGood_hltId", LepGood_hltId, &b_LepGood_hltId);
   fChain->SetBranchAddress("LepGood_eleCutIdSpring15_25ns_v1", LepGood_eleCutIdSpring15_25ns_v1, &b_LepGood_eleCutIdSpring15_25ns_v1);
   fChain->SetBranchAddress("LepGood_hasGainSwitchFlag", LepGood_hasGainSwitchFlag, &b_LepGood_hasGainSwitchFlag);
   fChain->SetBranchAddress("LepGood_dxy", LepGood_dxy, &b_LepGood_dxy);
   fChain->SetBranchAddress("LepGood_dz", LepGood_dz, &b_LepGood_dz);
   fChain->SetBranchAddress("LepGood_edxy", LepGood_edxy, &b_LepGood_edxy);
   fChain->SetBranchAddress("LepGood_edz", LepGood_edz, &b_LepGood_edz);
   fChain->SetBranchAddress("LepGood_ip3d", LepGood_ip3d, &b_LepGood_ip3d);
   fChain->SetBranchAddress("LepGood_sip3d", LepGood_sip3d, &b_LepGood_sip3d);
   fChain->SetBranchAddress("LepGood_convVeto", LepGood_convVeto, &b_LepGood_convVeto);
   fChain->SetBranchAddress("LepGood_lostHits", LepGood_lostHits, &b_LepGood_lostHits);
   fChain->SetBranchAddress("LepGood_relIso03", LepGood_relIso03, &b_LepGood_relIso03);
   fChain->SetBranchAddress("LepGood_relIso04", LepGood_relIso04, &b_LepGood_relIso04);
   fChain->SetBranchAddress("LepGood_miniRelIso", LepGood_miniRelIso, &b_LepGood_miniRelIso);
   fChain->SetBranchAddress("LepGood_relIsoAn04", LepGood_relIsoAn04, &b_LepGood_relIsoAn04);
   fChain->SetBranchAddress("LepGood_tightCharge", LepGood_tightCharge, &b_LepGood_tightCharge);
   fChain->SetBranchAddress("LepGood_mcMatchId", LepGood_mcMatchId, &b_LepGood_mcMatchId);
   fChain->SetBranchAddress("LepGood_mcMatchAny", LepGood_mcMatchAny, &b_LepGood_mcMatchAny);
   fChain->SetBranchAddress("LepGood_mcMatchTau", LepGood_mcMatchTau, &b_LepGood_mcMatchTau);
   fChain->SetBranchAddress("LepGood_mcPt", LepGood_mcPt, &b_LepGood_mcPt);
   fChain->SetBranchAddress("LepGood_mediumMuonId", LepGood_mediumMuonId, &b_LepGood_mediumMuonId);
   fChain->SetBranchAddress("LepGood_ICHEPsoftMuonId", LepGood_ICHEPsoftMuonId, &b_LepGood_ICHEPsoftMuonId);
   fChain->SetBranchAddress("LepGood_ICHEPmediumMuonId", LepGood_ICHEPmediumMuonId, &b_LepGood_ICHEPmediumMuonId);
   fChain->SetBranchAddress("LepGood_pdgId", LepGood_pdgId, &b_LepGood_pdgId);
   fChain->SetBranchAddress("LepGood_pt", LepGood_pt, &b_LepGood_pt);
   fChain->SetBranchAddress("LepGood_eta", LepGood_eta, &b_LepGood_eta);
   fChain->SetBranchAddress("LepGood_phi", LepGood_phi, &b_LepGood_phi);
   fChain->SetBranchAddress("LepGood_mass", LepGood_mass, &b_LepGood_mass);
   fChain->SetBranchAddress("LepGood_etaSc", LepGood_etaSc, &b_LepGood_etaSc);
   fChain->SetBranchAddress("LepGood_mcMatchPdgId", LepGood_mcMatchPdgId, &b_LepGood_mcMatchPdgId);
   fChain->SetBranchAddress("LepGood_mcPromptGamma", LepGood_mcPromptGamma, &b_LepGood_mcPromptGamma);
   fChain->SetBranchAddress("LepGood_triggerMatched", LepGood_triggerMatched, &b_LepGood_triggerMatched);
   fChain->SetBranchAddress("wgtsum", &wgtsum, &b_wgtsum);
   fChain->SetBranchAddress("kfact", &kfact, &b_kfact);
   Notify();
}

Bool_t XTrackNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void XTrackNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t XTrackNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef XTrackNtuple_cxx
