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
#include "TFile.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TLegend.h"

// Namespaces //
using namespace std;

// Typedefs //

//   lumi  events
typedef map< UInt_t , vector<ULong64_t> > MapLumiEvent;


// Declare functions //
Int_t plotEfficiency(TString plotDir, TString version, TString era, const UInt_t nDS, 
		     TString* nameDS, TString* titleDS,Int_t* colorDS, Int_t* styleDS);

Int_t processChain(Int_t nEvents, TChain* fChain, TString nameDS, Bool_t isMC, TH1F* hDen, TH1F* hNum, 
		   Bool_t doCutOnMetFlags, Bool_t doPrintBeforeCuts);


// MAIN //
Int_t loop(TString version="v0_test"   , Bool_t doReadChains=kTRUE, 
	   TString era="2018"          , TString sample="all",
	   Int_t nEvents=-1, 
	   Bool_t doCutOnMetFlags=kTRUE, Bool_t doPrintBeforeCuts=kTRUE,
	   Bool_t doConstantBinning=kFALSE)
{

  // Data Sets //
  const UInt_t nDS=3;
  TChain* chain  [nDS];
  TString nameDS [nDS] = {"data_1mu_"+era, "mc_wjets_"+era, "mc_zjets_"+era};
  TString titleDS[nDS] = {"SingleMuon "+era, "W+Jets "+era, "DYJetsM50 "+era};
  Int_t   colorDS[nDS] = {kBlack, kBlue, kRed};
  Int_t   styleDS[nDS] = {kFullTriangleUp, kOpenCircle, kOpenSquare};
  // Int_t* colorDS, Int_t* styleDS

  Bool_t  mcDS[   nDS] = {kFALSE, kTRUE, kTRUE};

  TString baseDir, subDir, subDirData, subDirMC, suffixDir, fullPath;
  vector<TString> dirDS[ nDS]; 
  // fullPath = baseDir + subDir + dirDS[iS][iDS] + suffixDir;

  // 2018 Samples
  if(era=="2018") {

    baseDir    = "/eos/cms/store/group/phys_exotica/xtracks/";
    subDirData = "/1Apr2019-Hadded-2018/";
    subDirMC   = "";
    suffixDir  = "/treeProducerXtracks/tree.root";

    dirDS[0] = {"SingleMuon_Run2018A", "SingleMuon_Run2018B", "SingleMuon_Run2018C"};

    dirDS[1] = { "23April2019_Samples2018_Hadded/WJets_HT100to200/",
		 "23April2019_Samples2018_Hadded/WJets_HT200to400/",
		 "23April2019_Samples2018_Hadded/WJets_HT400to600/",
		 "23April2019_Samples2018_Hadded/WJets_HT600to800/",
		 "23April2019_Samples2018_Hadded/WJets_HT800to1200/",
		 "23April2019_Samples2018_Hadded/WJets_HT1200to2500/",
		 "23April2019_Samples2018_Hadded/WJets_HT2500toInf/" };

    dirDS[2] = { "23April2019_Samples2018_Hadded/DYJetsM50_HT100to200/",
		 "23April2019_Samples2018_Hadded/DYJetsM50_HT200to400/",
		 "1Apr2019-Hadded-2018/DYJetsM50_HT400to600/", 
		 "23April2019_Samples2018_Hadded/DYJetsM50_HT600to800/",
		 "1Apr2019-Hadded-2018/DYJetsM50_HT800to1200/" };
  }

  // 2017 Samples
  else if(era=="2017") {

    baseDir    = "/eos/cms/store/"; 
    subDirData = "/cmst3/user/amassiro/CMG/Data-CR1L-NewGeometry-Calibrated-and-Smeared-Correct/";
    subDirMC   = "";
    suffixDir  = "";

    dirDS[0] = {"tree_SingleMuon_*.root"};

    dirDS[1] = { "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT100to200/treeProducerXtracks/tree.root",
		 "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT1200to2500/treeProducerXtracks/tree.root",
		 "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT200to400/treeProducerXtracks/tree.root",
		 "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT400to600/treeProducerXtracks/tree.root",
		 "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT600to800/treeProducerXtracks/tree.root",
		 "/group/phys_exotica/xtracks/23April2019_Samples2017_Hadded/WJets_HT800to1200/treeProducerXtracks/tree.root" };

    dirDS[2] = {"/cmst3/user/amassiro/CMG/MC-CR1L-NewGeometry-Calibrated-and-Smeared-Correct/tree_DY*.root"};

  }

  else {
    cout << "ERROR: choose era between 2018 and 2017. Exit." << endl;
    return -1;
  }

  // Output results //
  TString plotDir = "plots/trigger/";
  TString fmode = "read";
  if(doReadChains) fmode = "recreate";
  TFile *fout = 0;

  // Define Histograms //
  TH1F*   hDen[ nDS];
  TH1F*   hNum[ nDS];
  //TEfficiency* teff[nDS];

  const UInt_t nBinsMet = 26;
  Float_t bins_met[nBinsMet]  = { 050,  075,  100,  110,  120, 
				  130,  140,  150,  160,  170, 
				  180,  190,  200,  220,  245, 
				  255,  275,  300,  350,  400,  
				  500,  650,  800, 1000, 1300, 
				  2000};

  if(doReadChains) {
    for(UInt_t iS=0; iS<nDS; iS++) {

      if(nameDS[iS]!=sample && sample!="all") continue;

      // Define output root file to store histograms
      fout = new TFile(plotDir + "/" + version + "/histos_"+era+"_"+sample+".root", fmode);

      // Define histograms
      if(doConstantBinning) {
	hDen[iS] = new TH1F("hDen_"+nameDS[iS], "Denominator ("+titleDS[iS]+")",
			     nBinsMet-1, 50, 2050);
	hNum[iS] = new TH1F("hNum_"+nameDS[iS], "Numerator ("+titleDS[iS]+")",
			     nBinsMet-1, 50, 2050);
      }

      else {
	hDen[iS] = new TH1F("hDen_"+nameDS[iS], "Denominator ("+titleDS[iS]+")",
			     nBinsMet-1, bins_met);
	hNum[iS] = new TH1F("hNum_"+nameDS[iS], "Numerator ("+titleDS[iS]+")",
			     nBinsMet-1, bins_met);
      }

      hDen[iS]->SetXTitle("PFMET NoMu (GeV)");
      hNum[iS]->SetXTitle("PFMET NoMu (GeV)");

      // Process Data Sets //
      cout << endl << "-- Processing sample: " 
	   << nameDS[iS] << " " << titleDS[iS] << " --" << endl;
    
      chain[iS] = new TChain("tree");
      for(UInt_t iDS=0; iDS<dirDS[iS].size(); iDS++) {
	if(mcDS[iS]) subDir = subDirMC;
	else         subDir = subDirData;
	fullPath = baseDir + subDir + dirDS[iS][iDS] + suffixDir;
	cout << "--- add file: " << fullPath << endl;
	chain[iS]->Add(fullPath);
      }
      cout << endl;

      processChain(nEvents, chain[iS], nameDS[iS], mcDS[iS], hDen[iS], hNum[iS], doCutOnMetFlags, doPrintBeforeCuts);
      if( !fout->IsZombie() ) fout->cd();

      cout << "-- writing histogram: " << hDen[iS]->GetName() << " (" << hDen[iS]->GetEntries() << " entries)" << endl;
      hDen[iS]->Sumw2();
      hDen[iS]->Write();

      cout << "-- writing histogram: " << hNum[iS]->GetName() << " (" << hNum[iS]->GetEntries() << " entries)" << endl;
      hNum[iS]->Sumw2();
      hNum[iS]->Write();

      // Delete pointers
      if(hNum[iS]) delete hNum[iS];
      if(hDen[iS]) delete hDen[iS];
      delete chain[iS];

    }

    if( !fout->IsZombie() ) fout->Write();
    if( !fout->IsZombie() ) fout->Close();
    if(fout) delete fout;
  }

  // Plot efficiencies
  else {
    plotEfficiency(plotDir, version, era, nDS, nameDS, titleDS, colorDS, styleDS);
  }

  // END //
  return 0;
}


Int_t plotEfficiency(TString plotDir, TString version, TString era, const UInt_t nDS, 
		     TString* nameDS, TString* titleDS,Int_t* colorDS, Int_t* styleDS)
{

  // Open histogram file
  TFile*  fread[nDS];
  TString fpath[nDS];
  TH1F*   hDen[nDS];
  TH1F*   hNum[nDS];
  TEfficiency* teff[nDS];
  TString hTitle, theTitle, theAxes;

  // Prepare ratio plot
  TH1F* hRatio[nDS];  

  // Prepare canvas //
  TCanvas* c  = new TCanvas("c" , "c" , 100, 100, 600, 600);
  TCanvas* c1 = new TCanvas("c1", "c1", 100, 100, 600, 600);
  TCanvas* c2 = new TCanvas("c2", "c2", 100, 100, 600, 600);
  TCanvas* c3 = new TCanvas("c3", "c3", 100, 100, 600, 600);

  // Prepare legend for overlayed canvas //
  //TLegend *leg = new TLegend(0.58,0.58,0.73,0.68);
  TLegend *leg = new TLegend(0.43,0.58,0.73,0.68);
  leg->SetFillColor(kWhite);

  // Prepare legend for ratio plot
  //TLegend *leg2 = new TLegend(0.43,0.58,0.73,0.68);
  TLegend *leg2 = new TLegend(0.43,0.78,0.73,0.88);
  leg2->SetFillColor(kWhite);

  // Get histograms from individual root files
  for(UInt_t iS=0; iS<nDS; iS++) {

    fpath[iS] = plotDir + "/" + version + "/histos_"+era+"_"+nameDS[iS]+".root";
    fread[iS] = new TFile(fpath[iS] , "read");

    if(fread[iS]->IsZombie()) {
      cout << "-- ERROR: file not found (" << fpath[iS] << ")" << endl;
      continue;
    }

    hDen[iS] = (TH1F*) fread[iS]->Get("hDen_"+nameDS[iS]);
    hNum[iS] = (TH1F*) fread[iS]->Get("hNum_"+nameDS[iS]);

    if(!hDen[iS]) {
      cout << "-- ERROR: " << "hDen_"+nameDS[iS] << " not found! continue;" << endl;
      continue;
    }

    if(!hNum[iS]) {
      cout << "-- ERROR: " << "hNum_"+nameDS[iS] << " not found! continue;" << endl;
      continue;
    }    

  } // end loop: root files

  for(UInt_t iS=0; iS<nDS; iS++) {

    /// plot numerator and denominator for current process
    c1->cd();
    //
    hNum[iS]->SetLineColor(kRed);
    hNum[iS]->SetFillColor(kRed);
    hDen[iS]->SetLineColor(kBlue);
    //
    hDen[iS]->Draw("H");
    hNum[iS]->Draw("HSAME");
    //
    c1->Update();
    c1->Print(plotDir + "/" + version + "/hist_"+nameDS[iS]+".png", "png");

    /// plot efficiencies
    if( TEfficiency::CheckConsistency(*hNum[iS], *hDen[iS]) ) {

      teff[iS] = new TEfficiency(*hNum[iS], *hDen[iS]);
      theAxes   = ";PFMET NoMu (GeV);HLT PFMNoMu120 Efficiency";
      theTitle  = titleDS[iS] + theAxes;
      teff[iS]->SetTitle(theTitle);
      //teff[iS]->;
      teff[iS]->SetMarkerStyle(styleDS[iS]);
      teff[iS]->SetMarkerColor(colorDS[iS]);
      teff[iS]->SetLineColor  (colorDS[iS]);
      //teff[iS]->SetMarkerSize(1);

      c->cd();
      teff[iS]->Draw("");
      //
      c->Print(plotDir + "/" + version + "/eff_"+nameDS[iS]+".png", "png");
      c->Print(plotDir + "/" + version + "/eff_"+nameDS[iS]+".pdf", "pdf");
      gPad->SetLogx(kTRUE);
      c->Update();
      c->Print(plotDir + "/" + version + "/eff_"+nameDS[iS]+"_logx.png", "png");
      c->Print(plotDir + "/" + version + "/eff_"+nameDS[iS]+"_logx.pdf", "pdf");
      gPad->SetLogx(kFALSE);

      c2->cd();
      teff[iS]->SetTitle(theAxes);
      if(iS==0) teff[iS]->Draw("");
      else      teff[iS]->Draw("SAME");

      if(hNum[iS]->GetEntries()>0)
	leg->AddEntry( teff[iS] , titleDS[iS] , "L" );
      leg->Draw();
      c2->Update();

    } //endif: check TEff TH1 consistency

    /// plot efficiency ratios
    if(iS>0) {
      cout << "-- Preparing ratio plot (iS=" << iS << ")" << endl;
      if(hNum[0]) hRatio[iS] = (TH1F*) hNum[0]->Clone("hratio_"+nameDS[iS]);
      else continue;
      //
      //hRatio[iS]->Sumw2();
      //
      hRatio[iS]->SetTitle ("");
      hRatio[iS]->SetYTitle("HLT PFMNoMu120 Efficiency ratio");
      hRatio[iS]->SetStats(kFALSE);
      //
      hRatio[iS]->SetMarkerStyle(styleDS[iS]);
      hRatio[iS]->SetMarkerColor(colorDS[iS]);
      hRatio[iS]->SetLineColor  (colorDS[iS]);
      //
      hRatio[iS]->Divide  (hDen[0]  );
      hRatio[iS]->Multiply(hDen[iS]);
      hRatio[iS]->Divide  (hNum[iS]);
      //
      c3->cd();
      //
      if(iS==1) hRatio[iS]->Draw("P");
      else      hRatio[iS]->Draw("PSAME");
      //
      if(hRatio[iS]->GetEntries()>0)
	leg2->AddEntry( hRatio[iS] , titleDS[iS] , "L" );
      leg2->Draw();
      //
      c3->Update();
      //

    }

  } // end loop: samples

  c2->cd();
  c2->Print(plotDir + "/" + version + "/eff_"+era+"_overlayed.png", "png");
  c2->Print(plotDir + "/" + version + "/eff_"+era+"_overlayed.pdf", "pdf");
  gPad->SetLogx(kTRUE);
  c2->Update();
  c2->Print(plotDir + "/" + version + "/eff_"+era+"_overlayed_logx.png", "png");
  c2->Print(plotDir + "/" + version + "/eff_"+era+"_overlayed_logx.pdf", "pdf");
  gPad->SetLogx(kFALSE);

  c3->cd();
  c3->Print(plotDir + "/" + version + "/ratio_"+era+"_overlayed.png", "png");
  c3->Print(plotDir + "/" + version + "/ratio_"+era+"_overlayed.pdf", "pdf");
  gPad->SetLogx(kTRUE);
  c3->Update();
  c3->Print(plotDir + "/" + version + "/ratio_"+era+"_overlayed_logx.png", "png");
  c3->Print(plotDir + "/" + version + "/ratio_"+era+"_overlayed_logx.pdf", "pdf");
  gPad->SetLogx(kFALSE);

  for(UInt_t iS=0; iS<nDS; iS++) {
    delete teff[iS];
    if(hRatio[iS] && iS>0) 
      delete hRatio[iS];
    delete hNum[iS];
    delete hDen[iS];    
    fread[iS]->Close();  
    delete fread[iS];
  }

  delete c;
  delete c1;
  delete c2;
  delete c3;
  delete leg;
  delete leg2;

  // End
  return 0;
}


Int_t processChain(Int_t nEvents, TChain* fChain, TString nameDS, Bool_t isMC, 
		   TH1F* hDen, TH1F* hNum, Bool_t doCutOnMetFlags, Bool_t doPrintBeforeCuts)
{

  // Log file
  ofstream outlog("log_trigger_"+nameDS+".txt", ios::out); 
  ofstream outmap("map_trigger_"+nameDS+".txt", ios::out); 

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
  Int_t           HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
  Int_t           HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v;
  Int_t           HLT_BIT_HLT_IsoMu24_v;
  Int_t           HLT_BIT_HLT_IsoMu24_eta2p1_v;
  Int_t           HLT_BIT_HLT_IsoMu27_v;
  Int_t           HLT_BIT_HLT_Ele32_WPTight_Gsf_v;
  Int_t           HLT_BIT_HLT_Ele35_WPTight_Gsf_v;
  Int_t           HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v;
  //
  Int_t           HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  Int_t           HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
  Int_t           HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
  Int_t           HLT_BIT_HLT_IsoMu24;
  Int_t           HLT_BIT_HLT_IsoMu24_eta2p1;
  Int_t           HLT_BIT_HLT_IsoMu27;
  Int_t           HLT_BIT_HLT_Ele32_WPTight_Gsf;
  Int_t           HLT_BIT_HLT_Ele35_WPTight_Gsf;
  Int_t           HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT;
  //
  Int_t           HLT_SingleEl;
  Int_t           HLT_MET;
  Int_t           HLT_SingleMu;
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
  const UInt_t    NLEP=3;
  Int_t           LepGood_pdgId[NLEP];   //[nLepGood]   
  Float_t         LepGood_pt[NLEP];   //[nLepGood]      
  Float_t         LepGood_eta[NLEP];   //[nLepGood]     
  Float_t         LepGood_phi[NLEP];   //[nLepGood]     
  Int_t           LepGood_triggerMatched[NLEP];   //[nLepGood]                             
  Int_t           LepGood_charge[NLEP];   //[nLepGood]  
  Int_t           LepGood_tightId[NLEP];   //[nLepGood] 
  Int_t           LepGood_hltId[NLEP];   //[nLepGood]   
  Int_t           LepGood_eleCutIdSpring15_25ns_v1[NLEP];   //[nLepGood]                   
  Int_t           LepGood_hasGainSwitchFlag[NLEP];   //[nLepGood]                          
  Float_t         LepGood_relIso03[NLEP];   //[nLepGood]
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

  if(isMC) {
    fChain->SetBranchAddress("xsec", &xsec); // , &b_xsec);
    fChain->SetBranchAddress("puWeight", &puWeight); // , &b_puWeight);
    fChain->SetBranchAddress("nTrueInt", &nTrueInt); // , &b_nTrueInt);
    fChain->SetBranchAddress("genWeight", &genWeight); // , &b_genWeight);
    fChain->SetBranchAddress("nVert", &nVert); // , &b_nVert);
  }

  //
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v", &HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", &HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v", &HLT_BIT_HLT_IsoMu24_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1_v", &HLT_BIT_HLT_IsoMu24_eta2p1_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v", &HLT_BIT_HLT_IsoMu27_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_WPTight_Gsf_v", &HLT_BIT_HLT_Ele32_WPTight_Gsf_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele35_WPTight_Gsf_v", &HLT_BIT_HLT_Ele35_WPTight_Gsf_v); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v", &HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v); 
  //
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight); 
  fChain->SetBranchAddress("HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24", &HLT_BIT_HLT_IsoMu24); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu24_eta2p1", &HLT_BIT_HLT_IsoMu24_eta2p1); 
  fChain->SetBranchAddress("HLT_BIT_HLT_IsoMu27", &HLT_BIT_HLT_IsoMu27); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele32_WPTight_Gsf", &HLT_BIT_HLT_Ele32_WPTight_Gsf); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele35_WPTight_Gsf", &HLT_BIT_HLT_Ele35_WPTight_Gsf); 
  fChain->SetBranchAddress("HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT); 
  //
  fChain->SetBranchAddress("HLT_MET", &HLT_MET); // , &b_HLT_MET);
  fChain->SetBranchAddress("HLT_SingleMu", &HLT_SingleMu); // , &b_HLT_SingleMu);
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

  if(isMC) {
    fChain->SetBranchAddress("met_genPt", &met_genPt); // , &b_met_genPt);
    fChain->SetBranchAddress("met_genPhi", &met_genPhi); // , &b_met_genPhi);
    fChain->SetBranchAddress("met_genEta", &met_genEta); // , &b_met_genEta);
  }

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
  fChain->SetBranchAddress("LepGood_pdgId", LepGood_pdgId); // , &b_LepGood_pdgId);
  fChain->SetBranchAddress("LepGood_pt", LepGood_pt); // , &b_LepGood_pt);
  fChain->SetBranchAddress("LepGood_eta", LepGood_eta); // , &b_LepGood_eta);
  fChain->SetBranchAddress("LepGood_phi", LepGood_phi); // , &b_LepGood_phi);
  fChain->SetBranchAddress("LepGood_triggerMatched", LepGood_triggerMatched); // , &b_LepGood_triggerMatched);
  fChain->SetBranchAddress("LepGood_charge", LepGood_charge); // , &b_LepGood_charge);
  fChain->SetBranchAddress("LepGood_tightId", LepGood_tightId); // , &b_LepGood_tightId);
  fChain->SetBranchAddress("LepGood_hltId", LepGood_hltId); // , &b_LepGood_hltId);
  fChain->SetBranchAddress("LepGood_eleCutIdSpring15_25ns_v1", LepGood_eleCutIdSpring15_25ns_v1); // , &b_LepGood_eleCutIdSpring15_25ns_v1);
  fChain->SetBranchAddress("LepGood_hasGainSwitchFlag", LepGood_hasGainSwitchFlag); // , &b_LepGood_hasGainSwitchFlag);
  fChain->SetBranchAddress("LepGood_relIso03", LepGood_relIso03); // , &b_LepGood_relIso03);
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
  //

  // Utility variables //
  UInt_t  iMax=0;
  Float_t ptMax=0;
  Float_t dphi_jme=-999;

  // Counters
  const UInt_t nEVC = 9;
  UInt_t  nEventCount[nEVC];
  for(UInt_t i=0; i<nEVC; i++) nEventCount[i]=0;
  TString nameStep[nEVC] = { "Muon Trigger" , "Muon Offline", "MET Flags" , "Tau Veto" , 
			     "MET Cut" , "#Jets" , "#Tracks" , "Leading Jet" , "DeltaPhi(Jet,MET)" }; 

  Int_t   foundPathologicalJet = -1;
  Bool_t leadJetIsFwd = kFALSE;
  Bool_t passDeltaPhi = kTRUE;
  Bool_t passMuon = kFALSE;

  // < <index in array, isForward>
  vector< pair<UInt_t,Bool_t> > idxJetID;
  vector<UInt_t> idxTrackID;

  // Loop over the chain //
  UInt_t nEntries   = fChain->GetEntries();
  UInt_t nToProcess = nEvents;
  if(nEvents<0 || nEvents>(Int_t)nEntries) nToProcess = nEntries;

  cout << "- About to process: " << nToProcess << " entries." << endl;
  for(UInt_t iE=0; iE<nToProcess; iE++) {

    // Initialize HLT bit variables
    HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v = 0 ;
    HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = 0 ;
    HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = 0 ;
    HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v = 0 ;
    HLT_BIT_HLT_IsoMu24_v = 0 ;
    HLT_BIT_HLT_IsoMu24_eta2p1_v = 0 ;
    HLT_BIT_HLT_IsoMu27_v = 0 ;
    HLT_BIT_HLT_Ele32_WPTight_Gsf_v = 0 ;
    HLT_BIT_HLT_Ele35_WPTight_Gsf_v = 0 ;
    HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT_v = 0 ;
    //
    HLT_BIT_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 = 0 ;
    HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = 0 ;
    HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = 0 ;
    HLT_BIT_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = 0 ;
    HLT_BIT_HLT_IsoMu24 = 0 ;
    HLT_BIT_HLT_IsoMu24_eta2p1 = 0 ;
    HLT_BIT_HLT_IsoMu27 = 0 ;
    HLT_BIT_HLT_Ele32_WPTight_Gsf = 0 ;
    HLT_BIT_HLT_Ele35_WPTight_Gsf = 0 ;
    HLT_BIT_HLT_Ele115_CaloIdVT_GsfTrkIdT = 0 ;

    // Get entry #iE    
    fChain->GetEntry(iE);
    if(iE%10000==0) {
      cout << "-- Processing entry #" << iE << " / " << nEntries << endl;
    }

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

    // 4- Level-1 Selection Cuts //

    //// Single lepton selection
    if(!HLT_SingleMu) continue;
    nEventCount[0]++ ;

    ///// FIXME: add electron veto for singlemu data
    passMuon = kFALSE;
    for(UInt_t iL=0; iL<NLEP && (Int_t)iL<nLepGood; iL++) { 
      if(TMath::Abs(LepGood_pdgId[iL])==13 && LepGood_pt[iL]>25)
	passMuon = kTRUE;
    }
    if(!passMuon) continue;
    nEventCount[1]++ ;

    //// MET flags
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
    nEventCount[2]++ ;
    //// 

    //// Lepton veto
    //if(nLepGood>0) continue;  
    //nEventCount[2]++ ;

    if(nTauGood>0) continue;
    nEventCount[3]++ ;

    //// MET Cut
    //if(metNoMu_pt<=200)   continue;
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
    //if(!passDeltaPhi) continue;
    nEventCount[8]++ ;
    ////

    //// Map selected events up to the previous step
    mapEvents[lumi].push_back(evt);

    // Printouts after cuts //
    for(UInt_t iJ=0; iJ<NJET && (Int_t)iJ<nJet ; iJ++) {
      if(Jet_pt[iJ]>2000) {
	cout << "Entry #" << iE << " Jet_pt[" << iJ << "] = " << Jet_pt[iJ] << endl;
      }
    }


    // 5- Fill Histograms for trigger efficiencies
    hDen->Fill(metNoMu_pt);
    if(HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v) 
      hNum->Fill(metNoMu_pt);

  } // end loop: entries


  // Printout # accepted events
  cout << endl << "----------------" << endl << "ACCEPTED EVENTS" << endl;
  for(UInt_t i=0; i<nEVC; i++) 
    cout << "#Events passing step #" << i << " = " << nEventCount[i] 
	 << "\t(" << nameStep[i] << ")" << endl;

  // Printout # rejected events
  // cout << endl << "----------------" << endl << "REJECTED EVENTS" << endl;
  // for(UInt_t i=0; i<nEVC; i++) 
  //   cout << "#Events rejected by step #" << i << " = " << nEntries - nEventCount[i] 
  // 	 << "\t(" << nameStep[i] << ")" << endl;

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
