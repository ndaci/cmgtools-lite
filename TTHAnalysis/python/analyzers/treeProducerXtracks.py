#!/bin/env python
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer import *  
from PhysicsTools.Heppy.analyzers.core.autovars import *  
from PhysicsTools.Heppy.analyzers.objects.autophobj import  *
from math import *

##------------------------------------------  
## LEPTON
##------------------------------------------  

leptonTypeXtracks = NTupleObjectType("leptonXtracks", baseObjectTypes = [ leptonType ], variables = [
    NTupleVariable("etaSc", lambda x : x.superCluster().eta() if abs(x.pdgId())==11 else -100, help="Electron supercluster pseudorapidity"),
    NTupleVariable("mcMatchPdgId",  lambda x : x.mcLep.pdgId() if getattr(x,'mcLep',None)!=None else -99, int, mcOnly=True, help="Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z): pdgId of the matched gen-level lepton, zero if non-prompt or fake"),
    NTupleVariable("mcPromptGamma", lambda x : x.mcPho.isPromptFinalState() if getattr(x,"mcPho",None) else 0, int, mcOnly=True, help="Photon isPromptFinalState"),
    NTupleVariable("triggerMatched", lambda x : 1 if x.matchedTrgObjSingleLep else 0, int, mcOnly=False, help="Is the lepton matched to a single lepton trigger"),
])
leptonTypeXtracks.removeVariable("relIsoAn04")
leptonTypeXtracks.removeVariable("eleCutIdSpring15_25ns_v1")
leptonTypeXtracks.removeVariable("ICHEPsoftMuonId")
leptonTypeXtracks.removeVariable("ICHEPmediumMuonId")

##------------------------------------------  
## TAU
##------------------------------------------  

tauTypeXtracks = NTupleObjectType("tauXtracks",  baseObjectTypes = [ tauType ], variables = [
     NTupleVariable("idMVAdR03", lambda x : x.idMVAdR03, int, help="1,2,3,4,5,6 if the tau passes the very loose to very very tight WP of the IsolationMVArun2v1DBdR03oldDMwLT discriminator"),
])
tauTypeXtracks.removeVariable("idMVANewDM")
tauTypeXtracks.removeVariable("idCI3hit")


jetTypeXtracks = NTupleObjectType("jetXtracks",  baseObjectTypes = [ jetTypeExtra ], variables = [
    NTupleVariable("chHEF", lambda x : x.chargedHadronEnergyFraction(), float, mcOnly = False, help="chargedHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neHEF", lambda x : x.neutralHadronEnergyFraction(), float, mcOnly = False,help="neutralHadronEnergyFraction (relative to uncorrected jet energy)"),
])
jetTypeXtracks.removeVariable("btagCMVA")
jetTypeXtracks.removeVariable("nLeptons")

jetTypeXtracksFwd = NTupleObjectType("jetFwd",  baseObjectTypes = [ jetType ], variables = [
    NTupleVariable("chHEF", lambda x : x.chargedHadronEnergyFraction(), float, mcOnly = False, help="chargedHadronEnergyFraction (relative to uncorrected jet energy)"),
    NTupleVariable("neHEF", lambda x : x.neutralHadronEnergyFraction(), float, mcOnly = False,help="neutralHadronEnergyFraction (relative to uncorrected jet energy)"),
])
jetTypeXtracks.removeVariable("btagCMVA")
jetTypeXtracksFwd.removeVariable("nLeptons")

##------------------------------------------  
## MET
##------------------------------------------  
  
metTypeXtracks = NTupleObjectType("metXtracks", baseObjectTypes = [ metType ], variables = [
])

metTypeXtracksBasic = NTupleObjectType("metXtracksBasic", baseObjectTypes = [ fourVectorType ], variables = [
])

##------------------------------------------  
## IsoTrackDeDx
##------------------------------------------

def _isoDBeta(x):
    return x.chargedHadronIso() + max(x.neutralHadronIso()+x.photonIso()-0.5*x.puChargedHadronIso(), 0)

hitPatternType = NTupleObjectType("hitPatternType", baseObjectTypes = [], variables = [
    NTupleVariable("trackerLayers", lambda x : x.hitPattern().trackerLayersWithMeasurement(), int, help="Tracker Layers"),
    NTupleVariable("pixelLayers", lambda x : x.hitPattern().pixelLayersWithMeasurement(), int, help="Pixel Layers"),
    NTupleVariable("trackerHits", lambda x : x.hitPattern().numberOfValidTrackerHits(), int, help="Tracker hits"),
    NTupleVariable("pixelHits", lambda x : x.hitPattern().numberOfValidPixelHits(), int, help="Pixel hits"),
    NTupleVariable("missingInnerPixelHits", lambda x : x.hitPattern().numberOfLostPixelHits(1), int, help="Missing inner pixel hits"),
    NTupleVariable("missingOuterPixelHits", lambda x : x.hitPattern().numberOfLostPixelHits(2), int, help="Missing outer pixel hits"),
    NTupleVariable("missingInnerStripHits", lambda x : x.hitPattern().numberOfLostStripHits(1), int, help="Missing inner strips hits"),
    NTupleVariable("missingOuterStripHits", lambda x : x.hitPattern().numberOfLostStripHits(2), int, help="Missing outer strips hits"),
    NTupleVariable("missingInnerTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(1), int, help="Missing inner tracker hits"),
    NTupleVariable("missingOuterTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(2), int, help="Missing outer tracker hits"),
    NTupleVariable("missingMiddleTrackerHits", lambda x : x.hitPattern().numberOfLostTrackerHits(0), int, help="Missing tracker hits in the middle of the track"),
])
hitDeDxType = NTupleObjectType("hitDeDxType", baseObjectTypes = [], variables = [
    NTupleVariable("dedxByLayer0", lambda x : x.dedxByLayer[0]),
    NTupleVariable("dedxByLayer1", lambda x : x.dedxByLayer[1]),
    NTupleVariable("dedxByLayer2", lambda x : x.dedxByLayer[2]),
    NTupleVariable("dedxByLayer3", lambda x : x.dedxByLayer[3]),
    NTupleVariable("dedxByLayer4", lambda x : x.dedxByLayer[4]),
    NTupleVariable("dedxByLayer5", lambda x : x.dedxByLayer[5]),
    NTupleVariable("dedxByLayer6", lambda x : x.dedxByLayer[6]),
    NTupleVariable("dedxByLayer7", lambda x : x.dedxByLayer[7]),
    NTupleVariable("dedxByLayer8", lambda x : x.dedxByLayer[8]),
    NTupleVariable("dedxByLayer9", lambda x : x.dedxByLayer[9]),
    NTupleVariable("dedxByLayer10", lambda x : x.dedxByLayer[10]),
    NTupleVariable("dedxByLayer11", lambda x : x.dedxByLayer[11]),
    NTupleVariable("dedxByLayer12", lambda x : x.dedxByLayer[12]),
    NTupleVariable("dedxByLayer13", lambda x : x.dedxByLayer[13]),

    NTupleVariable("dedxUnSmearedByLayer0", lambda x : x.dedxUnSmearedByLayer[0]),
    NTupleVariable("dedxUnSmearedByLayer1", lambda x : x.dedxUnSmearedByLayer[1]),
    NTupleVariable("dedxUnSmearedByLayer2", lambda x : x.dedxUnSmearedByLayer[2]),
    NTupleVariable("dedxUnSmearedByLayer3", lambda x : x.dedxUnSmearedByLayer[3]),
    NTupleVariable("dedxUnSmearedByLayer4", lambda x : x.dedxUnSmearedByLayer[4]),
    NTupleVariable("dedxUnSmearedByLayer5", lambda x : x.dedxUnSmearedByLayer[5]),
    NTupleVariable("dedxUnSmearedByLayer6", lambda x : x.dedxUnSmearedByLayer[6]),
    NTupleVariable("dedxUnSmearedByLayer7", lambda x : x.dedxUnSmearedByLayer[7]),
    NTupleVariable("dedxUnSmearedByLayer8", lambda x : x.dedxUnSmearedByLayer[8]),
    NTupleVariable("dedxUnSmearedByLayer9", lambda x : x.dedxUnSmearedByLayer[9]),
    NTupleVariable("dedxUnSmearedByLayer10", lambda x : x.dedxUnSmearedByLayer[10]),
    NTupleVariable("dedxUnSmearedByLayer11", lambda x : x.dedxUnSmearedByLayer[11]),
    NTupleVariable("dedxUnSmearedByLayer12", lambda x : x.dedxUnSmearedByLayer[12]),
    NTupleVariable("dedxUnSmearedByLayer13", lambda x : x.dedxUnSmearedByLayer[13]),
    
    NTupleVariable("subDetIdByLayer0", lambda x : x.subDetIdByLayer[0], int),
    NTupleVariable("subDetIdByLayer1", lambda x : x.subDetIdByLayer[1], int),
    NTupleVariable("subDetIdByLayer2", lambda x : x.subDetIdByLayer[2], int),
    NTupleVariable("subDetIdByLayer3", lambda x : x.subDetIdByLayer[3], int),
    NTupleVariable("subDetIdByLayer4", lambda x : x.subDetIdByLayer[4], int),
    NTupleVariable("subDetIdByLayer5", lambda x : x.subDetIdByLayer[5], int),
    NTupleVariable("subDetIdByLayer6", lambda x : x.subDetIdByLayer[6], int),
    NTupleVariable("subDetIdByLayer7", lambda x : x.subDetIdByLayer[7], int),
    NTupleVariable("subDetIdByLayer8", lambda x : x.subDetIdByLayer[8], int),
    NTupleVariable("subDetIdByLayer9", lambda x : x.subDetIdByLayer[9], int),
    NTupleVariable("subDetIdByLayer10", lambda x : x.subDetIdByLayer[10], int),
    NTupleVariable("subDetIdByLayer11", lambda x : x.subDetIdByLayer[11], int),
    NTupleVariable("subDetIdByLayer12", lambda x : x.subDetIdByLayer[12], int),
    NTupleVariable("subDetIdByLayer13", lambda x : x.subDetIdByLayer[13], int),

    NTupleVariable("sizeXbyLayer0", lambda x : x.sizeXbyLayer[0], int),
    NTupleVariable("sizeXbyLayer1", lambda x : x.sizeXbyLayer[1], int),
    NTupleVariable("sizeXbyLayer2", lambda x : x.sizeXbyLayer[2], int),
    NTupleVariable("sizeXbyLayer3", lambda x : x.sizeXbyLayer[3], int),
    NTupleVariable("sizeXbyLayer4", lambda x : x.sizeXbyLayer[4], int),
    NTupleVariable("sizeXbyLayer5", lambda x : x.sizeXbyLayer[5], int),
    NTupleVariable("sizeXbyLayer6", lambda x : x.sizeXbyLayer[6], int),
    NTupleVariable("sizeXbyLayer7", lambda x : x.sizeXbyLayer[7], int),
    NTupleVariable("sizeXbyLayer8", lambda x : x.sizeXbyLayer[8], int),
    NTupleVariable("sizeXbyLayer9", lambda x : x.sizeXbyLayer[9], int),
    NTupleVariable("sizeXbyLayer10", lambda x : x.sizeXbyLayer[10], int),
    NTupleVariable("sizeXbyLayer11", lambda x : x.sizeXbyLayer[11], int),
    NTupleVariable("sizeXbyLayer12", lambda x : x.sizeXbyLayer[12], int),
    NTupleVariable("sizeXbyLayer13", lambda x : x.sizeXbyLayer[13], int),
                                                                                                         
    NTupleVariable("sizeYbyLayer0", lambda x : x.sizeYbyLayer[0], int),
    NTupleVariable("sizeYbyLayer1", lambda x : x.sizeYbyLayer[1], int),
    NTupleVariable("sizeYbyLayer2", lambda x : x.sizeYbyLayer[2], int),
    NTupleVariable("sizeYbyLayer3", lambda x : x.sizeYbyLayer[3], int),
    NTupleVariable("sizeYbyLayer4", lambda x : x.sizeYbyLayer[4], int),
    NTupleVariable("sizeYbyLayer5", lambda x : x.sizeYbyLayer[5], int),
    NTupleVariable("sizeYbyLayer6", lambda x : x.sizeYbyLayer[6], int),
    NTupleVariable("sizeYbyLayer7", lambda x : x.sizeYbyLayer[7], int),
    NTupleVariable("sizeYbyLayer8", lambda x : x.sizeYbyLayer[8], int),
    NTupleVariable("sizeYbyLayer9", lambda x : x.sizeYbyLayer[9], int),
    NTupleVariable("sizeYbyLayer10", lambda x : x.sizeYbyLayer[10], int),
    NTupleVariable("sizeYbyLayer11", lambda x : x.sizeYbyLayer[11], int),
    NTupleVariable("sizeYbyLayer12", lambda x : x.sizeYbyLayer[12], int),
    NTupleVariable("sizeYbyLayer13", lambda x : x.sizeYbyLayer[13], int),
    
    NTupleVariable("layerOrSideByLayer0", lambda x : x.layerOrSideByLayer[0], int),
    NTupleVariable("layerOrSideByLayer1", lambda x : x.layerOrSideByLayer[1], int),
    NTupleVariable("layerOrSideByLayer2", lambda x : x.layerOrSideByLayer[2], int),
    NTupleVariable("layerOrSideByLayer3", lambda x : x.layerOrSideByLayer[3], int),
    NTupleVariable("layerOrSideByLayer4", lambda x : x.layerOrSideByLayer[4], int),
    NTupleVariable("layerOrSideByLayer5", lambda x : x.layerOrSideByLayer[5], int),
    NTupleVariable("layerOrSideByLayer6", lambda x : x.layerOrSideByLayer[6], int),
    NTupleVariable("layerOrSideByLayer7", lambda x : x.layerOrSideByLayer[7], int),
    NTupleVariable("layerOrSideByLayer8", lambda x : x.layerOrSideByLayer[8], int),
    NTupleVariable("layerOrSideByLayer9", lambda x : x.layerOrSideByLayer[9], int),
    NTupleVariable("layerOrSideByLayer10", lambda x : x.layerOrSideByLayer[10], int),
    NTupleVariable("layerOrSideByLayer11", lambda x : x.layerOrSideByLayer[11], int),
    NTupleVariable("layerOrSideByLayer12", lambda x : x.layerOrSideByLayer[12], int),
    NTupleVariable("layerOrSideByLayer13", lambda x : x.layerOrSideByLayer[13], int),

    NTupleVariable("ladderOrBladeByLayer0", lambda x : x.ladderOrBladeByLayer[0], int),
    NTupleVariable("ladderOrBladeByLayer1", lambda x : x.ladderOrBladeByLayer[1], int),
    NTupleVariable("ladderOrBladeByLayer2", lambda x : x.ladderOrBladeByLayer[2], int),
    NTupleVariable("ladderOrBladeByLayer3", lambda x : x.ladderOrBladeByLayer[3], int),
    NTupleVariable("ladderOrBladeByLayer4", lambda x : x.ladderOrBladeByLayer[4], int),
    NTupleVariable("ladderOrBladeByLayer5", lambda x : x.ladderOrBladeByLayer[5], int),
    NTupleVariable("ladderOrBladeByLayer6", lambda x : x.ladderOrBladeByLayer[6], int),
    NTupleVariable("ladderOrBladeByLayer7", lambda x : x.ladderOrBladeByLayer[7], int),
    NTupleVariable("ladderOrBladeByLayer8", lambda x : x.ladderOrBladeByLayer[8], int),
    NTupleVariable("ladderOrBladeByLayer9", lambda x : x.ladderOrBladeByLayer[9], int),
    NTupleVariable("ladderOrBladeByLayer10", lambda x : x.ladderOrBladeByLayer[10], int),
    NTupleVariable("ladderOrBladeByLayer11", lambda x : x.ladderOrBladeByLayer[11], int),
    NTupleVariable("ladderOrBladeByLayer12", lambda x : x.ladderOrBladeByLayer[12], int),
    NTupleVariable("ladderOrBladeByLayer13", lambda x : x.ladderOrBladeByLayer[13], int),

    NTupleVariable("sideByLayer0", lambda x : x.sideByLayer[0], int),
    NTupleVariable("sideByLayer1", lambda x : x.sideByLayer[1], int),
    NTupleVariable("sideByLayer2", lambda x : x.sideByLayer[2], int),
    NTupleVariable("sideByLayer3", lambda x : x.sideByLayer[3], int),
    NTupleVariable("sideByLayer4", lambda x : x.sideByLayer[4], int),
    NTupleVariable("sideByLayer5", lambda x : x.sideByLayer[5], int),
    NTupleVariable("sideByLayer6", lambda x : x.sideByLayer[6], int),
    NTupleVariable("sideByLayer7", lambda x : x.sideByLayer[7], int),
    NTupleVariable("sideByLayer8", lambda x : x.sideByLayer[8], int),
    NTupleVariable("sideByLayer9", lambda x : x.sideByLayer[9], int),
    NTupleVariable("sideByLayer10", lambda x : x.sideByLayer[10], int),
    NTupleVariable("sideByLayer11", lambda x : x.sideByLayer[11], int),
    NTupleVariable("sideByLayer12", lambda x : x.sideByLayer[12], int),
    NTupleVariable("sideByLayer13", lambda x : x.sideByLayer[13], int),



    #NTupleVariable("moduleByLayer0", lambda x : x.moduleByLayer[0], int),        # not used so far
    #NTupleVariable("moduleByLayer1", lambda x : x.moduleByLayer[1], int),        # not used so far
    #NTupleVariable("moduleByLayer2", lambda x : x.moduleByLayer[2], int),        # not used so far
    #NTupleVariable("moduleByLayer3", lambda x : x.moduleByLayer[3], int),        # not used so far
    #NTupleVariable("moduleByLayer4", lambda x : x.moduleByLayer[4], int),        # not used so far
    #NTupleVariable("moduleByLayer5", lambda x : x.moduleByLayer[5], int),        # not used so far
    #NTupleVariable("moduleByLayer6", lambda x : x.moduleByLayer[6], int),        # not used so far
    #NTupleVariable("moduleByLayer7", lambda x : x.moduleByLayer[7], int),        # not used so far
    #NTupleVariable("moduleByLayer8", lambda x : x.moduleByLayer[8], int),        # not used so far
    #NTupleVariable("moduleByLayer9", lambda x : x.moduleByLayer[9], int),        # not used so far
    #NTupleVariable("moduleByLayer10", lambda x : x.moduleByLayer[10], int),      # not used so far
    #NTupleVariable("moduleByLayer11", lambda x : x.moduleByLayer[11], int),      # not used so far
    #NTupleVariable("moduleByLayer12", lambda x : x.moduleByLayer[12], int),      # not used so far
    #NTupleVariable("moduleByLayer13", lambda x : x.moduleByLayer[13], int),      # not used so far

    NTupleVariable("pixByLayer0", lambda x : x.pixByLayer[0], int),
    NTupleVariable("pixByLayer1", lambda x : x.pixByLayer[1], int),
    NTupleVariable("pixByLayer2", lambda x : x.pixByLayer[2], int),
    NTupleVariable("pixByLayer3", lambda x : x.pixByLayer[3], int),
    NTupleVariable("pixByLayer4", lambda x : x.pixByLayer[4], int),
    NTupleVariable("pixByLayer5", lambda x : x.pixByLayer[5], int),
    NTupleVariable("pixByLayer6", lambda x : x.pixByLayer[6], int),
    NTupleVariable("pixByLayer7", lambda x : x.pixByLayer[7], int),
    NTupleVariable("pixByLayer8", lambda x : x.pixByLayer[8], int),
    NTupleVariable("pixByLayer9", lambda x : x.pixByLayer[9], int),
    NTupleVariable("pixByLayer10", lambda x : x.pixByLayer[10], int),
    NTupleVariable("pixByLayer11", lambda x : x.pixByLayer[11], int),
    NTupleVariable("pixByLayer12", lambda x : x.pixByLayer[12], int),
    NTupleVariable("pixByLayer13", lambda x : x.pixByLayer[13], int),


])

isoTrackTypeDeDx = NTupleObjectType("isoTrackTypeDeDx", baseObjectTypes = [ particleType, hitPatternType, hitDeDxType ], variables = [
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("dxy",   lambda x : x.dxy(), help="d_{xy} with respect to PV, in cm (with sign)"),
    NTupleVariable("dz",    lambda x : x.dz() , help="d_{z} with respect to PV, in cm (with sign)"),
    NTupleVariable("edxy",  lambda x : x.dxyError(), help="#sigma(d_{xy}) with respect to PV, in cm"),
    NTupleVariable("edz", lambda x : x.dzError(), help="#sigma(d_{z}) with respect to PV, in cm"),    


    NTupleVariable("highPurity", lambda x : x.isHighPurityTrack(), int, help="High purity"),

    NTupleVariable("miniIsoCH",   lambda x : x.miniPFIsolation().chargedHadronIso(), help="Charged hadron mini-isolation"),
    NTupleVariable("miniIsoNH",   lambda x : x.miniPFIsolation().neutralHadronIso(), help="Neutral hadron mini-isolation"),
    NTupleVariable("miniIsoPH",   lambda x : x.miniPFIsolation().photonIso(), help="Photon mini-isolation"),
    NTupleVariable("miniIsoPU",   lambda x : x.miniPFIsolation().puChargedHadronIso(), help="Pileup charged hadron mini-isolation"),
    NTupleVariable("miniRelIso",   lambda x : _isoDBeta(x.miniPFIsolation())/x.pt(), help="mini-isolation (relative, PU-corrected using dBeta)"),
    NTupleVariable("dR03IsoCH",   lambda x : x.pfIsolationDR03().chargedHadronIso(), help="Charged hadron isolation dR=0.3"),
    NTupleVariable("dR03IsoNH",   lambda x : x.pfIsolationDR03().neutralHadronIso(), help="Neutral hadron isolation dR=0.3"),
    NTupleVariable("dR03IsoPH",   lambda x : x.pfIsolationDR03().photonIso(), help="Photon isolation = dR=0.3"),
    NTupleVariable("dR03IsoPU",   lambda x : x.pfIsolationDR03().puChargedHadronIso(), help="Pileup charged hadron isolation dR=0.3"),
    NTupleVariable("relIso03",   lambda x : _isoDBeta(x.pfIsolationDR03())/x.pt(), help="relative isolation dR=0.3 (PU-corrected using dBeta)"),
    NTupleVariable("caloEmEnergy",   lambda x : x.matchedCaloJetEmEnergy(), help="Energy in the ECAL behind the track"),
    NTupleVariable("caloHadEnergy",   lambda x : x.matchedCaloJetHadEnergy(), help="Energy in the HCAL behind the track"),
    NTupleVariable("channelsGoodECAL", lambda x : x.channelsGoodECAL, int, help="Flag set to 1 when the track extrapolates to all good ECAL channels"),
    NTupleVariable("channelsGoodHCAL", lambda x : x.channelsGoodHCAL, int, help="Flag set to 1 when the track extrapolates to all good HCAL channels"),

    NTupleVariable("awayJet_idx", lambda x : x.leadAwayJet.index if x.leadAwayJet else -1, int),
    NTupleVariable("awayJet_pt", lambda x : x.leadAwayJet.pt() if x.leadAwayJet else 0),
    NTupleVariable("awayJet_dr", lambda x : deltaR(x, x.leadAwayJet) if x.leadAwayJet else 0),

    NTupleVariable("awayNJet", lambda x : x.awayJetInfo['num'], int),
    NTupleVariable("awayHTJet", lambda x : x.awayJetInfo['ht']),

    NTupleVariable("awayMu_idx", lambda x : x.leadAwayMu.index if x.leadAwayMu else -1, int),
    NTupleVariable("awayMu_dr", lambda x : deltaR(x, x.leadAwayMu) if x.leadAwayMu else 0),
    NTupleVariable("awayMu_mll", lambda x : (x.leadAwayMu.p4()+x.p4()).M() if x.leadAwayMu else 0),

    NTupleVariable("awayEle_idx", lambda x : x.leadAwayEle.index if x.leadAwayEle else -1, int),
    NTupleVariable("awayEle_dr", lambda x : deltaR(x, x.leadAwayEle) if x.leadAwayEle else 0),
    NTupleVariable("awayEle_mll", lambda x : (x.leadAwayEle.p4()+x.p4()).M() if x.leadAwayEle else 0),

    NTupleVariable("closestMu_idx", lambda x : x.closestMu.index if x.closestMu else -1, int),
    NTupleVariable("closestEle_idx", lambda x : x.closestEle.index if x.closestEle else -1, int),
    NTupleVariable("closestTau_idx", lambda x : x.closestTau.index if x.closestTau else -1, int),

    NTupleVariable("trigLepton_idx", lambda x : x.trigLepton.index if getattr(x, 'trigLepton', None) else -1, int),

    NTupleVariable("myDeDx", lambda x : x.myDeDx),
                                                                                                         
    NTupleVariable("mcMatch", lambda x : x.mcMatch.index if x.mcMatch else -1, int, mcOnly=True),
    NTupleVariable("mcMatchAnyId", lambda x : x.mcMatchAny.pdgId()*(1+99*x.mcMatchAny.isDirectPromptTauDecayProductFinalState()) if x.mcMatchAny else 0, int, mcOnly=True, help="MC pdgId of the matched gen lepton, tau, photon or chargino (for leptons from tau, it's pdgId*100)"),
    NTupleVariable("mcMatchAnyPt", lambda x : x.mcMatchAny.pt() if x.mcMatchAny else 0, int, mcOnly=True, help="MC pt of the matched gen lepton, tau, photon or chargino"),
]
)

calibTrackTypeDeDx = NTupleObjectType("calibTrackTypeDeDx", baseObjectTypes = [ hitPatternType, hitDeDxType ], variables = [
    NTupleVariable("pt",   lambda x : x.pt()),
    NTupleVariable("eta",   lambda x : x.eta()),
    NTupleVariable("phi",   lambda x : x.phi()),
    NTupleVariable("p",   lambda x : x.p()),
    NTupleVariable("charge",   lambda x : x.charge(), int),
    NTupleVariable("dxy",   lambda x : x.dxy, help="d_{xy} with respect to PV, in cm (with sign)"),
    NTupleVariable("dz",    lambda x : x.dz , help="d_{z} with respect to PV, in cm (with sign)"),
    NTupleVariable("edxy",  lambda x : x.dxyError(), help="#sigma(d_{xy}) with respect to PV, in cm"),
    NTupleVariable("edz", lambda x : x.dzError(), help="#sigma(d_{z}) with respect to PV, in cm"),    
    NTupleVariable("highPurity", lambda x : x.quality(x.highPurity), int, help="High purity"),
    NTupleVariable("dedxPrescale", lambda x : x.dedxPrescale, int, help="Prescale factor applied when selecting this track to save dE/dx"),
    NTupleVariable("stripDeDx", lambda x : x.stripDeDx.dEdx(), help="dE/dx in the strip detector (harmonic2 estimator)"),
    NTupleVariable("stripDeDxErr", lambda x : x.stripDeDx.dEdxError(), help="dE/dx uncertainty in the strip detector (harmonic2 estimator)"),
    NTupleVariable("stripDeDxN", lambda x : x.stripDeDx.numberOfMeasurements(), help="number of dE/dx measurements in the strip detector"),
])


##------------------------------------------  
## genCharginoType
##------------------------------------------  
genCharginoType  = NTupleObjectType("genCharginoType", baseObjectTypes = [ genParticleWithMotherId ], variables = [
    NTupleVariable("beta", lambda x : x.p()/x.energy()),
    NTupleVariable("decayR", lambda x : x.decayPoint.R()),
    NTupleVariable("decayZ", lambda x : x.decayPoint.Z()),
])

## Tree Producer
treeProducer = cfg.Analyzer(
    AutoFillTreeProducer, name='treeProducerXtracks',
    vectorTree = True, saveTLorentzVectors = False,  defaultFloatType = 'F', PDFWeights = [],
    globalVariables = [
        NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
        NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),

        NTupleVariable("nJet30", lambda ev: sum([j.pt() > 30 for j in ev.cleanJets]), int, help="Number of jets with pt > 30, |eta|<2.4"),
        NTupleVariable("nJet30a", lambda ev: sum([j.pt() > 30 for j in ev.cleanJetsAll]), int, help="Number of jets with pt > 30, |eta|<4.7"),

        ## ------- lheHT, needed for merging HT binned samples
        NTupleVariable("lheHT", lambda ev : getattr(ev,"lheHT",-999), mcOnly=True, help="H_{T} computed from quarks and gluons in Heppy LHEAnalyzer"),
        NTupleVariable("lheHTIncoming", lambda ev : getattr(ev,"lheHTIncoming",-999), mcOnly=True, help="H_{T} computed from quarks and gluons in Heppy LHEAnalyzer (only LHE status<0 as mothers)"),
        NTupleVariable("lheVpt", lambda ev : getattr(ev,"lheV_pt",-999), mcOnly=True, help="p_{T} of the V boson at LHE level"),
        ],
    globalObjects = {
        "met"             : NTupleObject("met", metTypeXtracks, help="PF E_{T}^{miss}"),
        "metNoMu"         : NTupleObject("metNoMu", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons"),
        "met_jecUp"       : NTupleObject("met_jecUp", metTypeXtracks, help="PF E_{T}^{miss}, after type 1 corrections (JEC plus 1sigma)"),
        "metNoMu_jecUp"   : NTupleObject("metNoMu_jecUp", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons, after type 1 corrections (JEC plus 1sigma)"),
        "met_jecDown"     : NTupleObject("met_jecDown", metTypeXtracks, help="PF E_{T}^{miss}, after type 1 corrections (JEC minus 1sigma)"),
        "metNoMu_jecDown" : NTupleObject("metNoMu_jecUp", metTypeXtracksBasic, help="PF E_{T}^{miss} without muons, after type 1 corrections (JEC minus 1sigma)"),
        },
    collections = {
        ##--------------------------------------------------
        "selectedTaus"    : NTupleCollection("TauGood",  tauTypeXtracks, 8, help="Taus after the preselection"),
        "selectedLeptons" : NTupleCollection("LepGood",  leptonTypeXtracks, 8, help="Leptons after the preselection"),
        ##------------------------------------------------
        "cleanJets"       : NTupleCollection("Jet",     jetTypeXtracks, 15, help="Cental jets after full selection and cleaning, sorted by pt"),
        "cleanJetsFwd"    : NTupleCollection("JetFwd",  jetTypeXtracksFwd,  6, help="Forward jets after full selection and cleaning, sorted by pt"),
        ##------------------------------------------------
        "isoTracks"       : NTupleCollection("IsoTrack",  isoTrackTypeDeDx, 4, help="Isolated tracks"),
        ##------------------------------------------------
        "genCharginos"    : NTupleCollection("GenChargino",  genCharginoType, 4, mcOnly=True, help="Gen chargino"),
        },
    )

## Tree Producer
calibTreeProducer = cfg.Analyzer(
    AutoFillTreeProducer, name='treeProducerXtracks',
    vectorTree = True, saveTLorentzVectors = False,  defaultFloatType = 'F', PDFWeights = [],
    globalVariables = [
        NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
        NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),
        ],
    collections = {
       "tracks" : NTupleCollection("Track",  calibTrackTypeDeDx, 4, help="Isolated tracks"),
       }
   )

