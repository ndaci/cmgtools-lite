##########################################################
##       CONFIGURATION FOR SUSY STOP SOFT B TREES       ##
##########################################################
import PhysicsTools.HeppyCore.framework.config as cfg
import re

#-------- LOAD ALL ANALYZERS -----------
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
from CMGTools.RootTools.samples.configTools import *
from CMGTools.RootTools.samples.autoAAAconfig import *

from CMGTools.TTHAnalysis.analyzers.xtracks_modules_cff import *

lepAna.loose_muon_isoCut     = lambda muon : muon.relIso03 < 0.25
lepAna.loose_electron_isoCut = lambda elec : elec.relIso03 < 0.25

tauAna.loose_ptMin = 20
tauAna.loose_etaMax = 2.3
tauAna.loose_tauID = "decayModeFindingNewDMs"
tauAna.loose_vetoLeptons = False 

jetAna.addJECShifts = True
jetAna.jetPtOrUpOrDnSelection = True

jetAna.copyJetsByValue = True 
jetAna.calculateType1METCorrection = True
jetAnaScaleUp.calculateType1METCorrection = True
jetAnaScaleDown.calculateType1METCorrection = True

## early skimming on the isolated track
from CMGTools.TTHAnalysis.analyzers.ttHFastJetSkimmer import ttHFastJetSkimmer
fastJetSkim = cfg.Analyzer(ttHFastJetSkimmer, name="fastJetSkimmer",
    jets = 'slimmedJets',
    jetCut = lambda j : j.pt() > 80 and abs(j.eta()) < 2.4, # looser pt because of non-final JECs
    minJets = 1,
    )
from CMGTools.TTHAnalysis.analyzers.isoTrackFastSkimmer import isoTrackFastSkimmer
isoTrackFastSkim = cfg.Analyzer(isoTrackFastSkimmer, name="isoTrackFastSkim",
    cut = lambda t : (t.pt() > 50 and
                      abs(t.eta()) < 2.4 and
                      t.isHighPurityTrack() and
                      abs(t.dxy()) < 0.5 and abs(t.dz()) < 0.5 and
                      (t.miniPFIsolation().chargedHadronIso() < 1.0*t.pt() or t.pt() > 100))
)

## late skimming on jets and MET
from CMGTools.TTHAnalysis.analyzers.xtracksFilters import xtracksFilters
jetMETSkim = cfg.Analyzer( xtracksFilters, name='xtracksSkim',
    jets       = "cleanJets",
    jetPtCuts  = [ 90., ],
    metCut     =  0,
    metNoMuCut =  80,
)

## Full DeDx analyzer
from CMGTools.TTHAnalysis.analyzers.isoTrackDeDxAnalyzer import isoTrackDeDxAnalyzer
isoTrackDeDxAna = cfg.Analyzer(isoTrackDeDxAnalyzer, name="isoTrackDeDxAna",
    doDeDx = "94XMiniAODv1-Hack", 
        # for 94X MiniAOD v2, just set it to True
        # for 94X MiniAOD v1, you have two options
        #  - set it to False, and have no DeDx
        #  - set it to "94XMiniAODv1-Hack" and follow step (1) of https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3586/1/1/1/1.html 
    )

## Tree Producer
from CMGTools.TTHAnalysis.analyzers.treeProducerXtracks import *

## Sample production and setup
### Trigger
from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import *
triggerFlagsAna.triggerBits = {
    'SingleMu' : triggers_1mu_iso,
    'SingleEl' : triggers_1e_iso + triggers_1e_noniso,
    'MET'      : triggers_SOS_highMET
}
triggerFlagsAna.unrollbits = True

### MC
from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17MiniAOD import *
mcSamples = [ W3JetsToLNu_LO ]
autoAAA(mcSamples)
for c in mcSamples:
    c.triggers = triggers_1mu_iso[:]
    c.triggers += triggers_1e_iso + triggers_1e_noniso
    c.triggers += triggers_SOS_highMET

## Data
from CMGTools.RootTools.samples.samples_13TeV_DATA2017 import *
selectedComponents = mcSamples # + dataSamples

#-------- SEQUENCE -----------
sequence = cfg.Sequence( [
    lheWeightAna,
    pileUpAna,
    skimAnalyzer,
    jsonAna,
    triggerAna,

    fastJetSkim,
    isoTrackFastSkim,

    genAna,
    genHFAna,

    vertexAna,
    lepAna,
    tauAna,
    jetAna,
    jetAnaScaleUp,
    jetAnaScaleDown,
    metAna,
    metAnaScaleUp,
    metAnaScaleDown,
    jetMETSkim,

    isoTrackDeDxAna,

    triggerFlagsAna,
    eventFlagsAna,

    treeProducer,
])

#-------- HOW TO RUN -----------
test = getHeppyOption('test')
if test == "1":
    print "The test wil use %s " % selectedComponents[0].name
    selectedComponents = doTest1(selectedComponents[0], sequence=sequence, cache=True )
    print "The test wil use file %s " % selectedComponents[0].files[0]
elif test == "1S":
    comp = selectedComponents[0]
    comp.name = "Signal"
    comp.files = [ '/media/Disk1/avartak/CMS/Samples/DisappearingTracks/Wino_M_300_cTau_10/Wino_M_300_cTau_10.run17.ev17000.MiniAODv2.root' ]
    selectedComponents = doTest1(comp, sequence=sequence, cache=False )
    print "The test wil use file %s " % comp.files[0]
    fastJetSkim.minJets = 0
    fastMETSkim.metCut = 0
    isoTrackDeDxAna.doDeDx = True
    comp.triggers = []
elif test in ('2','3','5s'):
    doTestN(test,selectedComponents)

printSummary(selectedComponents)

config = autoConfig(selectedComponents, sequence) #, xrd_aggressive=-1)
