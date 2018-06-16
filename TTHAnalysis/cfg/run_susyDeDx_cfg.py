##########################################################
##       CONFIGURATION FOR SUSY STOP SOFT B TREES       ##
##########################################################
import PhysicsTools.HeppyCore.framework.config as cfg
import re, sys

#-------- LOAD ALL ANALYZERS -----------
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
from CMGTools.RootTools.samples.configTools import *
from CMGTools.RootTools.samples.autoAAAconfig import *

from CMGTools.TTHAnalysis.analyzers.xtracks_modules_cff import *

## Configuration that depends on the region
region = getHeppyOption("region","sr") # use 'sr' or "cr1l"

fastJetSkim.minJets = 1 if region =="sr" else 0
eventSkim.region = region

## Sample production and setup
### Trigger
from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import *
triggerFlagsAna.triggerBits = {
    'SingleMu' : triggers_1mu_iso,
    'SingleEl' : triggers_1e_iso + triggers_1e_noniso,
    'MET'      : triggers_SOS_highMET
}

### MC
from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17MiniAOD import *
Top = [ TTLep_pow, TTSemi_pow, T_tch, TBar_tch, T_tWch_noFullyHad, TBar_tWch_noFullyHad ]

Wino_M_300_cTau_3  = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_3",  "Wino_M_300_cTau_3",  "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_10", "Wino_M_300_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_30 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_30", "Wino_M_300_cTau_30", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_500_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_500_cTau_10", "Wino_M_500_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Winos = [ Wino_M_300_cTau_3  , Wino_M_300_cTau_10 , Wino_M_300_cTau_30 , Wino_M_500_cTau_10 ]

if region == "sr":   
    mcSamples = ([ W3JetsToLNu_LO, W1JetsToLNu_LO, W2JetsToLNu_LO, W4JetsToLNu_LO ]
                 + DYJetsToLLM50HT + DYJetsToLLM4to50HT
                 + ZvvLOHT 
                 + Top )
    mcSignals = Winos
    mcTriggers = triggers_SOS_highMET[:] 
elif region == "cr1l": 
    mcSamples = [ DYJetsToLL_M50, WJetsToLNu_LO ] + Top
    mcTriggers = triggers_1mu_iso + triggers_1e_iso + triggers_1e_noniso
    mcSignals = []

autoAAA(mcSamples)
cropToLumi(mcSamples, 10*41.7)
for c in mcSamples + mcSignals:
    c.triggers = mcTriggers

## Data
from CMGTools.RootTools.samples.samples_13TeV_DATA2017 import *
if region == "sr":   
    datasetsAndTriggers = [ ("MET", triggers_SOS_highMET) ]
elif region == "cr1l":   
    datasetsAndTriggers = [ ("SingleMuon", triggers_1mu_iso),
                            ("SingleElectron", triggers_1e_iso + triggers_1e_noniso) ]
json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
dataSamples = []; vetoTriggers = []
for (pdname, trigs) in datasetsAndTriggers:
    for d in dataSamples_17Nov2017:
        if pdname in d.name:
            d.json = json
            d.triggers = trigs[:]
            d.vetoTriggers = vetoTriggers[:]
            dataSamples.append(d)
    vetoTriggers += trigs

run = getHeppyOption("run","all")
if run == "all":    selectedComponents = mcSamples + dataSamples + mcSignals
elif run == "data": selectedComponents = dataSamples
elif run == "mc":   selectedComponents = mcSamples
elif run == "sig":  selectedComponents = mcSignals


#-------- SEQUENCE -----------
sequence = cfg.Sequence( xtracks_sequence )

#-------- HOW TO RUN -----------
test = getHeppyOption('test')
if test == "1":
    print "The test wil use %s " % selectedComponents[0].name
    selectedComponents = doTest1(selectedComponents[0], sequence=sequence, cache=True )
    print "The test wil use file %s " % selectedComponents[0].files[0]
elif test == "1S":
    comp = selectedComponents[0]
    comp.name = "Signal"
    comp.files = [ '/eos/cms/store/cmst3/user/gpetrucc/SusyWithDeDx/Wino_M_500_cTau_10.merged/Wino_M_500_cTau_10.MiniAODv2_job1.root' ]
    selectedComponents = doTest1(comp, sequence=sequence, cache=False )
    print "The test wil use file %s " % comp.files[0]
    fastJetSkim.minJets = 0
    isoTrackDeDxAna.doDeDx = True
    comp.triggers = []
elif test in ('2','3','5s'):
    doTestN(test,selectedComponents)

printSummary(selectedComponents)
if getHeppyOption("justSummary",False): sys.exit(0)
config = autoConfig(selectedComponents, sequence)
