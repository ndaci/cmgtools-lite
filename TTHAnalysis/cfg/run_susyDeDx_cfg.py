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
region = getHeppyOption("region","sr") # use 'sr', 'cr1l'

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
from CMGTools.RootTools.samples.samples_13TeV_RunIIFall17MiniAODv2 import *

Wino_M_300_cTau_3  = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_3",  "Wino_M_300_cTau_3",  "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_10", "Wino_M_300_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_300_cTau_30 = kreator.makeMCComponentFromEOS("Wino_M_300_cTau_30", "Wino_M_300_cTau_30", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_500_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_500_cTau_10", "Wino_M_500_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_650_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_650_cTau_10", "Wino_M_650_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_800_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_800_cTau_10", "Wino_M_800_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Wino_M_1000_cTau_10 = kreator.makeMCComponentFromEOS("Wino_M_1000_cTau_10", "Wino_M_1000_cTau_10", "/store/cmst3/user/gpetrucc/SusyWithDeDx/%s.merged", ".*root", 1)
Winos = [ Wino_M_300_cTau_3  , Wino_M_300_cTau_10 , Wino_M_300_cTau_30 , Wino_M_500_cTau_10, Wino_M_650_cTau_10, Wino_M_800_cTau_10, Wino_M_1000_cTau_10 ]

Top = [ TTLep, TTHad, TTSemi] + Ts
VV  = [ WW, WZ, ZZ ]

Zll = [
    DYJetsM50_HT100to200,     DYJetsM50_HT100to200e,
    DYJetsM50_HT200to400,     DYJetsM50_HT200to400e,
    DYJetsM50_HT400to600,     DYJetsM50_HT400to600e,
    DYJetsM50_HT600to800,
    DYJetsM50_HT800to1200,
    DYJetsM50_HT1200to2500,
    DYJetsM50_HT2500toInf,
]

if region == "sr":   
    mcSamples =  ( QCD
                 + Ws
                 + Zvv
                 + Zll 
                 + VV
                 + Top)
    mcSignals = Winos
    mcTriggers = triggers_SOS_highMET[:] 
elif region == "cr1l": 
    mcSamples =  ( QCD
                 + Ws
                 + Zll 
                 + VV
                 + Top)
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
    for d in dataSamples_31Mar2018:
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
elif run == "aod":
    from CMGTools.RootTools.samples.samples_13TeV_DATA2017_AOD import dataSamples_17Nov2017_AOD
    selectedComponents = [ d for d in dataSamples_17Nov2017_AOD 
                           if "DoubleMuon" in d.name or "ZeroBias" in d.name ]

if run == "sig" or run == "mc" : isoTrackDeDxAna.doDeDx = True

if run == "data":
    for c in  selectedComponents:
       c.splitFactor = len(c.files)/3

if run == "mc":
    for c in  selectedComponents:
       c.splitFactor = len(c.files)/2

if run == "aod":
    prescaleComponents(selectedComponents, int(getHeppyOption("prescale","1")))
    configureSplittingFromTime(selectedComponents, 50, 2.5, maxFiles=4)

#-------- SEQUENCE -----------
sequence = cfg.Sequence( xtracks_sequence )
if run == "aod":
    sequence = cfg.Sequence( xtracks_sequence_AOD )


#-------- HOW TO RUN -----------
test = getHeppyOption('test')
if test == "1":
    if getHeppyOption("sample"):
        selected = [ s for s in selectedComponents if getHeppyOption("sample") == s.name ]
        if not selected:
            print "Sample %s not found. Known samples are: %s" % [ getHeppyOption("sample"), ", ".join(sorted(s.name for s in selectedComponents)) ]
            sys.exit(1)
        selectedComponents = selected
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
elif test == "1A":
    comp = dataSamples[0]
    comp.name = "AOD"
    comp.files = [ '/eos/cms/store/cmst3/user/gpetrucc/Run2017D_ZeroBias_17Nov2017_AOD.root' ]
    comp.triggers = []
    sequence = cfg.Sequence( xtracks_sequence_AOD )
    selectedComponents = doTest1(comp, sequence=sequence, cache=False )
    print "The test wil use file %s " % comp.files[0]
elif test in ('2','3','5s'):
    doTestN(test,selectedComponents)

printSummary(selectedComponents)
if getHeppyOption("justSummary",False): sys.exit(0)
config = autoConfig(selectedComponents, sequence)
