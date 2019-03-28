# COMPONENT CREATOR
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

# QCD
QCD_HT100to200    = kreator.makeMCComponent("QCD_HT100to200"  , "/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"      , "CMS", ".*root", 2.463e+07*1.13073)
QCD_HT200to300    = kreator.makeMCComponent("QCD_HT200to300"  , "/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"           , "CMS", ".*root", 1.553e+06*1.1056 )
QCD_HT300to500    = kreator.makeMCComponent("QCD_HT300to500"  , "/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"           , "CMS", ".*root", 347500*1.01094   )
QCD_HT500to700    = kreator.makeMCComponent("QCD_HT500to700"  , "/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_old_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"   , "CMS", ".*root", 29930*1.0568     )
QCD_HT700to1000   = kreator.makeMCComponent("QCD_HT700to1000" , "/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"          , "CMS", ".*root", 6370*1.06782     )
QCD_HT1000to1500  = kreator.makeMCComponent("QCD_HT1000to1500", "/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"         , "CMS", ".*root", 1100*1.09636     )
QCD_HT1500to2000  = kreator.makeMCComponent("QCD_HT1500to2000", "/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"         , "CMS", ".*root", 98.71            )
QCD_HT2000toInf   = kreator.makeMCComponent("QCD_HT2000toInf" , "/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_old_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"  , "CMS", ".*root", 20.2             )

QCD = [
    QCD_HT100to200,
    QCD_HT200to300,
    QCD_HT300to500,
    QCD_HT500to700,
    QCD_HT700to1000,
    QCD_HT1000to1500,
    QCD_HT1500to2000,
    QCD_HT2000toInf,
]

# W + Jets
WJets_HT100to200   = kreator.makeMCComponent("WJets_HT100to200"  , "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"  , "CMS", ".*root", 1345.0)
WJets_HT200to400   = kreator.makeMCComponent("WJets_HT200to400"  , "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  , "CMS", ".*root",  359.7)
WJets_HT400to600   = kreator.makeMCComponent("WJets_HT400to600"  , "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  , "CMS", ".*root",  48.91)
WJets_HT600to800   = kreator.makeMCComponent("WJets_HT600to800"  , "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  , "CMS", ".*root",  12.05)
WJets_HT800to1200  = kreator.makeMCComponent("WJets_HT800to1200" , "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" , "CMS", ".*root",  5.501)
WJets_HT1200to2500 = kreator.makeMCComponent("WJets_HT1200to2500", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "CMS", ".*root",  1.329)
WJets_HT2500toInf  = kreator.makeMCComponent("WJets_HT2500toInf" , "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM" , "CMS", ".*root",  3.216e-2)

Ws = [ 
    WJets_HT100to200,
    WJets_HT200to400,
    WJets_HT400to600,
    WJets_HT600to800,
    WJets_HT800to1200,
    WJets_HT1200to2500,
    #WJets_HT2500toInf,
]

# DY + Jets
DYJetsM50               = kreator.makeMCComponent("DYJetsM50"              , "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                , "CMS", ".*root", 1921.8*3, fracNegWeights=0.16)
DYJetsM50e              = kreator.makeMCComponent("DYJetsM50e"             , "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"           , "CMS", ".*root", 1921.8*3, fracNegWeights=0.16)

DYJetsM50_HT100to200    = kreator.makeMCComponent("DYJetsM50_HT100to200"   , "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"     ,  "CMS", ".*root", 147.40)
DYJetsM50_HT100to200e   = kreator.makeMCComponent("DYJetsM50_HT100to200e"  , "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",  "CMS", ".*root", 147.40)
DYJetsM50_HT200to400    = kreator.makeMCComponent("DYJetsM50_HT200to400"   , "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"     ,  "CMS", ".*root",  40.99)
DYJetsM50_HT200to400e   = kreator.makeMCComponent("DYJetsM50_HT200to400e"  , "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",  "CMS", ".*root",  40.99)
DYJetsM50_HT400to600    = kreator.makeMCComponent("DYJetsM50_HT400to600"   , "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"     ,  "CMS", ".*root",  5.678)
DYJetsM50_HT400to600e   = kreator.makeMCComponent("DYJetsM50_HT400to600e"  , "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",  "CMS", ".*root",  5.678)
DYJetsM50_HT600to800    = kreator.makeMCComponent("DYJetsM50_HT600to800"   , "/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"     ,  "CMS", ".*root",  1.367)
DYJetsM50_HT800to1200   = kreator.makeMCComponent("DYJetsM50_HT800to1200"  , "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    ,  "CMS", ".*root", 0.6304)
DYJetsM50_HT1200to2500  = kreator.makeMCComponent("DYJetsM50_HT1200to2500" , "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"   ,  "CMS", ".*root", 0.1514)
DYJetsM50_HT2500toInf   = kreator.makeMCComponent("DYJetsM50_HT2500toInf"  , "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    ,  "CMS", ".*root", 3.565e-3)

DY = [
    DYJetsM50,                DYJetsM50e,  
    DYJetsM50_HT100to200,     DYJetsM50_HT100to200e,
    DYJetsM50_HT200to400,     DYJetsM50_HT200to400e,
    DYJetsM50_HT400to600,     DYJetsM50_HT400to600e,
    DYJetsM50_HT600to800,
    DYJetsM50_HT800to1200,
    DYJetsM50_HT1200to2500,
    DYJetsM50_HT2500toInf,
]

# Z(vv) + Jets
ZvvJets_HT100to200      = kreator.makeMCComponent("ZvvJets_HT100to200"  ,    "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  ,  "CMS", ".*root", 280.35  )
ZvvJets_HT200to400      = kreator.makeMCComponent("ZvvJets_HT200to400"  ,    "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  ,  "CMS", ".*root", 77.67   )
ZvvJets_HT400to600      = kreator.makeMCComponent("ZvvJets_HT400to600"  ,    "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  ,  "CMS", ".*root", 10.73   )
ZvvJets_HT600to800      = kreator.makeMCComponent("ZvvJets_HT600to800"  ,    "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  ,  "CMS", ".*root", 2.559   )
ZvvJets_HT800to1200     = kreator.makeMCComponent("ZvvJets_HT800to1200" ,    "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" ,  "CMS", ".*root", 1.1796  )
ZvvJets_HT1200to2500    = kreator.makeMCComponent("ZvvJets_HT1200to2500",    "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",  "CMS", ".*root", 0.28833 )
ZvvJets_HT2500toInf     = kreator.makeMCComponent("ZvvJets_HT2500toInf" ,    "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" ,  "CMS", ".*root", 0.006945)

Zvv = [
    ZvvJets_HT100to200,
    ZvvJets_HT200to400,
    ZvvJets_HT400to600,
    ZvvJets_HT600to800,
    ZvvJets_HT800to1200
    #ZvvJets_HT1200to2500,
    #ZvvJets_HT2500toInf
]

# TT
TTJets = kreator.makeMCComponent("TTJets", "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"              , "CMS", ".*root", 831.76, fracNegWeights=0.319   )
TTLep  = kreator.makeMCComponent("TTLep" , "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"       , "CMS", ".*root", 831.76*((3*0.108)**2)          )
TTHad  = kreator.makeMCComponent("TTHad" , "/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    , "CMS", ".*root", 831.76*((1-3*0.108)**2)        )
TTSemi = kreator.makeMCComponent("TTSemi", "/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "CMS", ".*root", 831.76*2*(3*0.108)*(1-3*0.108) )

TT = [ TTJets, TTLep, TTHad, TTSemi]

# Single top
T_tch     = kreator.makeMCComponent("T_tch"    , "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    , "CMS", ".*root", 136.02)
TBar_tch  = kreator.makeMCComponent("TBar_tch" , "/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "CMS", ".*root",  80.95) 
T_tWch    = kreator.makeMCComponent("T_tWch"   , "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"               , "CMS", ".*root",  19.55)
TBar_tWch = kreator.makeMCComponent("TBar_tWch", "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"      , "CMS", ".*root",  19.55)

Ts = [
    T_tch, 
    TBar_tch,
    T_tWch, 
    TBar_tWch
]

# Diboson
WW         = kreator.makeMCComponent("WW"        , "/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                                     , "CMS", ".*root", 63.21 * 1.82)
WZ         = kreator.makeMCComponent("WZ"        , "/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                                     , "CMS", ".*root", 47.13)
ZZ         = kreator.makeMCComponent("ZZ"        , "/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                                     , "CMS", ".*root", 16.523)

WWTo2L2Nu  = kreator.makeMCComponent("WWTo2L2Nu" , "/WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM", "CMS", ".*root", 10.481 )
WWToLNuQQ  = kreator.makeMCComponent("WWToLNuQQ" , "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"               , "CMS", ".*root", 43.53  )
WWToLNuQQe = kreator.makeMCComponent("WWToLNuQQe", "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"          , "CMS", ".*root", 43.53  )

WZTo3LNu   = kreator.makeMCComponent("WZTo3LNu"  , "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                  , "CMS", ".*root", 5.063, fracNegWeights=0.189 )
WZToLNu2Q  = kreator.makeMCComponent("WZToLNu2Q" , "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"               , "CMS", ".*root", 10.71, fracNegWeights=0.204 )
WZTo2L2Q   = kreator.makeMCComponent("WZTo2L2Q"  , "/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                  , "CMS", ".*root", 5.595, fracNegWeights=0.204 )

ZZTo4L     = kreator.makeMCComponent("ZZTo4L"    , "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                                  , "CMS", ".*root", 1.256)
ZZTo4Le    = kreator.makeMCComponent("ZZTo4Le"   , "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"                             , "CMS", ".*root", 1.256)
ZZTo2L2Q   = kreator.makeMCComponent("ZZTo2L2Q"  , "/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"                  , "CMS", ".*root", 3.220)


DiBosons = [
    WW,
    WWTo2L2Nu,
    WWToLNuQQ,
    WWToLNuQQe,
    WZ,
    WZTo3LNu,
    WZToLNu2Q,
    WZTo2L2Q,
    ZZ,
    ZZTo4L, 
    ZZTo4Le,
    ZZTo2L2Q,
]

# ----------------------------- summary ----------------------------------------


mcSamples = QCD + Ws + DY + Zvv + TT + Ts + DiBosons


samples = mcSamples

# ---------------------------------------------------------------------

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples, localobjs=locals())
