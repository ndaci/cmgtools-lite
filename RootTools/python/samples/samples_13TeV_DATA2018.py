# COMPONENT CREATOR
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

# ----------------------------- 2018 pp run  ----------------------------------------
#
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2018Analysis
#

json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

run_range = (315252, 325175)
label = "_runs%s_%s" % (run_range[0], run_range[1])

# ----------------------------- Run2018A ----------------------------------------

JetHT_Run2018A = kreator.makeDataComponent("JetHT_Run2018A", "/JetHT/Run2018A-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
MET_Run2018A = kreator.makeDataComponent("MET_Run2018A", "/MET/Run2018A-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2018A = kreator.makeDataComponent("SingleMuon_Run2018A", "/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD", "CMS", ".*root", json)
EGamma_Run2018A = kreator.makeDataComponent("EGamma_Run2018A", "/EGamma/Run2018A-17Sep2018-v2/MINIAOD", "CMS", ".*root", json)
MuonEG_Run2018A = kreator.makeDataComponent("MuonEG_Run2018A", "/MuonEG/Run2018A-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
DoubleMuon_Run2018A = kreator.makeDataComponent("DoubleMuon_Run2018A", "/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD", "CMS", ".*root", json)

dataSamples_Run2018A = [JetHT_Run2018A, MET_Run2018A, SingleMuon_Run2018A, EGamma_Run2018A, MuonEG_Run2018A, DoubleMuon_Run2018A]


# ----------------------------- Run2018B ----------------------------------------

JetHT_Run2018B = kreator.makeDataComponent("JetHT_Run2018B", "/JetHT/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
MET_Run2018B = kreator.makeDataComponent("MET_Run2018B", "/MET/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2018B = kreator.makeDataComponent("SingleMuon_Run2018B", "/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
EGamma_Run2018B = kreator.makeDataComponent("EGamma_Run2018B", "/EGamma/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
MuonEG_Run2018B = kreator.makeDataComponent("MuonEG_Run2018B", "/MuonEG/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
DoubleMuon_Run2018B = kreator.makeDataComponent("DoubleMuon_Run2018B", "/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)

dataSamples_Run2018B = [JetHT_Run2018B, MET_Run2018B, SingleMuon_Run2018B, EGamma_Run2018B, MuonEG_Run2018B, DoubleMuon_Run2018B]


# ----------------------------- Run2018C ----------------------------------------

JetHT_Run2018C = kreator.makeDataComponent("JetHT_Run2018C", "/JetHT/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
MET_Run2018C = kreator.makeDataComponent("MET_Run2018C", "/MET/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2018C = kreator.makeDataComponent("SingleMuon_Run2018C", "/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
EGamma_Run2018C = kreator.makeDataComponent("EGamma_Run2018C", "/EGamma/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
MuonEG_Run2018C = kreator.makeDataComponent("MuonEG_Run2018C", "/MuonEG/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)
DoubleMuon_Run2018C = kreator.makeDataComponent("DoubleMuon_Run2018C", "/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD", "CMS", ".*root", json)

dataSamples_Run2018C = [JetHT_Run2018C, MET_Run2018C, SingleMuon_Run2018C, EGamma_Run2018C, MuonEG_Run2018C, DoubleMuon_Run2018C]


# ----------------------------- Run2018D ----------------------------------------

JetHT_Run2018D = kreator.makeDataComponent("JetHT_Run2018D", "/JetHT/Run2018D-PromptReco-v1/MINIAOD", "CMS", ".*root", json)
MET_Run2018D = kreator.makeDataComponent("MET_Run2018D", "/MET/Run2018D-PromptReco-v1/MINIAOD", "CMS", ".*root", json)
SingleMuon_Run2018D = kreator.makeDataComponent("SingleMuon_Run2018D", "/SingleMuon/Run2018D-PromptReco-v2/MINIAOD", "CMS", ".*root", json)
EGamma_Run2018D = kreator.makeDataComponent("EGamma_Run2018D", "/EGamma/Run2018D-PromptReco-v2/MINIAOD", "CMS", ".*root", json)
MuonEG_Run2018D = kreator.makeDataComponent("MuonEG_Run2018D", "/MuonEG/Run2018D-PromptReco-v2/MINIAOD", "CMS", ".*root", json)
DoubleMuon_Run2018D = kreator.makeDataComponent("DoubleMuon_Run2018D", "/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD", "CMS", ".*root", json)

dataSamples_Run2018D = [JetHT_Run2018D, MET_Run2018D, SingleMuon_Run2018D, EGamma_Run2018D, MuonEG_Run2018D, DoubleMuon_Run2018D]



# ------------------------------------------------------------------------------------
# Summary of 1Apr2019
dataSamples_1Apr2019 = dataSamples_Run2018A + dataSamples_Run2018B + dataSamples_Run2018C + dataSamples_Run2018D 



dataSamples = dataSamples_1Apr2019
samples = dataSamples

# ---------------------------------------------------------------------

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples)

