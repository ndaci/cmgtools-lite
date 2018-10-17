# COMPONENT CREATOR
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

# ----------------------------- 2017 pp run  ----------------------------------------

json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'

datasets = {}
datasetsByEra = {}
datasetsByPD  = {}
for PD in  "SingleMuon", "SingleElectron", "SinglePhoton", "DoubleMuon", "DoubleEG", "MuonEG", "Tau", "MET", "JetHT", "HTMHT", "ZeroBias", "MuOnia":
    if PD not in datasetsByPD: datasetsByPD[PD] = []
    for Era in ["Run2017"+X for X in "BCDEF"]:
        if Era not in datasetsByEra: datasetsByEra[Era] = []
        dataset = kreator.makeDataComponent("%s_%s_17Nov2017_AOD"%(PD,Era), "/%s/%s-17Nov2017-v1/AOD" % (PD,Era), "CMS", ".*root", json)
        datasets[dataset.name] = dataset
        datasetsByEra[Era].append(dataset)
        datasetsByPD[PD].append(dataset)
        locals()[dataset.name] = dataset
del locals()['dataset'] # avoids warning in checkdecl

dataSamples_17Nov2017_AOD = datasets.values()
samples = dataSamples_17Nov2017_AOD
# ---------------------------------------------------------------------

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples,localobjs=locals())
