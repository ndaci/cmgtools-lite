import glob
import pickle
import sys
from ROOT import *
from array import array


args = sys.argv[1:]
if len(args) != 2:
    print "Need two command line arguments -- sample folder and sample string"
    quit()

sample = args[0]
path   = args[1]

sumwgt = 0.0

pckfile  = path+sample+"/skimAnalyzerCount/SkimReport.pck"
pckobj   = pickle.load(open(pckfile,'r'))
counters = dict(pckobj)
sumwgt  += counters['Sum Weights']

wgtsum = array( 'f', [ 0 ] )
wgtsum[0] = sumwgt
print wgtsum[0]

treefile = TFile(path+sample+"/treeProducerXtracks/tree.root", "update")
tree     = treefile.Get("tree")
entries  = tree.GetEntries()
tree.Branch('wgtsum', wgtsum, 'wgtsum/F')

for i in xrange(entries):
    tree.Fill()

treefile.Write()
treefile.Close()

