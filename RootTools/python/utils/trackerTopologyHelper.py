import ROOT

def loadTrackerTopology(filename, toponame):
    ROOT.gSystem.Load("libCondFormatsGeometryObjects.so");
    fIn = ROOT.TFile.Open(filename)
    if not fIn: raise RuntimeError("Couldn't open %s" % filename)
    ptp = fIn.Get(toponame)
    if not ptp: raise RuntimeError("Couldn't get %s from %s" % (toponame, filename))

    ROOT.gSystem.Load("libDataFormatsTrackerCommon.so");
    ROOT.gInterpreter.ProcessLine("#include <DataFormats/TrackerCommon/interface/TrackerTopology.h>")

    pxbVals_ = ROOT.TrackerTopology.PixelBarrelValues()
    pxfVals_ = ROOT.TrackerTopology.PixelEndcapValues()
    tecVals_ = ROOT.TrackerTopology.TECValues()
    tibVals_ = ROOT.TrackerTopology.TIBValues()
    tidVals_ = ROOT.TrackerTopology.TIDValues()
    tobVals_ = ROOT.TrackerTopology.TOBValues()

    #from https:##github.com/cms-sw/cmssw/blob/CMSSW_9_4_7/Geometry/TrackerNumberingBuilder/plugins/TrackerTopologyEP.cc
    pxbVals_.layerStartBit_ = ptp.vitems[0].vpars[0]; ## 16
    pxbVals_.ladderStartBit_ = ptp.vitems[0].vpars[1]; ## 8
    pxbVals_.moduleStartBit_ = ptp.vitems[0].vpars[2]; ## 2
    pxbVals_.layerMask_ = ptp.vitems[0].vpars[3]; ## 0xF
    pxbVals_.ladderMask_ = ptp.vitems[0].vpars[4]; ## 0xFF
    pxbVals_.moduleMask_ = ptp.vitems[0].vpars[5]; ## 0x3F
    
    pxfVals_.sideStartBit_ = ptp.vitems[1].vpars[0];
    pxfVals_.diskStartBit_ = ptp.vitems[1].vpars[1];
    pxfVals_.bladeStartBit_ = ptp.vitems[1].vpars[2];
    pxfVals_.panelStartBit_ = ptp.vitems[1].vpars[3];
    pxfVals_.moduleStartBit_ = ptp.vitems[1].vpars[4];
    pxfVals_.sideMask_ = ptp.vitems[1].vpars[5];
    pxfVals_.diskMask_ = ptp.vitems[1].vpars[6];
    pxfVals_.bladeMask_ = ptp.vitems[1].vpars[7];
    pxfVals_.panelMask_ = ptp.vitems[1].vpars[8];
    pxfVals_.moduleMask_ = ptp.vitems[1].vpars[9];
    
    ## TEC: 6
    tecVals_.sideStartBit_ = ptp.vitems[5].vpars[0];
    tecVals_.wheelStartBit_ = ptp.vitems[5].vpars[1];
    tecVals_.petal_fw_bwStartBit_ = ptp.vitems[5].vpars[2];
    tecVals_.petalStartBit_ = ptp.vitems[5].vpars[3];
    tecVals_.ringStartBit_ = ptp.vitems[5].vpars[4];
    tecVals_.moduleStartBit_ = ptp.vitems[5].vpars[5];
    tecVals_.sterStartBit_ = ptp.vitems[5].vpars[6];
    tecVals_.sideMask_ = ptp.vitems[5].vpars[7];
    tecVals_.wheelMask_ = ptp.vitems[5].vpars[8];
    tecVals_.petal_fw_bwMask_ = ptp.vitems[5].vpars[9];
    tecVals_.petalMask_ = ptp.vitems[5].vpars[10];
    tecVals_.ringMask_ = ptp.vitems[5].vpars[11];
    tecVals_.moduleMask_ = ptp.vitems[5].vpars[12];
    tecVals_.sterMask_ = ptp.vitems[5].vpars[13];
    
    ## TIB: 3
    tibVals_.layerStartBit_ = ptp.vitems[2].vpars[0];
    tibVals_.str_fw_bwStartBit_ = ptp.vitems[2].vpars[1];
    tibVals_.str_int_extStartBit_ = ptp.vitems[2].vpars[2];
    tibVals_.strStartBit_ = ptp.vitems[2].vpars[3];
    tibVals_.moduleStartBit_ = ptp.vitems[2].vpars[4];
    tibVals_.sterStartBit_ = ptp.vitems[2].vpars[5];
    tibVals_.layerMask_ = ptp.vitems[2].vpars[6];
    tibVals_.str_fw_bwMask_ = ptp.vitems[2].vpars[7];
    tibVals_.str_int_extMask_ = ptp.vitems[2].vpars[8];
    tibVals_.strMask_ = ptp.vitems[2].vpars[9];
    tibVals_.moduleMask_ = ptp.vitems[2].vpars[10];
    tibVals_.sterMask_ = ptp.vitems[2].vpars[11];
    
    ## TID: 4
    tidVals_.sideStartBit_= ptp.vitems[3].vpars[0];
    tidVals_.wheelStartBit_= ptp.vitems[3].vpars[1];
    tidVals_.ringStartBit_= ptp.vitems[3].vpars[2];
    tidVals_.module_fw_bwStartBit_= ptp.vitems[3].vpars[3];
    tidVals_.moduleStartBit_= ptp.vitems[3].vpars[4];
    tidVals_.sterStartBit_= ptp.vitems[3].vpars[5];
    tidVals_.sideMask_= ptp.vitems[3].vpars[6];
    tidVals_.wheelMask_= ptp.vitems[3].vpars[7];
    tidVals_.ringMask_= ptp.vitems[3].vpars[8];
    tidVals_.module_fw_bwMask_= ptp.vitems[3].vpars[9];
    tidVals_.moduleMask_= ptp.vitems[3].vpars[10];
    tidVals_.sterMask_= ptp.vitems[3].vpars[11];
    
    ## TOB: 5
    tobVals_.layerStartBit_ = ptp.vitems[4].vpars[0];
    tobVals_.rod_fw_bwStartBit_= ptp.vitems[4].vpars[1];
    tobVals_.rodStartBit_= ptp.vitems[4].vpars[2];
    tobVals_.moduleStartBit_= ptp.vitems[4].vpars[3];
    tobVals_.sterStartBit_= ptp.vitems[4].vpars[4];
    tobVals_.layerMask_= ptp.vitems[4].vpars[5];
    tobVals_.rod_fw_bwMask_= ptp.vitems[4].vpars[6];
    tobVals_.rodMask_= ptp.vitems[4].vpars[7];
    tobVals_.moduleMask_= ptp.vitems[4].vpars[8];
    tobVals_.sterMask_= ptp.vitems[4].vpars[9];

    return ROOT.TrackerTopology( pxbVals_, pxfVals_, tecVals_, tibVals_, tidVals_, tobVals_ )

if __name__ == "__main__":
    import sys
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libFWCoreFWLite.so");
    ROOT.gSystem.Load("libDataFormatsFWLite.so");
    ROOT.FWLiteEnabler.enable()
    
    print sys.argv
    if len(sys.argv) >= 3:
        print "Reading %s from file %s" % (sys.argv[2], sys.argv[1])
    topo = loadTrackerTopology(sys.argv[1], sys.argv[2])
