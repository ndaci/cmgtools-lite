#include <algorithm>
#include <iostream>
#include <TH1.h>
#include <TFile.h>
// Histograms for k-factors and SFs
TH1 *_zkfact = nullptr, *_wkfact = nullptr;

float weight_ZvvFromZll(float metNoMu_pt) { 
    // placeholder
    return 9; 
}

void _loadKFactors() {
    TFile *kfactfile = TFile::Open("../../macros/xtracks/KFactors/kfactors.root");
    TH1 *zknum = (TH1*) kfactfile->Get("EWKcorr/Z");
    TH1 *zkden = (TH1*) kfactfile->Get("ZJets_LO/inv_pt");
    TH1 *wknum = (TH1*) kfactfile->Get("EWKcorr/W");
    TH1 *wkden = (TH1*) kfactfile->Get("WJets_LO/inv_pt");
    _zkfact = (TH1*) zknum->Clone("_zkfact"); _zkfact->SetDirectory(nullptr); 
    _wkfact = (TH1*) wknum->Clone("_wkfact"); _wkfact->SetDirectory(nullptr); 
    _zkfact->Divide(zkden);
    _wkfact->Divide(wkden);
}

float zkfactor(float vpt) {
    if (!_zkfact) _loadKFactors();
    if (vpt <= 0) {
        static int _errors = 0;
        if (++_errors < 5) {
            std::cerr << "ERROR: requesting z k-factor at invalid vpt " << vpt << std::endl;
        }
        return 1;
    }
    vpt = std::min(std::max(150.1f, vpt), 1199.9f); // crop range
    return _zkfact->GetBinContent(_zkfact->FindBin(vpt));
}

float wkfactor(float vpt) {
    if (!_wkfact) _loadKFactors();
    if (vpt <= 0) {
        static int _errors = 0;
        if (++_errors < 5) {
            std::cerr << "ERROR: requesting w k-factor at invalid vpt " << vpt << std::endl;
        }
        return 1;
    }
    vpt = std::min(std::max(150.1f, vpt), 1199.9f); // crop range
    return _wkfact->GetBinContent(_wkfact->FindBin(vpt));
}


void functionsXtracks() {
}
