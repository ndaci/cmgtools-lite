#include <algorithm>
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
// Histograms for k-factors and SFs
TH1 *_zkfact = nullptr, *_wkfact = nullptr;
// Histograms for muon SFs
TH2 *_muSF_LooseId =nullptr, *_muSF_TightId = nullptr;

// Utility function
float _readHisto2D(float x, float y, TH2 *h2, float nsigma = 0, bool warnIfZero=false) {
    int ix = std::min(std::max(1, h2->GetXaxis()->FindBin(x)), h2->GetNbinsX());
    int iy = std::min(std::max(1, h2->GetYaxis()->FindBin(y)), h2->GetNbinsY());
    float val = h2->GetBinContent(ix,iy);
    if (nsigma != 0) val += nsigma * h2->GetBinError(ix,iy);
    if (warnIfZero && val <= 0) {
        static int _errors = 0;
        if (++_errors < 5) {
            std::cerr << "ERROR: got a bad value " << val << " from " << h2->GetName() << " at x = " << x << ", y = " << y << ", nsigma = " << nsigma << std::endl;
        }
    }
    return val;
}

// ===== Code for CR->SR extrapolation etc. ======

float weight_ZvvFromZll(float metNoMu_pt) { 
    // placeholder
    return 9; 
}

// ====== Code for k-factors =====

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

/// ========== Code for scale factors =========

void _loadMuSF() {
   TFile *idSF = TFile::Open("../../macros/xtracks/MuonIDIsoTrig/Muon_ID_SF_2017.root");
   _muSF_TightId = (TH2*) idSF->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone("_muSF_TightId");
   _muSF_LooseId = (TH2*) idSF->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone("_muSF_LooseId");
}
float muonSF(float pt, float eta, int tightId, float nsigma=0) {
    if (!_muSF_TightId) _loadMuSF();
    // if it passes the tight we apply the tight id SF
    // otherwise we apply the loose id SF (this is not really correct, but it's just an example)
    if (tightId) return _readHisto2D(pt, std::abs(eta), _muSF_TightId, nsigma, /*warnIfZero=*/true); 
    else         return _readHisto2D(pt, std::abs(eta), _muSF_LooseId, nsigma, /*warnIfZero=*/true); 
}

// =========
// leave this dummy at the end
void functionsXtracks() {
}
