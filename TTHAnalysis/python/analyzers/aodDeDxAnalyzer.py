import os 

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle

from CMGTools.RootTools.utils.trackerTopologyHelper import loadTrackerTopology

class aodDeDxAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(aodDeDxAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)
        self.topology = loadTrackerTopology(os.path.expandvars(cfg_ana.trackerTopology[0]), cfg_ana.trackerTopology[1])

    def declareHandles(self):
        super(aodDeDxAnalyzer, self).declareHandles()
        self.handles['tracks'] = AutoHandle("generalTracks","std::vector<reco::Track>")            
        self.handles['dedxA']  = AutoHandle("dedxHitInfo","edm::Association<reco::DeDxHitInfoCollection>")            
        self.handles['dedxP'] = AutoHandle("dedxHitInfo:prescale","edm::ValueMap<int>")            
        self.handles['strip'] = AutoHandle("dedxHarmonic2", "edm::ValueMap<reco::DeDxData>")

    def beginLoop(self,setup):
        super(aodDeDxAnalyzer,self).beginLoop(setup)
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('all tracks')
        count.register('all pixel tracks')
        count.register('accepted tracks')


    def process(self, event):
        self.readCollections( event.input )

        getref = self.handles['dedxA'].product().get
        prescale = self.handles['dedxP'].product()
        getstrip = self.handles['strip'].product().get

        pixelChargeToEnergyCoefficient = 3.61e-6
 
        tracks = []
        for itk, tk in enumerate(self.handles['tracks'].product()):
            self.counters.counter('events').inc('all tracks')
            if tk.hitPattern().numberOfValidPixelHits() < 3: continue
            self.counters.counter('events').inc('all pixel tracks')
            if tk.pt() < 0.5 or tk.pt() >= 10: continue
            ref = getref(itk)
            if ref.isNull(): continue

            tk.dedxByLayer          = [0 for i in xrange(20)]
            tk.subDetIdByLayer      = [0 for i in xrange(20)]
            tk.layerOrSideByLayer   = [0 for i in xrange(20)]
            tk.ladderOrBladeByLayer = [0 for i in xrange(20)]
            tk.moduleByLayer        = [0 for i in xrange(20)]
            tk.sizeXbyLayer         = [0 for i in xrange(20)]
            tk.sizeYbyLayer         = [0 for i in xrange(20)]
 
            tracks.append(tk)
            tk.dedx = ref.get()
            tk.dxy = tk.dxy(event.goodVertices[0].position())
            tk.dz  = tk.dz(event.goodVertices[0].position())
            tk.dedxPrescale = prescale.get(ref.id(), ref.key())
            tk.stripDeDx = getstrip(itk)

            ifound = -1
            for ih in xrange(tk.dedx.size()):
                pixelCluster = tk.dedx.pixelCluster(ih)
                if pixelCluster:
                    ifound += 1
                    tk.dedxByLayer[ifound]     = pixelCluster.charge()/tk.dedx.pathlength(ih) * pixelChargeToEnergyCoefficient
                    tk.sizeXbyLayer[ifound]    = pixelCluster.sizeX()
                    tk.sizeYbyLayer[ifound]    = pixelCluster.sizeY()
                    detid = tk.dedx.detId(ih)
                    tk.subDetIdByLayer[ifound] = detid.subdetId()
                    tk.layerOrSideByLayer[ifound] = self.topology.layer(detid) 
                    if detid.subdetId() == 1: 
                        tk.ladderOrBladeByLayer[ifound] *= self.topology.pxbLadder(detid)
                    if detid.subdetId() == 2: 
                        tk.layerOrSideByLayer[ifound] *= 2*self.topology.side(detid)-3 # side is 2 for eta > 0, 1 for eta < 0 -> map to +1, -1
                        tk.ladderOrBladeByLayer[ifound] *= self.topology.pxfBlade(detid)
                    tk.moduleByLayer[ifound] = self.topology.module(detid) 
                 
            self.counters.counter('events').inc('accepted tracks')

        if not tracks:
            return False

        event.tracks = tracks
        return True
