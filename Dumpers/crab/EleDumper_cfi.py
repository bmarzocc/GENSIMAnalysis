import FWCore.ParameterSet.Config as cms

elephodumper = cms.EDAnalyzer("ElePhoDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll",""),
    pileupSummary                   = cms.InputTag("slimmedAddPileupInfo",""),
    vertexCollection                = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    genParticleCollection           = cms.InputTag("prunedGenParticles",""),
    ebRechitCollection              = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeRechitCollection              = cms.InputTag("reducedEgamma","reducedEERecHits"),
    patElectronCollection           = cms.InputTag("slimmedElectrons",""), 
    patPhotonCollection             = cms.InputTag("slimmedPhotons",""), 
    
    isMC                            = cms.bool(True),  #isMC
    isDY                            = cms.bool(True),  #isDY
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)
    
    savePhotons                     = cms.bool(False),  #save patPhotons and patMET information
    saveElectrons                   = cms.bool(True),  #save patElectrons and patMET information

    egmCutBasedElectronIDVeto       = cms.string('cutBasedElectronID-Fall17-94X-V2-veto'),  #cutBasedEleID veto
    egmCutBasedElectronIDloose      = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'),  #cutBasedEleID loose  
    egmCutBasedElectronIDmedium     = cms.string('cutBasedElectronID-Fall17-94X-V2-medium'),  #cutBasedEleID medium 
    egmCutBasedElectronIDtight      = cms.string('cutBasedElectronID-Fall17-94X-V2-tight'),  #cutBasedEleID tight
    egmMVAElectronIDloose           = cms.string('mvaEleID-Fall17-iso-V2-wpLoose'),  #mvaEleID loose 
    egmMVAElectronIDmedium          = cms.string('mvaEleID-Fall17-iso-V2-wp90'),  #mvaEleID medium 
    egmMVAElectronIDtight           = cms.string('mvaEleID-Fall17-iso-V2-wp80'),  #mvaEleID tight  
    egmMVAElectronIDlooseNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wpLoose'),  #mvaEleIDNoIso loose 
    egmMVAElectronIDmediumNoIso     = cms.string('mvaEleID-Fall17-noIso-V2-wp90'),  #mvaEleIDNoIso medium 
    egmMVAElectronIDtightNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wp80'),  #mvaEleIDNoIso tight  
    heepElectronID                  = cms.string('heepElectronID-HEEPV70'),  #mvaEleIDNoIso tight
    egmCutBasedPhotonIDloose        = cms.string('cutBasedPhotonID-Fall17-94X-V2-loose'),  #cutBasedPhoID loose  
    egmCutBasedPhotonIDmedium       = cms.string('cutBasedPhotonID-Fall17-94X-V2-medium'),  #cutBasedPhoID medium 
    egmCutBasedPhotonIDtight        = cms.string('cutBasedPhotonID-Fall17-94X-V2-tight'),  #cutBasedPhoID tight 
    egmMVAPhotonIDmedium            = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),  #mvaPhoID medium 
    egmMVAPhotonIDtight             = cms.string('mvaPhoID-RunIIFall17-v2-wp80'),  #mvaPhoID tight     
)
