import FWCore.ParameterSet.Config as cms

simanalysis = cms.EDAnalyzer("simAnalysis",

    genParticles = cms.InputTag("genParticles",""),
    pcaloHitsEB  = cms.InputTag("g4SimHits", "EcalHitsEB", "SIM"), 
    pcaloHitsEE  = cms.InputTag("g4SimHits", "EcalHitsEE", "SIM")  
)
