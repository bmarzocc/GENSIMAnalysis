# GENSIM-Analysis

Install:

    scram project CMSSW_10_6_20
    cd CMSSW_10_6_20/src/
    cmsenv
    git cms-init
    git cms-checkout-topic bmarzocc:10_6_20_ParticleGuns
    git clone git@github.com:bmarzocc/GENSIM-Analysis.git
    scram b -j 5

Run:

    cd GENSIM-Analysis/Dumpers/crab
    cmsRun Photons_closeEcal_E1-1000_simAnalysis_cfi.py

