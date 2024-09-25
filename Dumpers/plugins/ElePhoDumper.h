#ifndef RecoSimStudies_Dumpers_ElePhoDumper_H
#define RecoSimStudies_Dumpers_ElePhoDumper_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

class ElePhoDumper : public edm::EDAnalyzer
{
      public:
         explicit ElePhoDumper(const edm::ParameterSet&);
	 ~ElePhoDumper();
  
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      void setTree(TTree* tree);
      void setVectors(int patElectron_size, int patPhoton_size);
      int getGenStatusFlag(const reco::GenParticle* genParticle);
      void printGenStatusFlag(const reco::GenParticle* genParticle);
      int passPreselections(const pat::Photon* photon);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<pat::Electron> > patElectronToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > patPhotonToken_;
      
      edm::Service<TFileService> iFile;
      
      // ----------config inputs-------------------
      bool isMC_;
      bool isDY_;
      bool doCompression_;
      int nBits_;
      bool savePhotons_; 
      bool saveElectrons_; 
      
      std::string egmCutBasedElectronIDVeto_;
      std::string egmCutBasedElectronIDloose_;
      std::string egmCutBasedElectronIDmedium_;
      std::string egmCutBasedElectronIDtight_; 
      std::string egmMVAElectronIDloose_;
      std::string egmMVAElectronIDmedium_;
      std::string egmMVAElectronIDtight_; 
      std::string egmMVAElectronIDlooseNoIso_;
      std::string egmMVAElectronIDmediumNoIso_;
      std::string egmMVAElectronIDtightNoIso_; 
      std::string heepElectronID_; 
      std::string egmCutBasedPhotonIDloose_;
      std::string egmCutBasedPhotonIDmedium_;
      std::string egmCutBasedPhotonIDtight_;  
      std::string egmMVAPhotonIDmedium_;
      std::string egmMVAPhotonIDtight_; 
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      long int eventId;
      int lumiId;
      int runId; 
      double truePU;
      double obsPU;
      int nVtx;
      float rho; 
      int genParticle_size; 
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<int> genParticle_statusFlag; 
      std::vector<int> genParticle_statusFlag_isLastCopyBeforeFSR;            
      std::vector<int> genParticle_statusFlag_isLastCopy;  
      std::vector<int> genParticle_statusFlag_isFirstCopy;   
      std::vector<int> genParticle_statusFlag_fromHardProcessBeforeFSR;  
      std::vector<int> genParticle_statusFlag_isDirectHardProcessTauDecayProduct; 
      std::vector<int> genParticle_statusFlag_isHardProcessTauDecayProduct;   
      std::vector<int> genParticle_statusFlag_fromHardProcess;  
      std::vector<int> genParticle_statusFlag_isHardProcess;
      std::vector<int> genParticle_statusFlag_isDirectHadronDecayProduct;  
      std::vector<int> genParticle_statusFlag_isDirectPromptTauDecayProduct;   
      std::vector<int> genParticle_statusFlag_isDirectTauDecayProduct;   
      std::vector<int> genParticle_statusFlag_isPromptTauDecayProduct;   
      std::vector<int> genParticle_statusFlag_isTauDecayProduct;    
      std::vector<int> genParticle_statusFlag_isDecayedLeptonHadron;   
      std::vector<int> genParticle_statusFlag_isPrompt;
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;
      int patElectron_size;
      std::vector<int> patElectron_index; 
      std::vector<int> patElectron_classification;
      std::vector<int> patElectron_scNPFClusters;
      std::vector<int> patElectron_charge;
      std::vector<bool> patElectron_isEB;
      std::vector<bool> patElectron_isEE;
      std::vector<bool> patElectron_isEBEEGap;
      std::vector<bool> patElectron_isEBEtaGap;
      std::vector<bool> patElectron_isEBPhiGap;
      std::vector<bool> patElectron_isEEDeeGap;
      std::vector<bool> patElectron_isEERingGap;
      std::vector<bool> patElectron_isEcalDriven;
      std::vector<bool> patElectron_isTrackerDriven;
      std::vector<bool> patElectron_passConversionVeto;
      std::vector<float> patElectron_eta;
      std::vector<float> patElectron_phi;
      std::vector<float> patElectron_p;
      std::vector<float> patElectron_pt;
      std::vector<float> patElectron_pIn;
      std::vector<float> patElectron_pOut;
      std::vector<float> patElectron_trackFbrem;
      std::vector<float> patElectron_superClusterFbrem;
      std::vector<float> patElectron_energy;
      std::vector<float> patElectron_energyErr;
      std::vector<float> patElectron_ecalEnergy;
      std::vector<float> patElectron_ecalEnergyErr;
      std::vector<float> patElectron_et;
      std::vector<float> patElectron_scEnergy;
      std::vector<float> patElectron_scRawEnergy;
      std::vector<float> patElectron_scRawESEnergy;
      std::vector<float> patElectron_scEt;
      std::vector<float> patElectron_scPhiWidth;
      std::vector<float> patElectron_scEtaWidth;
      std::vector<float> patElectron_scEoP;
      std::vector<float> patElectron_scEta;
      std::vector<float> patElectron_scPhi;
      std::vector<float> patElectron_scSwissCross;
      std::vector<float> patElectron_scEMax;
      std::vector<float> patElectron_scR9; 
      std::vector<float> patElectron_scSigmaIEtaIEta;
      std::vector<float> patElectron_scSigmaIEtaIPhi;
      std::vector<float> patElectron_scSigmaIPhiIPhi;
      std::vector<float> patElectron_full5x5_scEMax;
      std::vector<float> patElectron_full5x5_scR9;  
      std::vector<float> patElectron_full5x5_scSigmaIEtaIEta;
      std::vector<float> patElectron_full5x5_scSigmaIEtaIPhi;
      std::vector<float> patElectron_full5x5_scSigmaIPhiIPhi;
      std::vector<float> patElectron_HoE;
      std::vector<int> patElectron_egmCutBasedElectronIDVeto;
      std::vector<int> patElectron_egmCutBasedElectronIDloose;
      std::vector<int> patElectron_egmCutBasedElectronIDmedium;
      std::vector<int> patElectron_egmCutBasedElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDloose;
      std::vector<int> patElectron_egmMVAElectronIDmedium;
      std::vector<int> patElectron_egmMVAElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDlooseNoIso;
      std::vector<int> patElectron_egmMVAElectronIDmediumNoIso;
      std::vector<int> patElectron_egmMVAElectronIDtightNoIso;
      std::vector<int> patElectron_heepElectronID;
      int patPhoton_size;
      std::vector<int> patPhoton_index; 
      std::vector<int> patPhoton_scNPFClusters;
      std::vector<bool> patPhoton_isEB;
      std::vector<bool> patPhoton_isEE;
      std::vector<bool> patPhoton_isEBEEGap;
      std::vector<bool> patPhoton_isEBEtaGap;
      std::vector<bool> patPhoton_isEBPhiGap;
      std::vector<bool> patPhoton_isEEDeeGap;
      std::vector<bool> patPhoton_isEERingGap;
      std::vector<bool> patPhoton_passElectronVeto;
      std::vector<bool> patPhoton_hasPixelSeed;
      std::vector<bool> patPhoton_hasConversionTracks;
      std::vector<int> patPhoton_nConversions;
      std::vector<int> patPhoton_nConversionsOneLeg;  
      std::vector<float> patPhoton_eta;
      std::vector<float> patPhoton_phi;
      std::vector<float> patPhoton_energy; 
      std::vector<float> patPhoton_energyErr;
      std::vector<float> patPhoton_ecalEnergy;
      std::vector<float> patPhoton_ecalEnergyErr;
      std::vector<float> patPhoton_et;
      std::vector<float> patPhoton_mt;
      std::vector<float> patPhoton_dphiMET;     
      std::vector<float> patPhoton_scEnergy;  
      std::vector<float> patPhoton_scRawEnergy;  
      std::vector<float> patPhoton_scRawESEnergy;
      std::vector<float> patPhoton_scEt;
      std::vector<float> patPhoton_scEtaWidth;
      std::vector<float> patPhoton_scPhiWidth;    
      std::vector<float> patPhoton_scEta;
      std::vector<float> patPhoton_scPhi;
      std::vector<float> patPhoton_scSwissCross;
      std::vector<float> patPhoton_scE2x2;
      std::vector<float> patPhoton_scE3x3; 
      std::vector<float> patPhoton_scE5x5; 
      std::vector<float> patPhoton_scEMax;
      std::vector<float> patPhoton_scR9;
      std::vector<float> patPhoton_scSigmaIEtaIEta;
      std::vector<float> patPhoton_scSigmaIEtaIPhi;
      std::vector<float> patPhoton_scSigmaIPhiIPhi;
      std::vector<float> patPhoton_full5x5_scEMax;
      std::vector<float> patPhoton_full5x5_scR9;
      std::vector<float> patPhoton_full5x5_scSigmaIEtaIEta;
      std::vector<float> patPhoton_full5x5_scSigmaIEtaIPhi;
      std::vector<float> patPhoton_full5x5_scSigmaIPhiIPhi;
      std::vector<float> patPhoton_HoE;
      std::vector<int> patPhoton_egmCutBasedPhotonIDloose;
      std::vector<int> patPhoton_egmCutBasedPhotonIDmedium;
      std::vector<int> patPhoton_egmCutBasedPhotonIDtight;
      std::vector<int> patPhoton_egmMVAPhotonIDmedium;
      std::vector<int> patPhoton_egmMVAPhotonIDtight;  
      std::vector<int> patPhoton_passPreselections;
  
};

#endif
