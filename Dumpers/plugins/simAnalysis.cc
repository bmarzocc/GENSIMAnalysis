// -*- C++ -*-
//
// Package:    simAnalysis/simAnalysis
// Class:      simAnalysis
//
/**\class simAnalysis simAnalysis.cc simAnalysis/simAnalysis/plugins/simAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rajdeep Mohan Chatterjee
//         Created:  Fri, 16 Oct 2020 14:03:21 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "RecoCaloTools/Navigation/interface/CaloRectangle.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TVirtualFitter.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include <iostream>
#include <fstream>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class simAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit simAnalysis(const edm::ParameterSet&);
  ~simAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static const int offset = 0x3;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void setTree(TTree* tree);
  bool passPreselections(const pat::Photon* photon);
  std::vector<DetId> getRefinedSCIds(reco::SuperClusterRef sc);
  double getMean(std::vector<double>* vals, std::vector<double>* weights);
  double getError(std::vector<double>* vals, std::vector<double>* weights);
  void fillAvgProfile();
  // ----------member ---------------------------
  std::vector<TH1D*> Profiles;
  TH1D* maxProfile; 
  TH1D* sumProfile; 
  TH1D* avgProfile; 
  TF1* f;
  TTree* tree;

  long long int runId;
  long long int eventId;
  long long int lumiId;
  int genParticle_pdgId;
  int genParticle_status; 
  float genParticle_energy;
  float genParticle_pt;
  float genParticle_eta;
  float genParticle_phi; 
  float simEnergy;
  float simEnergy_withES; 
  float simR9;
  float simR9_withES;
  float recoPhoton_energy;
  float recoPhoton_scEnergy; 
  float recoPhoton_pt;
  float recoPhoton_eta;
  float recoPhoton_phi;
  float recoPhoton_scEta;
  float recoPhoton_scPhi;
  float recoPhoton_r9;
  float recoPhoton_full5x5_r9;
  int recoPhoton_nPFClusters;
  bool recoPhoton_hasConversionTracks;
  int recoPhoton_nConversions;  
  int recoPhoton_nConversionsOneLeg;  
  int recoPhoton_egmCutBasedPhotonIDloose;
  int recoPhoton_egmCutBasedPhotonIDmedium;
  int recoPhoton_egmCutBasedPhotonIDtight;
  int recoPhoton_egmMVAPhotonIDmedium;
  int recoPhoton_egmMVAPhotonIDtight; 
  int recoPhoton_passPreselections;
  float recoElectron_energy;
  float recoElectron_scEnergy;
  float recoElectron_pt;
  float recoElectron_eta;
  float recoElectron_phi;
  float recoElectron_scEta;
  float recoElectron_scPhi;
  float recoElectron_r9;
  float recoElectron_full5x5_r9;
  int recoElectron_nBrems;
  float recoElectron_superClusterFbrem;
  float recoElectron_trackFbrem; 
  bool recoElectron_passConversionVeto;
  int recoElectron_nConversions;  
  int recoElectron_nConversionsOneLeg;  
  int recoElectron_egmCutBasedElectronIDVeto;  
  int recoElectron_egmCutBasedElectronIDloose;  
  int recoElectron_egmCutBasedElectronIDmedium;  
  int recoElectron_egmCutBasedElectronIDtight;  
  int recoElectron_egmMVAElectronIDloose;  
  int recoElectron_egmMVAElectronIDmedium;  
  int recoElectron_egmMVAElectronIDtight;  
  int recoElectron_egmMVAElectronIDlooseNoIso;  
  int recoElectron_egmMVAElectronIDmediumNoIso;  
  int recoElectron_egmMVAElectronIDtightNoIso;  
  int recoElectron_heepElectronID;  
  std::vector<float> xtalMatrix_energy;
  std::vector<int> xtalMatrix_ieta;
  std::vector<int> xtalMatrix_iphi; 
  std::vector<int> xtalMatrix_iz; 
  float xtalMatrix_totEnergy;
  float xtalMatrix_e3x3;
  float xtalMax_binMaxVal;
  float xtalMax_tMax;
  float xtalMax_tMaxError;
  float xtalMax_tMax_fitPar0;
  float xtalMax_tMax_fitPar1;
  float xtalMax_tMax_fitPar2;
  float xtalMax_tMax_fitParError0;
  float xtalMax_tMax_fitParError1;
  float xtalMax_tMax_fitParError2; 
  int xtalMax_tMax_fitStatus;
  float xtalMax_energy;
  int xtalMax_ieta;
  int xtalMax_iphi;
  int xtalMax_iz;
  std::vector<float> xtalMax_showerProfile;
  std::vector<float> xtalMax_showerProfileError;
  std::vector<float> avgMatrix_showerProfile;
  std::vector<float> avgMatrix_showerProfileError;
  std::vector<float> sumMatrix_showerProfile;
  std::vector<float> sumMatrix_showerProfileError;

  std::map<long long int,float> energyMapEB;
  std::map<long long int,float> energyMapEE;
  std::map<long long int,float> energyMapES; 
  std::vector<DetId> xtalMatrix;
  std::vector<DetId> xtalMatrix_3x3;
  std::vector<DetId> xtalRefinedSC;

  const CaloTopology* topology;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* eeGeometry;
  edm::ESHandle<CaloTopology> caloTopology;

  edm::Handle<std::vector<reco::GenParticle> > genParticle_Handle;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticle_Token;

  edm::Handle<edm::PCaloHitContainer> pCaloHits_EB_Handle;
  edm::EDGetTokenT<edm::PCaloHitContainer> pCaloHits_EB_Token;
  
  edm::Handle<edm::PCaloHitContainer> pCaloHits_EE_Handle;
  edm::EDGetTokenT<edm::PCaloHitContainer> pCaloHits_EE_Token;

  edm::Handle<edm::PCaloHitContainer> pCaloHits_ES_Handle;
  edm::EDGetTokenT<edm::PCaloHitContainer> pCaloHits_ES_Token; 

  const std::vector<pat::Electron>* electrons;
  edm::Handle<std::vector<pat::Electron> > electron_Handle;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_Token_;

  const std::vector<pat::Photon>* photons;
  edm::Handle<std::vector<pat::Photon> > photon_Handle;
  edm::EDGetTokenT<std::vector<pat::Photon> > photon_Token_;

  TRandom3 *gRandom;
  double depthUnit_;
  std::vector<double> showerProfileBinning_; 
  bool saveShowerProfile_; 
  bool saveShowerProfileRefinedSC_;
  int precision_;
  bool runOnReco_;
  bool savePhotons_;
  int matrixSize_; 
  int dim_;
  std::vector<double>* vals;
  std::vector<double>* weights;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
simAnalysis::simAnalysis(const edm::ParameterSet& iConfig)
{
  usesResource("TFileService");
  genParticle_Token           = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
  pCaloHits_EB_Token          = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("pcaloHitsEB"));
  pCaloHits_EE_Token          = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("pcaloHitsEE"));
  pCaloHits_ES_Token          = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("pcaloHitsES"));
  electron_Token_             = consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));
  photon_Token_               = consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"));
  depthUnit_                  = iConfig.getParameter<double>("depthUnit");    
  showerProfileBinning_       = iConfig.getParameter<std::vector<double>>("showerProfileBinning"); 
  saveShowerProfile_          = iConfig.getParameter<bool>("saveShowerProfile"); 
  saveShowerProfileRefinedSC_ = iConfig.getParameter<bool>("saveShowerProfileRefinedSC");  
  precision_                  = iConfig.getParameter<int>("precision");  
  runOnReco_                  = iConfig.getParameter<bool>("runOnReco");  
  matrixSize_                 = iConfig.getParameter<int>("matrixSize");
  dim_ = matrixSize_ * matrixSize_;

  gRandom = new TRandom3(); 
  vals = new std::vector<double>;
  weights = new std::vector<double>;
  xtalMatrix.resize(dim_);
  xtalMatrix_3x3.resize(9);
  
  if(saveShowerProfileRefinedSC_) dim_ = 1000; 
  Profiles.resize(dim_);

  edm::Service<TFileService> fs;
  
  for(int iXtal=0; iXtal<dim_; iXtal++){
      Profiles[iXtal] = new TH1D(Form("Profiles_%d",iXtal),Form("Profiles_%d",iXtal),int(showerProfileBinning_.at(0)),showerProfileBinning_.at(1), showerProfileBinning_.at(2)); 
      Profiles[iXtal]->Sumw2(); 
  }
  maxProfile = new TH1D("maxProfile","maxProfile",int(showerProfileBinning_.at(0)),showerProfileBinning_.at(1), showerProfileBinning_.at(2));
  maxProfile->Sumw2(); 
  sumProfile = new TH1D("sumProfile","sumProfile",int(showerProfileBinning_.at(0)),showerProfileBinning_.at(1), showerProfileBinning_.at(2));
  sumProfile->Sumw2(); 
  avgProfile = new TH1D("avgProfile","avgProfile",int(showerProfileBinning_.at(0)),showerProfileBinning_.at(1), showerProfileBinning_.at(2));
  avgProfile->Sumw2(); 

  tree = fs->make<TTree>("caloTree","caloTree"); 
  setTree(tree);
        
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

simAnalysis::~simAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void simAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{

  iSetup.get<CaloTopologyRecord>().get(caloTopology);
  topology = caloTopology.product();

  iEvent.getByToken(genParticle_Token, genParticle_Handle); //GenParticles
  const std::vector<reco::GenParticle>* genParticles = genParticle_Handle.product();
  
  iEvent.getByToken(pCaloHits_EB_Token, pCaloHits_EB_Handle); //PCaloHits EB
  iEvent.getByToken(pCaloHits_EE_Token, pCaloHits_EE_Handle); //PCaloHits EE
  iEvent.getByToken(pCaloHits_ES_Token, pCaloHits_ES_Handle); //PCaloHits ES
  edm::PCaloHitContainer::const_iterator caloHitsItr;

  if(runOnReco_) {
     iEvent.getByToken(photon_Token_,photon_Handle);
     iEvent.getByToken(electron_Token_,electron_Handle);
     photons = photon_Handle.product();
     electrons = electron_Handle.product();
  }
  
  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);
  const CaloGeometry *geometry = caloGeometry.product();
  ebGeometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  eeGeometry = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

  runId = iEvent.id().run();
  eventId = iEvent.id().event();
  lumiId = iEvent.id().luminosityBlock();

  EBDetId* DidEB;
  EEDetId* DidEE;
  ESDetId* DidES;
  GlobalPoint cell;
  GlobalPoint cell_seed;

  //GenParticles
  for(auto &p : *genParticles){
      genParticle_pdgId = p.pdgId();
      genParticle_status = p.status();
      genParticle_energy = p.energy();
      genParticle_pt = p.pt();
      genParticle_eta = p.eta();
      genParticle_phi = p.phi();
  }
  
  recoPhoton_energy = -999.;
  recoPhoton_scEnergy = -999.;
  recoPhoton_pt = -999.;
  recoPhoton_eta = -999.;
  recoPhoton_phi = -999.;
  recoPhoton_scEta = -999.;
  recoPhoton_scPhi = -999.;
  recoPhoton_r9 = -999.;
  recoPhoton_full5x5_r9 = -999.;
  recoPhoton_nPFClusters = -999;
  recoPhoton_hasConversionTracks = false;
  recoPhoton_nConversions = -1;  
  recoPhoton_nConversionsOneLeg = -1; 
  recoPhoton_egmCutBasedPhotonIDloose = -1; 
  recoPhoton_egmCutBasedPhotonIDmedium = -1; 
  recoPhoton_egmCutBasedPhotonIDtight = -1; 
  recoPhoton_egmMVAPhotonIDmedium = -1; 
  recoPhoton_egmMVAPhotonIDtight = -1; 
  recoPhoton_passPreselections = -1; 
  recoElectron_energy = -999.;
  recoElectron_scEnergy = -999.;
  recoElectron_pt = -999.;
  recoElectron_eta = -999.;
  recoElectron_phi = -999.;
  recoElectron_scEta = -999.;
  recoElectron_scPhi = -999.;
  recoElectron_r9 = -999.;
  recoElectron_full5x5_r9 = -999.;
  recoElectron_nBrems = -999;
  recoElectron_superClusterFbrem = -999.;
  recoElectron_trackFbrem = -999.;
  recoElectron_passConversionVeto = false;
  recoElectron_nConversions = -1;  
  recoElectron_nConversionsOneLeg = -1; 
  recoElectron_egmCutBasedElectronIDVeto = -1; 
  recoElectron_egmCutBasedElectronIDloose = -1; 
  recoElectron_egmCutBasedElectronIDmedium = -1; 
  recoElectron_egmCutBasedElectronIDtight = -1; 
  recoElectron_egmMVAElectronIDloose = -1; 
  recoElectron_egmMVAElectronIDmedium = -1; 
  recoElectron_egmMVAElectronIDtight = -1; 
  recoElectron_egmMVAElectronIDlooseNoIso = -1; 
  recoElectron_egmMVAElectronIDmediumNoIso = -1; 
  recoElectron_egmMVAElectronIDtightNoIso = -1; 
  recoElectron_heepElectronID = -1; 

  int iMatched=-1;
  if(runOnReco_)
  {    
      double dRmin = 999.;
      if(abs(genParticle_pdgId)==22){
         int iPho=0;
         for(const auto& iPhoton : *photons){
             double dR = deltaR(iPhoton.superCluster()->seed()->eta(),iPhoton.superCluster()->seed()->phi(),genParticle_eta,genParticle_phi);
             if(dR<dRmin){
                dRmin = dR;
                iMatched = iPho;
             }
             iPho++;   
         }  
         if(iMatched>=0){
            reco::PhotonCoreRef phoCoreRef = photons->at(iMatched).photonCore(); 
            recoPhoton_energy = photons->at(iMatched).energy();
            recoPhoton_scEnergy = photons->at(iMatched).superCluster()->energy();
            recoPhoton_pt = photons->at(iMatched).et();
            recoPhoton_eta = photons->at(iMatched).eta();
            recoPhoton_phi = photons->at(iMatched).phi();
            recoPhoton_scEta = photons->at(iMatched).superCluster()->eta();
            recoPhoton_scPhi = photons->at(iMatched).superCluster()->phi();
            recoPhoton_r9 = photons->at(iMatched).r9();
            recoPhoton_full5x5_r9 = photons->at(iMatched).full5x5_r9(); 
            recoPhoton_nPFClusters = photons->at(iMatched).superCluster()->clustersSize();
            recoPhoton_hasConversionTracks = photons->at(iMatched).hasConversionTracks(); 
            recoPhoton_nConversions = phoCoreRef->conversions().size();  
            recoPhoton_nConversionsOneLeg = phoCoreRef->conversionsOneLeg().size(); 
            recoPhoton_egmCutBasedPhotonIDloose = photons->at(iMatched).photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
            recoPhoton_egmCutBasedPhotonIDmedium = photons->at(iMatched).photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
            recoPhoton_egmCutBasedPhotonIDtight = photons->at(iMatched).photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
            recoPhoton_egmMVAPhotonIDmedium = photons->at(iMatched).photonID("mvaPhoID-RunIIFall17-v1p1-wp90");
            recoPhoton_egmMVAPhotonIDtight = photons->at(iMatched).photonID("mvaPhoID-RunIIFall17-v2-wp80");
            recoPhoton_passPreselections = passPreselections(&photons->at(iMatched));
         }
      }else if(abs(genParticle_pdgId)==11){ 
         int iEle=0; 
         for(const auto& iElectron : *electrons){
             double dR = deltaR(iElectron.superCluster()->seed()->eta(),iElectron.superCluster()->seed()->phi(),genParticle_eta,genParticle_phi);
             if(dR<dRmin){
                dRmin = dR;
                iMatched = iEle;
             }
             iEle++;  
         }
         if(iMatched>=0){
            reco::GsfElectronCoreRef eleCoreRef = electrons->at(iMatched).core();
            recoElectron_energy = electrons->at(iMatched).energy();
            recoElectron_scEnergy = electrons->at(iMatched).superCluster()->energy();
            recoElectron_pt = electrons->at(iMatched).pt();
            recoElectron_eta = electrons->at(iMatched).eta();
            recoElectron_phi = electrons->at(iMatched).phi();
            recoElectron_scEta = electrons->at(iMatched).superCluster()->eta();
            recoElectron_scPhi = electrons->at(iMatched).superCluster()->phi();
            recoElectron_r9 = electrons->at(iMatched).r9();
            recoElectron_full5x5_r9 = electrons->at(iMatched).full5x5_r9(); 
            recoElectron_nBrems = electrons->at(iMatched).numberOfBrems();
            recoElectron_superClusterFbrem = electrons->at(iMatched).superClusterFbrem();
            recoElectron_trackFbrem = electrons->at(iMatched).trackFbrem();
            recoElectron_passConversionVeto = electrons->at(iMatched).passConversionVeto(); 
            recoElectron_nConversions = eleCoreRef->conversions().size();  
            recoElectron_nConversionsOneLeg = eleCoreRef->conversionsOneLeg().size();  
            recoElectron_egmCutBasedElectronIDVeto = electrons->at(iMatched).electronID("cutBasedElectronID-Fall17-94X-V2-veto");
            recoElectron_egmCutBasedElectronIDloose = electrons->at(iMatched).electronID("cutBasedElectronID-Fall17-94X-V2-loose");
            recoElectron_egmCutBasedElectronIDmedium = electrons->at(iMatched).electronID("cutBasedElectronID-Fall17-94X-V2-medium");
            recoElectron_egmCutBasedElectronIDtight = electrons->at(iMatched).electronID("cutBasedElectronID-Fall17-94X-V2-tight");
            recoElectron_egmMVAElectronIDloose = electrons->at(iMatched).electronID("mvaEleID-Fall17-iso-V2-wpLoose");
            recoElectron_egmMVAElectronIDmedium = electrons->at(iMatched).electronID("mvaEleID-Fall17-iso-V2-wp90");
            recoElectron_egmMVAElectronIDtight = electrons->at(iMatched).electronID("mvaEleID-Fall17-iso-V2-wp80");
            recoElectron_egmMVAElectronIDlooseNoIso = electrons->at(iMatched).electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
            recoElectron_egmMVAElectronIDmediumNoIso = electrons->at(iMatched).electronID("mvaEleID-Fall17-noIso-V2-wp90");
            recoElectron_egmMVAElectronIDtightNoIso = electrons->at(iMatched).electronID("mvaEleID-Fall17-noIso-V2-wp80");
            recoElectron_heepElectronID = electrons->at(iMatched).electronID("heepElectronID-HEEPV70");
            
         } 
      }
  } 

  simEnergy = 0.;
  simEnergy_withES = 0.;
  simR9 = 0.;
  simR9_withES = 0.;
  xtalMatrix_e3x3 = 0.; 
  xtalMatrix_totEnergy = 0.;

  energyMapEB.clear();
  energyMapEE.clear();
  energyMapES.clear();
  xtalMatrix_energy.clear();
  xtalMatrix_ieta.clear();
  xtalMatrix_iphi.clear(); 
  xtalMatrix_iz.clear();
  xtalMax_showerProfile.clear();
  xtalMax_showerProfileError.clear();
  sumMatrix_showerProfile.clear();
  sumMatrix_showerProfileError.clear();
  avgMatrix_showerProfile.clear();
  avgMatrix_showerProfileError.clear();
  for(int iXtal=0; iXtal<dim_; iXtal++)
      Profiles[iXtal]->Reset();

  xtalMax_showerProfile.resize(int(showerProfileBinning_.at(0)));
  xtalMax_showerProfileError.resize(int(showerProfileBinning_.at(0)));
  sumMatrix_showerProfile.resize(int(showerProfileBinning_.at(0)));
  sumMatrix_showerProfileError.resize(int(showerProfileBinning_.at(0)));
  avgMatrix_showerProfile.resize(int(showerProfileBinning_.at(0)));
  avgMatrix_showerProfileError.resize(int(showerProfileBinning_.at(0)));
 
  for(caloHitsItr = pCaloHits_EB_Handle->begin(); caloHitsItr != pCaloHits_EB_Handle->end(); caloHitsItr++){
      DidEB = new EBDetId(caloHitsItr->id());
      if(energyMapEB.find(DidEB->rawId()) == energyMapEB.end()){
         energyMapEB[DidEB->rawId()] = caloHitsItr->energy();
      }else{
         energyMapEB[DidEB->rawId()] += caloHitsItr->energy();
      }
  }

  for(caloHitsItr = pCaloHits_EE_Handle->begin(); caloHitsItr != pCaloHits_EE_Handle->end(); caloHitsItr++){
      DidEE = new EEDetId(caloHitsItr->id());
      if(energyMapEE.find(DidEE->rawId()) == energyMapEE.end()){
         energyMapEE[DidEE->rawId()] = caloHitsItr->energy();
      }else{
         energyMapEE[DidEE->rawId()] += caloHitsItr->energy();
      }
  }
  
  for(caloHitsItr = pCaloHits_ES_Handle->begin(); caloHitsItr != pCaloHits_ES_Handle->end(); caloHitsItr++){
      DidES = new ESDetId(caloHitsItr->id());
      if(energyMapES.find(DidES->rawId()) == energyMapES.end()){
         energyMapES[DidES->rawId()] = caloHitsItr->energy();
      }else{
         energyMapES[DidES->rawId()] += caloHitsItr->energy();
      }
  }

  std::map<long long int,float>::iterator maxValEB = std::max_element(energyMapEB.begin(),energyMapEB.end(),[] (const std::pair<long long int,float>& a, const std::pair<long long int,float>& b)->bool{ return a.second < b.second; });
  std::map<long long int,float>::iterator maxValEE = std::max_element(energyMapEE.begin(),energyMapEE.end(),[] (const std::pair<long long int,float>& a, const std::pair<long long int,float>& b)->bool{ return a.second < b.second; });
  
  bool isEB = false;
  if(maxValEB->second > maxValEE->second) isEB = true;

  long long int rawId = 0;
  float energyMax = -999.;
  if(isEB){ 
     rawId = maxValEB->first;
     energyMax = maxValEB->second;
  }else{ 
     rawId = maxValEE->first;
     energyMax = maxValEE->second; 
  }
  xtalMax_energy = energyMax;
  
  int widthMatrix = int((matrixSize_-1.)/2.); 
  xtalMatrix = EcalClusterTools::matrixDetId(topology,DetId(rawId),{-widthMatrix, widthMatrix, -widthMatrix, widthMatrix}); 
  xtalMatrix_3x3 = EcalClusterTools::matrixDetId(topology,DetId(rawId),{-1, 1, -1, 1}); 
  if(saveShowerProfileRefinedSC_ && runOnReco_){ 
     if(abs(genParticle_pdgId)==22 && iMatched>=0){ 
        xtalRefinedSC = getRefinedSCIds(photons->at(iMatched).superCluster());
     }else if(abs(genParticle_pdgId)==11 && iMatched>=0){ 
        xtalRefinedSC = getRefinedSCIds(electrons->at(iMatched).superCluster());
     }
  }

  if(isEB){ //ECAL barrel 
  
     for(caloHitsItr = pCaloHits_EB_Handle->begin(); caloHitsItr != pCaloHits_EB_Handle->end(); caloHitsItr++)
     {
         DidEB = new EBDetId(caloHitsItr->id()); 
         if(DidEB->rawId() != rawId) continue; 
         
         xtalMax_ieta = DidEB->ieta();
         xtalMax_iphi = DidEB->iphi();
         xtalMax_iz = 0; 
         
         int temp = caloHitsItr->depth(); 
         if((temp == 1) || (temp == 2)) continue;
         float depth = (caloHitsItr->depth()>>3)/depthUnit_;
         maxProfile->Fill(depth, caloHitsItr->energy());
     } 
     
     for(auto const& xtal : energyMapEB){
         simEnergy += xtal.second;
         simEnergy_withES += xtal.second;
     }
     
     for(unsigned int i=0; i<xtalMatrix_3x3.size(); i++){ 
         float xtal_energy = 0.;
         if(energyMapEB.find(xtalMatrix_3x3[i].rawId()) != energyMapEB.end()){
            xtal_energy = energyMapEB[xtalMatrix_3x3[i].rawId()]; 
         }
         xtalMatrix_e3x3 += xtal_energy; 
     }   
        
     if(saveShowerProfileRefinedSC_ && runOnReco_ && iMatched>=0){

        //std::cout << "xtalRefinedSC.size(): " << xtalRefinedSC.size() << " - " << electrons->at(iMatched).full5x5_r9() << std::endl;
        for(unsigned int i=0; i<xtalRefinedSC.size(); i++)
        {
            xtalMatrix_totEnergy += energyMapEB[xtalRefinedSC[i].rawId()];  
            xtalMatrix_ieta.push_back(EBDetId(xtalRefinedSC[i].rawId()).ieta());
            xtalMatrix_iphi.push_back(EBDetId(xtalRefinedSC[i].rawId()).iphi());
            xtalMatrix_iz.push_back(EBDetId(xtalRefinedSC[i].rawId()).zside());
            xtalMatrix_energy.push_back(energyMapEB[xtalRefinedSC[i].rawId()]);
            for(caloHitsItr = pCaloHits_EB_Handle->begin(); caloHitsItr != pCaloHits_EB_Handle->end(); caloHitsItr++)
            {
                DidEB = new EBDetId(caloHitsItr->id());
                if(DidEB->rawId() != xtalRefinedSC[i].rawId()) continue;

                int temp = caloHitsItr->depth(); 
                if((temp == 1) || (temp == 2)) continue;
                float depth = (caloHitsItr->depth()>>3)/depthUnit_;
                
                Profiles[i]->Fill(depth, caloHitsItr->energy());
                sumProfile->Fill(depth, caloHitsItr->energy());              
            }
        } 
        
     }else{
     
        for(unsigned int i=0; i<xtalMatrix.size(); i++)
        {
            xtalMatrix_totEnergy += energyMapEB[xtalMatrix[i].rawId()];  
            xtalMatrix_ieta.push_back(EBDetId(xtalMatrix[i].rawId()).ieta());
            xtalMatrix_iphi.push_back(EBDetId(xtalMatrix[i].rawId()).iphi());
            xtalMatrix_iz.push_back(EBDetId(xtalMatrix[i].rawId()).zside());
            xtalMatrix_energy.push_back(energyMapEB[xtalMatrix[i].rawId()]);
            for(caloHitsItr = pCaloHits_EB_Handle->begin(); caloHitsItr != pCaloHits_EB_Handle->end(); caloHitsItr++)
            {
                DidEB = new EBDetId(caloHitsItr->id());           
                if(DidEB->rawId() != xtalMatrix[i].rawId()) continue;
             
                int temp = caloHitsItr->depth(); 
                if((temp == 1) || (temp == 2)) continue;
                float depth = (caloHitsItr->depth()>>3)/depthUnit_;
                
                Profiles[i]->Fill(depth, caloHitsItr->energy());
                sumProfile->Fill(depth, caloHitsItr->energy()); 
           } 
        } 
         
     }
     
  }else{ //ECAL endcaps
  
     for(caloHitsItr = pCaloHits_EE_Handle->begin(); caloHitsItr != pCaloHits_EE_Handle->end(); caloHitsItr++)
     {
         DidEE = new EEDetId(caloHitsItr->id()); 
         if(DidEE->rawId() != rawId) continue; 
          
         xtalMax_ieta = DidEE->ix();
         xtalMax_iphi = DidEE->iy();
         if(DidEE->zside()>0) xtalMax_iz = 1;
         else  xtalMax_iz = -1;
                      
         int temp = caloHitsItr->depth(); 
         if((temp == 1) || (temp == 2)) continue;
         float depth = (caloHitsItr->depth()>>3)/depthUnit_;
         maxProfile->Fill(depth, caloHitsItr->energy());
     } 
     
     for(auto const& xtal : energyMapEE){
         simEnergy += xtal.second;
         simEnergy_withES += xtal.second;
     }
     for(auto const& xtal : energyMapES)
         simEnergy_withES += xtal.second;

     for(unsigned int i=0; i<xtalMatrix_3x3.size(); i++){
         float xtal_energy = 0.;
         if(energyMapEE.find(xtalMatrix_3x3[i].rawId()) != energyMapEE.end()){
            xtal_energy = energyMapEE[xtalMatrix_3x3[i].rawId()]; 
         }
         xtalMatrix_e3x3 += xtal_energy;
     }
     
     if(saveShowerProfileRefinedSC_ && runOnReco_ && iMatched>=0){
     
        //std::cout << "xtalRefinedSC.size(): " << xtalRefinedSC.size() << " - " << electrons->at(iMatched).full5x5_r9() << std::endl;
        for(unsigned int i=0; i<xtalRefinedSC.size(); i++)
        {
            int iz = 0;
            if(EEDetId(xtalRefinedSC[i].rawId()).zside()>0) iz = 1;
            else  iz = -1;
         
            xtalMatrix_totEnergy += energyMapEE[xtalRefinedSC[i].rawId()];  
            xtalMatrix_ieta.push_back(EEDetId(xtalRefinedSC[i].rawId()).ix());
            xtalMatrix_iphi.push_back(EEDetId(xtalRefinedSC[i].rawId()).iy());
            xtalMatrix_iz.push_back(iz);
            xtalMatrix_energy.push_back(energyMapEE[xtalRefinedSC[i].rawId()]);
            for(caloHitsItr = pCaloHits_EE_Handle->begin(); caloHitsItr != pCaloHits_EE_Handle->end(); caloHitsItr++)
            {
                DidEE = new EEDetId(caloHitsItr->id());
                if(DidEE->rawId() != xtalRefinedSC[i].rawId()) continue;

                int temp = caloHitsItr->depth(); 
                if((temp == 1) || (temp == 2)) continue;
                float depth = (caloHitsItr->depth()>>3)/depthUnit_;
                
                Profiles[i]->Fill(depth, caloHitsItr->energy());
                sumProfile->Fill(depth, caloHitsItr->energy());              
            }
        } 
        
     }else{
     
        for(unsigned int i=0; i<xtalMatrix.size(); i++)
        {
            int iz = 0;
            if(EEDetId(xtalMatrix[i].rawId()).zside()>0) iz = 1;
            else  iz = -1;
             
            xtalMatrix_totEnergy += energyMapEE[xtalMatrix[i].rawId()];  
            xtalMatrix_ieta.push_back(EEDetId(xtalMatrix[i].rawId()).ix());
            xtalMatrix_iphi.push_back(EEDetId(xtalMatrix[i].rawId()).iy());
            xtalMatrix_iz.push_back(iz);
            xtalMatrix_energy.push_back(energyMapEE[xtalMatrix[i].rawId()]);
            for(caloHitsItr = pCaloHits_EE_Handle->begin(); caloHitsItr != pCaloHits_EE_Handle->end(); caloHitsItr++)
            {
                DidEE = new EEDetId(caloHitsItr->id());       
                if(DidEE->rawId() != xtalMatrix[i].rawId()) continue;
             
                int temp = caloHitsItr->depth(); 
                if((temp == 1) || (temp == 2)) continue;
                float depth = (caloHitsItr->depth()>>3)/depthUnit_;
                
                Profiles[i]->Fill(depth, caloHitsItr->energy());
                sumProfile->Fill(depth, caloHitsItr->energy()); 
           } 
        } 
         
     }
  }
  
  simR9 = xtalMatrix_e3x3/simEnergy;
  simR9_withES = xtalMatrix_e3x3/simEnergy_withES; 
  
  if(saveShowerProfile_){
     for(int bin=1; bin<=maxProfile->GetNbinsX(); bin++){
         float binContent = maxProfile->GetBinContent(bin);
         float binError = maxProfile->GetBinError(bin);
         binContent = std::round(binContent * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);
         binError = std::round(binError * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);    
         xtalMax_showerProfile[bin-1] = binContent;
         xtalMax_showerProfileError[bin-1] = binError;
     }
     for(int bin=1; bin<=sumProfile->GetNbinsX(); bin++){
         float binContent = sumProfile->GetBinContent(bin);
         float binError = sumProfile->GetBinError(bin);
         binContent = std::round(binContent * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);
         binError = std::round(binError * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);    
         sumMatrix_showerProfile[bin-1] = binContent;
         sumMatrix_showerProfileError[bin-1] = binError;
     }
    
     fillAvgProfile(); 
     for(int bin=1; bin<=avgProfile->GetNbinsX(); bin++){
         float binContent = avgProfile->GetBinContent(bin);
         float binError = avgProfile->GetBinError(bin);
         binContent = std::round(binContent * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);
         binError = std::round(binError * TMath::Power(10.,precision_)) / TMath::Power(10.,precision_);    
         avgMatrix_showerProfile[bin-1] = binContent;
         avgMatrix_showerProfileError[bin-1] = binError;
     }
  }
  
  int binmax = maxProfile->GetMaximumBin(); 
  float x = maxProfile->GetXaxis()->GetBinCenter(binmax);
  float x_min = x-1*maxProfile->GetRMS();
  float x_max = x+1*maxProfile->GetRMS();
  f = new TF1("f","[0]*TMath::Power(x,[1]*[2])*TMath::Exp(-1*[2]*x)",x_min,x_max);
  f->SetLineColor(kRed+1);
  f->SetLineWidth(2);
  f->SetParameter(0,0.001);
  f->SetParLimits(0,0.,10000.); 
  f->SetParameter(1,x);
  f->SetParLimits(1,0.,showerProfileBinning_.at(2));  
  f->SetParameter(2,0.005*depthUnit_); 
  f->SetParLimits(2,0.,10000.);
  
  //TVirtualFitter::SetDefaultFitter("Minuit2");
  int status = 0;
  TFitResultPtr rp;  
  rp = maxProfile->Fit("f", "QMERS+");
  while(status != 0)
  {
      rp = maxProfile->Fit("f", "QMERS+");
      status = rp;
      if(status == 0) break;
  }
  
  x_min = x-2*maxProfile->GetRMS();
  x_max = x+2*maxProfile->GetRMS();
  f->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));
  f->SetParLimits(0,f->GetParameter(0)-0.5*f->GetParError(0),f->GetParameter(0)+0.5*f->GetParError(0)); 
  f->SetParLimits(1,f->GetParameter(1)-0.1*f->GetParError(1),f->GetParameter(1)+0.1*f->GetParError(1)); 
  f->SetParLimits(2,f->GetParameter(2)-0.5*f->GetParError(2),f->GetParameter(2)+0.5*f->GetParError(2)); 
  status = 0;
  rp = maxProfile->Fit("f", "QMERS+");
  while(status != 0)
  {
      rp = maxProfile->Fit("f", "QMERS+");
      status = rp;
      if(status == 0) break;
  }
 
  xtalMax_tMax = f->GetParameter(1);  
  xtalMax_tMaxError = f->GetParError(1);  
  xtalMax_tMax_fitPar0 = f->GetParameter(0);
  xtalMax_tMax_fitPar1 = f->GetParameter(1);
  xtalMax_tMax_fitPar2 = f->GetParameter(2);
  xtalMax_tMax_fitParError0 = f->GetParError(0);
  xtalMax_tMax_fitParError1 = f->GetParError(1);
  xtalMax_tMax_fitParError2 = f->GetParError(2);
  xtalMax_tMax_fitStatus = status;

  tree->Fill(); 
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void simAnalysis::beginJob() {
  // please remove this method if not needed
  //Profiles = new TH1D("Profiles","Profiles",260,0., 26.);
  //Profiles = fs->make<TH1D>("Profiles","Profiles",260,0., 26.);  
  //Profiles->Sumw2(); 
}

// ------------ method called once each job just after ending the event loop  ------------
void simAnalysis::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void simAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //Fesc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //Fescriptions.addDefault(desc);
}

void simAnalysis::setTree(TTree* tree)
{
   tree->Branch("eventId",&eventId,"eventId/L");
   tree->Branch("lumiId",&lumiId,"lumiId/L");
   tree->Branch("runId",&runId,"runId/L");
   tree->Branch("genParticle_pdgId", &genParticle_pdgId, "genParticle_pdgId/I");
   tree->Branch("genParticle_status", &genParticle_status, "genParticle_status/F"); 
   tree->Branch("genParticle_energy", &genParticle_energy, "genParticle_energy/F");
   tree->Branch("genParticle_pt", &genParticle_pt, "genParticle_pt/F");
   tree->Branch("genParticle_eta", &genParticle_eta, "genParticle_eta/F");
   tree->Branch("genParticle_phi", &genParticle_phi, "genParticle_phi/F");
   tree->Branch("simEnergy", &simEnergy, "simEnergy/F");
   tree->Branch("simEnergy_withES", &simEnergy_withES, "simEnergy_withES/F");
   tree->Branch("simR9", &simR9, "simR9/F");
   tree->Branch("simR9_withES", &simR9_withES, "simR9_withES/F");
   //tree->Branch("xtalMatrix_ieta","std::vector<int>",&xtalMatrix_ieta);
   //tree->Branch("xtalMatrix_iphi","std::vector<int>",&xtalMatrix_iphi);
   //tree->Branch("xtalMatrix_iz","std::vector<int>",&xtalMatrix_iz);
   //tree->Branch("xtalMatrix_energy","std::vector<float>",&xtalMatrix_energy);
   tree->Branch("xtalMatrix_totEnergy", &xtalMatrix_totEnergy, "xtalMatrix_totEnergy/F");
   tree->Branch("xtalMatrix_e3x3", &xtalMatrix_e3x3, "xtalMatrix_e3x3/F");
   if(saveShowerProfile_){ 
      tree->Branch("xtalMax_showerProfile","std::vector<float>",&xtalMax_showerProfile);
      tree->Branch("xtalMax_showerProfileError","std::vector<float>",&xtalMax_showerProfileError);
      tree->Branch("sumMatrix_showerProfile","std::vector<float>",&sumMatrix_showerProfile);
      tree->Branch("sumMatrix_showerProfileError","std::vector<float>",&sumMatrix_showerProfileError);
      tree->Branch("avgMatrix_showerProfile","std::vector<float>",&avgMatrix_showerProfile);
      tree->Branch("avgMatrix_showerProfileError","std::vector<float>",&avgMatrix_showerProfileError);
   }
   tree->Branch("xtalMax_ieta",&xtalMax_ieta,"xtalMax_ieta/I");
   tree->Branch("xtalMax_iphi",&xtalMax_iphi,"xtalMax_iphi/I");
   tree->Branch("xtalMax_iz",&xtalMax_iz,"xtalMax_iz/I");
   tree->Branch("xtalMax_energy",&xtalMax_energy,"xtalMax_energy/F");
   tree->Branch("xtalMax_binMaxVal",&xtalMax_binMaxVal,"xtalMax_binMaxVal/F");
   tree->Branch("xtalMax_tMax",&xtalMax_tMax,"xtalMax_tMax/F");
   tree->Branch("xtalMax_tMaxError",&xtalMax_tMaxError,"xtalMax_tMaxError/F");
   tree->Branch("xtalMax_tMax_fitPar0",&xtalMax_tMax_fitPar0,"xtalMax_tMax_fitPar0/F"); 
   tree->Branch("xtalMax_tMax_fitPar1",&xtalMax_tMax_fitPar1,"xtalMax_tMax_fitPar1/F"); 
   tree->Branch("xtalMax_tMax_fitPar2",&xtalMax_tMax_fitPar2,"xtalMax_tMax_fitPar2/F"); 
   tree->Branch("xtalMax_tMax_fitParError0",&xtalMax_tMax_fitParError0,"xtalMax_tMax_fitParError0/F"); 
   tree->Branch("xtalMax_tMax_fitParError1",&xtalMax_tMax_fitParError1,"xtalMax_tMax_fitParError1/F"); 
   tree->Branch("xtalMax_tMax_fitParError2",&xtalMax_tMax_fitParError2,"xtalMax_tMax_fitParError2/F");  
   tree->Branch("xtalMax_tMax_fitStatus",&xtalMax_tMax_fitStatus,"xtalMax_tMax_fitStatus/I");  
   if(runOnReco_){
      tree->Branch("recoPhoton_energy",&recoPhoton_energy,"recoPhoton_energy/F");
      tree->Branch("recoPhoton_scEnergy",&recoPhoton_scEnergy,"recoPhoton_scEnergy/F");
      tree->Branch("recoPhoton_pt",&recoPhoton_pt,"recoPhoton_pt/F");
      tree->Branch("recoPhoton_eta",&recoPhoton_eta,"recoPhoton_eta/F");
      tree->Branch("recoPhoton_phi",&recoPhoton_phi,"recoPhoton_phi/F");
      tree->Branch("recoPhoton_scEta",&recoPhoton_scEta,"recoPhoton_scEta/F");
      tree->Branch("recoPhoton_scPhi",&recoPhoton_scPhi,"recoPhoton_scPhi/F");
      tree->Branch("recoPhoton_r9",&recoPhoton_r9,"recoPhoton_r9/F");
      tree->Branch("recoPhoton_full5x5_r9",&recoPhoton_full5x5_r9,"recoPhoton_full5x5_r9/F");
      tree->Branch("recoPhoton_nPFClusters",&recoPhoton_nPFClusters,"recoPhoton_nPFClusters/I");
      tree->Branch("recoPhoton_hasConversionTracks",&recoPhoton_hasConversionTracks,"recoPhoton_hasConversionTracks/O");
      tree->Branch("recoPhoton_nConversions",&recoPhoton_nConversions,"recoPhoton_nConversions/I");
      tree->Branch("recoPhoton_nConversionsOneLeg",&recoPhoton_nConversionsOneLeg,"recoPhoton_nConversionsOneLeg/I");
      tree->Branch("recoPhoton_egmCutBasedPhotonIDloose",&recoPhoton_egmCutBasedPhotonIDloose,"recoPhoton_egmCutBasedPhotonIDloose/I");
      tree->Branch("recoPhoton_egmCutBasedPhotonIDmedium",&recoPhoton_egmCutBasedPhotonIDmedium,"recoPhoton_egmCutBasedPhotonIDmedium/I");
      tree->Branch("recoPhoton_egmCutBasedPhotonIDtight",&recoPhoton_egmCutBasedPhotonIDtight,"recoPhoton_egmCutBasedPhotonIDtight/I");
      tree->Branch("recoPhoton_egmMVAPhotonIDmedium",&recoPhoton_egmMVAPhotonIDmedium,"recoPhoton_egmMVAPhotonIDmedium/I");
      tree->Branch("recoPhoton_egmMVAPhotonIDtight",&recoPhoton_egmMVAPhotonIDtight,"recoPhoton_egmMVAPhotonIDtight/I");
      tree->Branch("recoPhoton_passPreselections",&recoPhoton_passPreselections,"recoPhoton_passPreselections/I");
      tree->Branch("recoElectron_energy",&recoElectron_energy,"recoElectron_energy/F");
      tree->Branch("recoElectron_scEnergy",&recoElectron_scEnergy,"recoElectron_scEnergy/F");
      tree->Branch("recoElectron_pt",&recoElectron_pt,"recoElectron_pt/F");
      tree->Branch("recoElectron_eta",&recoElectron_eta,"recoElectron_eta/F");
      tree->Branch("recoElectron_phi",&recoElectron_phi,"recoElectron_phi/F");
      tree->Branch("recoElectron_scEta",&recoElectron_scEta,"recoElectron_scEta/F");
      tree->Branch("recoElectron_scPhi",&recoElectron_scPhi,"recoElectron_scPhi/F");
      tree->Branch("recoElectron_r9",&recoElectron_r9,"recoElectron_r9/F");
      tree->Branch("recoElectron_full5x5_r9",&recoElectron_full5x5_r9,"recoElectron_full5x5_r9/F");
      tree->Branch("recoElectron_nBrems",&recoElectron_nBrems,"recoElectron_nBrems/I");
      tree->Branch("recoElectron_superClusterFbrem",&recoElectron_superClusterFbrem,"recoElectron_superClusterFbrem/F");
      tree->Branch("recoElectron_trackFbrem",&recoElectron_trackFbrem,"recoElectron_trackFbrem/F");
      tree->Branch("recoElectron_passConversionVeto",&recoElectron_passConversionVeto,"recoElectron_passConversionVeto/O");
      tree->Branch("recoElectron_nConversions",&recoElectron_nConversions,"recoElectron_nConversions/I");
      tree->Branch("recoElectron_nConversionsOneLeg",&recoElectron_nConversionsOneLeg,"recoElectron_nConversionsOneLeg/I");
      tree->Branch("recoElectron_egmCutBasedElectronIDVeto",&recoElectron_egmCutBasedElectronIDVeto,"recoElectron_egmCutBasedElectronIDVeto/I");  
      tree->Branch("recoElectron_egmCutBasedElectronIDloose",&recoElectron_egmCutBasedElectronIDloose,"recoElectron_egmCutBasedElectronIDloose/I");    
      tree->Branch("recoElectron_egmCutBasedElectronIDmedium",&recoElectron_egmCutBasedElectronIDmedium,"recoElectron_egmCutBasedElectronIDmedium/I");    
      tree->Branch("recoElectron_egmCutBasedElectronIDtight",&recoElectron_egmCutBasedElectronIDtight,"recoElectron_egmCutBasedElectronIDtight/I");    
      tree->Branch("recoElectron_egmMVAElectronIDloose",&recoElectron_egmMVAElectronIDloose,"recoElectron_egmMVAElectronIDloose/I");    
      tree->Branch("recoElectron_egmMVAElectronIDmedium",&recoElectron_egmMVAElectronIDmedium,"recoElectron_egmMVAElectronIDmedium/I");    
      tree->Branch("recoElectron_egmMVAElectronIDtight",&recoElectron_egmMVAElectronIDtight,"recoElectron_egmMVAElectronIDtight/I");    
      tree->Branch("recoElectron_egmMVAElectronIDlooseNoIso",&recoElectron_egmMVAElectronIDlooseNoIso,"recoElectron_egmMVAElectronIDlooseNoIso/I");    
      tree->Branch("recoElectron_egmMVAElectronIDmediumNoIso",&recoElectron_egmMVAElectronIDmediumNoIso,"recoElectron_egmMVAElectronIDmediumNoIso/I");    
      tree->Branch("recoElectron_egmMVAElectronIDtightNoIso",&recoElectron_egmMVAElectronIDtightNoIso,"recoElectron_egmMVAElectronIDtightNoIso/I");    
      tree->Branch("recoElectron_heepElectronID",&recoElectron_heepElectronID,"recoElectron_heepElectronID/I");  
 
   }
}

bool simAnalysis::passPreselections(const pat::Photon* photon)
{
   bool pass = false;

   float pt = photon->pt();
   float scEta = photon->superCluster()->eta(); 
   float r9 = photon->full5x5_r9();
   float hoe = photon->hadronicOverEm(); 
   float sieie = photon->full5x5_sigmaIetaIeta();
   float pfPhoIso03 = photon->photonIso(); 
   float trackIso03 = photon->trkSumPtHollowConeDR03(); 
   float egCHIso = photon->chargedHadronIso();
   bool passEleVeto = photon->passElectronVeto();

   if( (( abs(scEta)<1.479 && r9>0.5 && r9>0.85 ) || 
        ( abs(scEta)<1.479 && r9>0.5 && r9<=0.85 && sieie<0.015 && pfPhoIso03<4.0 && trackIso03<6.0 ) ||
        ( abs(scEta)>1.479 && r9>0.8 && r9>0.90 ) ||
        ( abs(scEta)>1.479 && r9>0.8 && r9<=0.90 && sieie<0.035 && pfPhoIso03<4.0 && trackIso03<6.0 )) 
       && (( r9>0.8 && egCHIso<20. ) || ( egCHIso/pt<0.3 && pt>14. && hoe<0.15 ))
       && ( passEleVeto>0 && hoe<0.08 ) 
   ) pass = true;
   
   return pass;
}   

std::vector<DetId> simAnalysis::getRefinedSCIds(reco::SuperClusterRef sc)
{
   std::vector<DetId> detIds;
   std::vector<DetId> detIds_final;
   for(reco::CaloCluster_iterator iBC = sc->clustersBegin(); iBC != sc->clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &clrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < clrechits.size(); i++){ 
            if(std::find(detIds.begin(),detIds.end(),clrechits.at(i).first)==detIds.end()) detIds.push_back(clrechits.at(i).first);
        }
   } 
 
   detIds_final.resize(detIds.size()); 
   for(unsigned int i=0; i<detIds.size(); i++)
       detIds_final[i] = detIds.at(i);

   return detIds_final;
}  

double simAnalysis::getMean(std::vector<double>* vals, std::vector<double>* weights)
{
   double sumV=0.;
   double sumW=0.;
   for(unsigned int i=0; i<vals->size(); i++){ 
       sumW += weights->at(i); 
       sumV += weights->at(i)*vals->at(i);
   }
   double mean = 0.;
   if(sumW!=0.) mean = sumV/sumW;
   return mean;
}

double simAnalysis::getError(std::vector<double>* vals, std::vector<double>* weights)
{
   int M=0;
   double sumV=0.;
   double sumW=0.;
   double mean=getMean(vals,weights);
   for(unsigned int i=0; i<vals->size(); i++){ 
       if(weights->at(i)!=0) M++;
       sumW += weights->at(i); 
       sumV += weights->at(i)*(vals->at(i)-mean)*(vals->at(i)-mean);
   }

   double error = 0.;
   double norm = (double(M)-1.)/double(M);
   if(M!=0 && M!=1 && sumW!=0.) error = sqrt(sumV/(norm*sumW));
   return (error/TMath::Sqrt(vals->size()));
}

void simAnalysis::fillAvgProfile()
{
   for(int bin=1; bin<=avgProfile->GetNbinsX(); bin++)
   {
       vals->clear();
       weights->clear();
       for(int i=0; i<dim_; i++)
       { 
           if(Profiles[i]->GetEntries()==0) continue; 
           if(Profiles[i]->GetBinContent(bin)!=0. && Profiles[i]->GetBinError(bin)!=0.){    
              vals->push_back(Profiles[i]->GetBinContent(bin));
              weights->push_back(Profiles[i]->GetBinError(bin));
           } 
       } 
       double mean = 0.;
       double error = 0.; 
       if(vals->size()!=0 && weights->size()!=0){ 
          mean = getMean(vals,weights);
          error = getError(vals,weights); 
       }
       avgProfile->SetBinContent(bin,mean);
       avgProfile->SetBinError(bin,error);
   }   
}

//Fefine this as a plug-in
DEFINE_FWK_MODULE(simAnalysis);

