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
#include "GENSIMAnalysis/Dumpers/plugins/ElePhoDumper.h"

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
#include "TCanvas.h"
#include "TMath.h"
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

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
ElePhoDumper::ElePhoDumper(const edm::ParameterSet& iConfig)
{
   pileupSummaryToken_            = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   patElectronToken_              = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectronCollection"));
   patPhotonToken_                = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("patPhotonCollection")); 
   
   nBits_                         = iConfig.getParameter<int>("nBits");
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   isMC_                          = iConfig.getParameter<bool>("isMC"); 
   isDY_                          = iConfig.getParameter<bool>("isDY"); 
   savePhotons_                   = iConfig.getParameter<bool>("savePhotons"); 
   saveElectrons_                 = iConfig.getParameter<bool>("saveElectrons");

   egmCutBasedElectronIDVeto_     = iConfig.getParameter<std::string>("egmCutBasedElectronIDVeto"); 
   egmCutBasedElectronIDloose_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDloose");  
   egmCutBasedElectronIDmedium_   = iConfig.getParameter<std::string>("egmCutBasedElectronIDmedium"); 
   egmCutBasedElectronIDtight_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDtight");   
   egmMVAElectronIDloose_         = iConfig.getParameter<std::string>("egmMVAElectronIDloose");  
   egmMVAElectronIDmedium_        = iConfig.getParameter<std::string>("egmMVAElectronIDmedium"); 
   egmMVAElectronIDtight_         = iConfig.getParameter<std::string>("egmMVAElectronIDtight");   
   egmMVAElectronIDlooseNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDlooseNoIso");  
   egmMVAElectronIDmediumNoIso_   = iConfig.getParameter<std::string>("egmMVAElectronIDmediumNoIso"); 
   egmMVAElectronIDtightNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDtightNoIso");
   heepElectronID_                = iConfig.getParameter<std::string>("heepElectronID");   
   egmCutBasedPhotonIDloose_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDloose");  
   egmCutBasedPhotonIDmedium_     = iConfig.getParameter<std::string>("egmCutBasedPhotonIDmedium"); 
   egmCutBasedPhotonIDtight_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDtight"); 
   egmMVAPhotonIDmedium_          = iConfig.getParameter<std::string>("egmMVAPhotonIDmedium"); 
   egmMVAPhotonIDtight_           = iConfig.getParameter<std::string>("egmMVAPhotonIDtight");
  
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

ElePhoDumper::~ElePhoDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void ElePhoDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   //MC-only info and collections
   truePU=-1.;
   obsPU=-1.;
   edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   edm::Handle<std::vector<reco::GenParticle> > genParticles; 

   if(isMC_){    
      ev.getByToken(pileupSummaryToken_, PupInfo);
      if(PupInfo.isValid()) 
      {
         for(auto &pu : *PupInfo){
           if(pu.getBunchCrossing() == 0 ){
              truePU = pu.getTrueNumInteractions();
              obsPU = pu.getPU_NumInteractions();
              break;
           } 
         } 
      }else{
         std::cerr << "Analyze --> PupInfo not found" << std::endl;
      }

      ev.getByToken(genToken_,genParticles);
      if(!genParticles.isValid()) {
         std::cerr << "Analyze --> genParticles not found" << std::endl; 
         return;
      }
    
   }

   //Other collections
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }
   
   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
       std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
       return;
   } 
    
   edm::Handle<std::vector<pat::Photon> > patPhoton;
   ev.getByToken(patPhotonToken_,patPhoton);
   if(!patPhoton.isValid()) {
      std::cerr << "Analyze --> patPhotons not found" << std::endl; 
      return;
   }
   
   edm::Handle<std::vector<pat::Electron> > patElectron;
   ev.getByToken(patElectronToken_,patElectron);
   if (!patElectron.isValid()) {
       std::cerr << "Analyze --> patElectrons not found" << std::endl; 
       return;
   }
   
   edm::ESHandle<CaloTopology> caloTopology;
   iSetup.get<CaloTopologyRecord>().get(caloTopology);
   const CaloTopology* topology = caloTopology.product();
   
   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   nVtx = vertices->size();
   rho = *(rhos.product());

   if(isMC_ && !isDY_) genParticle_size = (*(genParticles.product())).size();
   patElectron_size = (*(patElectron.product())).size();
   patPhoton_size = (*(patPhoton.product())).size();
   setVectors(patElectron_size,patPhoton_size);
   
   int nGenParts = 0;
   if(isMC_){ 
      const std::vector<reco::GenParticle>& genParts = *(genParticles.product());
      for(unsigned int iGen=0; iGen<genParts.size(); iGen++)
      {
          if(isDY_ && abs(genParts.at(iGen).pdgId())==11 && genParts.at(iGen).status()==1 && genParts.at(iGen).energy()>0.) nGenParts++;
      }  
      if(isDY_ && nGenParts<2) return;
      for(unsigned int iGen=0; iGen<genParts.size(); iGen++)
      {
          if(isDY_ && (abs(genParts.at(iGen).pdgId())!=11 || genParts.at(iGen).status()!=1 || genParts.at(iGen).energy()<=0.)) continue;
          genParticle_pdgId.push_back(genParts.at(iGen).pdgId()); 
          genParticle_status.push_back(genParts.at(iGen).status()); 
          genParticle_statusFlag.push_back(getGenStatusFlag(&genParts.at(iGen)));
          genParticle_statusFlag_isLastCopyBeforeFSR.push_back(genParts.at(iGen).statusFlags().isLastCopyBeforeFSR());               
          genParticle_statusFlag_isLastCopy.push_back(genParts.at(iGen).statusFlags().isLastCopy());    
          genParticle_statusFlag_isFirstCopy.push_back(genParts.at(iGen).statusFlags().isFirstCopy());    
          genParticle_statusFlag_fromHardProcessBeforeFSR.push_back(genParts.at(iGen).statusFlags().isFirstCopy());    
          genParticle_statusFlag_isDirectHardProcessTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isDirectHardProcessTauDecayProduct());    
          genParticle_statusFlag_isHardProcessTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isHardProcessTauDecayProduct());    
          genParticle_statusFlag_fromHardProcess.push_back(genParts.at(iGen).statusFlags().fromHardProcess());    
          genParticle_statusFlag_isHardProcess.push_back(genParts.at(iGen).statusFlags().isHardProcess());    
          genParticle_statusFlag_isDirectHadronDecayProduct.push_back(genParts.at(iGen).statusFlags().isDirectHadronDecayProduct());    
          genParticle_statusFlag_isDirectPromptTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isDirectPromptTauDecayProduct());    
          genParticle_statusFlag_isDirectTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isDirectTauDecayProduct());    
          genParticle_statusFlag_isPromptTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isPromptTauDecayProduct());    
   	  genParticle_statusFlag_isTauDecayProduct.push_back(genParts.at(iGen).statusFlags().isTauDecayProduct());    
   	  genParticle_statusFlag_isDecayedLeptonHadron.push_back(genParts.at(iGen).statusFlags().isDecayedLeptonHadron());    
   	  genParticle_statusFlag_isPrompt.push_back(genParts.at(iGen).statusFlags().isPrompt());  
          genParticle_energy.push_back(genParts.at(iGen).energy()); 
          genParticle_pt.push_back(genParts.at(iGen).pt());
          genParticle_eta.push_back(genParts.at(iGen).eta());
          genParticle_phi.push_back(genParts.at(iGen).phi());    
      } 
      if(isDY_) genParticle_size = nGenParts;  
   }
          
   //save electron info
   if(saveElectrons_)
   { 
      int iEle=0;
      for(const auto& iElectron : *(patElectron.product())){ 
 
          reco::SuperClusterRef scRef = iElectron.superCluster();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
          
          double swissCross = -999.;
          double e3x3 = -999.;
          //double e2x2 = -999.;
          //double e5x5 = -999.;  
          double eMax = -999.;
          double full5x5_e3x3 = -999.;
          //double full5x5_e2x2 = -999.;
          //double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 
         
          reco::GsfElectron::ShowerShape eleSS = iElectron.showerShape();
          reco::GsfElectron::ShowerShape full5x5_eleSS = iElectron.full5x5_showerShape();
          const std::vector<std::pair<DetId,float> > &hits= iElectron.superCluster()->hitsAndFractions();
          if(iElectron.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             //e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             //e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));             
             //full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             //full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             //full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iElectron.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             //e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             //e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             //full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             //full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
          }
              
          patElectron_index[iEle] = iEle;
          patElectron_classification[iEle] = iElectron.classification();
          patElectron_scNPFClusters[iEle] = scRef->clusters().size();
          patElectron_charge[iEle] = iElectron.charge(); 
          patElectron_isEB[iEle] = iElectron.isEB(); 
          patElectron_isEE[iEle] = iElectron.isEE();  
          patElectron_isEBEEGap[iEle] = iElectron.isEBEEGap();  
          patElectron_isEBEtaGap[iEle] = iElectron.isEBEtaGap();  
          patElectron_isEBPhiGap[iEle] = iElectron.isEBPhiGap();  
          patElectron_isEEDeeGap[iEle] = iElectron.isEEDeeGap();  
          patElectron_isEERingGap[iEle] = iElectron.isEERingGap();   
          patElectron_isEcalDriven[iEle] = iElectron.ecalDrivenSeed(); 
          patElectron_isTrackerDriven[iEle] = iElectron.trackerDrivenSeed(); 
          patElectron_passConversionVeto[iEle] = iElectron.passConversionVeto(); 
          patElectron_eta[iEle] = reduceFloat(iElectron.p4().eta(),nBits_);
          patElectron_phi[iEle] = reduceFloat(iElectron.p4().phi(),nBits_);
          patElectron_p[iEle] = reduceFloat(iElectron.trackMomentumAtVtx().R(),nBits_);
          patElectron_pt[iEle] = reduceFloat(TMath::Sqrt(iElectron.trackMomentumAtVtx().Perp2()),nBits_);
          patElectron_pIn[iEle] = reduceFloat(iElectron.trackMomentumAtVtx().R(),nBits_);
          patElectron_pOut[iEle] = reduceFloat(iElectron.trackMomentumOut().R(),nBits_);
          patElectron_trackFbrem[iEle] = reduceFloat(iElectron.trackFbrem(),nBits_);
          patElectron_superClusterFbrem[iEle] = reduceFloat(iElectron.superClusterFbrem(),nBits_); 
          patElectron_energy[iEle] = reduceFloat(iElectron.energy(),nBits_);
          patElectron_energyErr[iEle] = reduceFloat(iElectron.p4Error(reco::GsfElectron::P4_COMBINATION),nBits_);
          patElectron_ecalEnergy[iEle] = reduceFloat(iElectron.ecalEnergy(),nBits_);
          patElectron_ecalEnergyErr[iEle] = reduceFloat(iElectron.ecalEnergyError(),nBits_);
          patElectron_et[iEle] = reduceFloat(iElectron.p4().Et(),nBits_);
          patElectron_HoE[iEle] = reduceFloat(iElectron.hadronicOverEm(),nBits_); 
          patElectron_scEta[iEle] = reduceFloat(scRef->eta(),nBits_);
          patElectron_scPhi[iEle] = reduceFloat(scRef->phi(),nBits_);
          patElectron_scEnergy[iEle] = reduceFloat(scRef->energy(),nBits_); 
          patElectron_scRawEnergy[iEle] = reduceFloat(scRef->rawEnergy(),nBits_); 
          patElectron_scRawESEnergy[iEle] = reduceFloat(scRef->preshowerEnergy(),nBits_); 
          patElectron_scEt[iEle] = reduceFloat(scRef->energy()*(Rt/R),nBits_);
          patElectron_scPhiWidth[iEle] = reduceFloat(scRef->phiWidth(),nBits_);
          patElectron_scEtaWidth[iEle] = reduceFloat(scRef->etaWidth(),nBits_); 
          patElectron_scEoP[iEle] = reduceFloat(scRef->energy()/iElectron.trackMomentumAtVtx().R(),nBits_); 
          patElectron_scSwissCross[iEle] = reduceFloat(swissCross,nBits_); 
          patElectron_scR9[iEle] = reduceFloat(e3x3/scRef->rawEnergy(),nBits_); 
          patElectron_scEMax[iEle] = reduceFloat(eMax,nBits_); 
          patElectron_scSigmaIEtaIEta[iEle] = reduceFloat(eleSS.sigmaIetaIeta,nBits_);
          patElectron_scSigmaIEtaIPhi[iEle] = reduceFloat(eleSS.sigmaIetaIphi,nBits_);
          patElectron_scSigmaIPhiIPhi[iEle] = reduceFloat(eleSS.sigmaIphiIphi,nBits_);  
          patElectron_full5x5_scR9[iEle] = reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_);
          patElectron_full5x5_scEMax[iEle] = reduceFloat(full5x5_eMax,nBits_);
          patElectron_full5x5_scSigmaIEtaIEta[iEle] = reduceFloat(full5x5_eleSS.sigmaIetaIeta,nBits_);
          patElectron_full5x5_scSigmaIEtaIPhi[iEle] = reduceFloat(full5x5_eleSS.sigmaIetaIphi,nBits_);
          patElectron_full5x5_scSigmaIPhiIPhi[iEle] = reduceFloat(full5x5_eleSS.sigmaIphiIphi,nBits_);  
          patElectron_egmCutBasedElectronIDVeto[iEle] = iElectron.electronID(egmCutBasedElectronIDVeto_.c_str());
          patElectron_egmCutBasedElectronIDloose[iEle] = iElectron.electronID(egmCutBasedElectronIDloose_.c_str());
          patElectron_egmCutBasedElectronIDmedium[iEle] = iElectron.electronID(egmCutBasedElectronIDmedium_.c_str());
          patElectron_egmCutBasedElectronIDtight[iEle] = iElectron.electronID(egmCutBasedElectronIDtight_.c_str());
          patElectron_egmMVAElectronIDloose[iEle] = iElectron.electronID(egmMVAElectronIDloose_.c_str());
          patElectron_egmMVAElectronIDmedium[iEle] = iElectron.electronID(egmMVAElectronIDmedium_.c_str());
          patElectron_egmMVAElectronIDtight[iEle] = iElectron.electronID(egmMVAElectronIDtight_.c_str());
          patElectron_egmMVAElectronIDlooseNoIso[iEle] = iElectron.electronID(egmMVAElectronIDlooseNoIso_.c_str());
          patElectron_egmMVAElectronIDmediumNoIso[iEle] = iElectron.electronID(egmMVAElectronIDmediumNoIso_.c_str());
          patElectron_egmMVAElectronIDtightNoIso[iEle] = iElectron.electronID(egmMVAElectronIDtightNoIso_.c_str());
          patElectron_heepElectronID[iEle] = iElectron.electronID(heepElectronID_.c_str());
          
          iEle++; 
      }
   }    
   
   //save photon info
   if(savePhotons_)
   {    
      int iPho=0;
      for(const auto& iPhoton : *(patPhoton.product())){ 
 
          reco::SuperClusterRef scRef = iPhoton.superCluster();
          reco::PhotonCoreRef phoCoreRef = iPhoton.photonCore();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

          double swissCross = -999.;
          double e3x3 = -999.;
          //double e2x2 = -999.;
          //double e5x5 = -999.;  
          double eMax = -999.;
          double full5x5_e3x3 = -999.;
          //double full5x5_e2x2 = -999.;
          //double full5x5_e5x5 = -999.;  
          double full5x5_eMax = -999.; 
         
          reco::Photon::ShowerShape phoSS = iPhoton.showerShapeVariables(); 
          reco::Photon::ShowerShape full5x5_phoSS = iPhoton.full5x5_showerShapeVariables();    
          const std::vector<std::pair<DetId,float> > &hits= iPhoton.superCluster()->hitsAndFractions();
          if(iPhoton.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             //e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             //e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
             //full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             //full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iPhoton.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             //e2x2 = EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             //e5x5 = EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             //full5x5_e2x2 = noZS::EcalClusterTools::e2x2( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             //full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
          }

          patPhoton_index[iPho] = iPho;
          patPhoton_scNPFClusters[iPho] = scRef->clusters().size();
          patPhoton_isEB[iPho] = iPhoton.isEB();
          patPhoton_isEE[iPho] = iPhoton.isEE();  
          patPhoton_isEBEEGap[iPho] = iPhoton.isEBEEGap();  
          patPhoton_isEBEtaGap[iPho] = iPhoton.isEBEtaGap();  
          patPhoton_isEBPhiGap[iPho] = iPhoton.isEBPhiGap();  
          patPhoton_isEEDeeGap[iPho] = iPhoton.isEEDeeGap();  
          patPhoton_isEERingGap[iPho] = iPhoton.isEERingGap();   
          patPhoton_passElectronVeto[iPho] = iPhoton.passElectronVeto(); 
          patPhoton_hasPixelSeed[iPho] = iPhoton.hasPixelSeed();  
          patPhoton_hasConversionTracks[iPho] = iPhoton.hasConversionTracks();
          patPhoton_nConversions[iPho] = phoCoreRef->conversions().size();
          patPhoton_nConversionsOneLeg[iPho] = phoCoreRef->conversionsOneLeg().size();
          patPhoton_eta[iPho] = reduceFloat(iPhoton.p4().eta(),nBits_);
          patPhoton_phi[iPho] = reduceFloat(iPhoton.p4().phi(),nBits_);
          patPhoton_energy[iPho] = reduceFloat(iPhoton.energy(),nBits_);
          patPhoton_energyErr[iPho] = reduceFloat(iPhoton.getCorrectedEnergyError(reco::Photon::regression2),nBits_);
          patPhoton_ecalEnergy[iPho] = reduceFloat(iPhoton.energyCorrections().phoEcalEnergy,nBits_);
          patPhoton_ecalEnergyErr[iPho] = reduceFloat(iPhoton.energyCorrections().phoEcalEnergyError,nBits_); 
          patPhoton_et[iPho] = reduceFloat(iPhoton.p4().Et(),nBits_);
          patPhoton_HoE[iPho] = reduceFloat(iPhoton.hadronicOverEm(),nBits_); 
          patPhoton_scEta[iPho] = reduceFloat(scRef->eta(),nBits_);
          patPhoton_scPhi[iPho] = reduceFloat(scRef->phi(),nBits_);
          patPhoton_scEnergy[iPho] = reduceFloat(scRef->energy(),nBits_); 
          patPhoton_scRawEnergy[iPho] = reduceFloat(scRef->rawEnergy(),nBits_); 
          patPhoton_scRawESEnergy[iPho] = reduceFloat(scRef->preshowerEnergy(),nBits_); 
          patPhoton_scEt[iPho] = reduceFloat(scRef->energy()*(Rt/R),nBits_);
          patPhoton_scPhiWidth[iPho] = reduceFloat(scRef->phiWidth(),nBits_);
          patPhoton_scEtaWidth[iPho] = reduceFloat(scRef->etaWidth(),nBits_); 
          patPhoton_scSwissCross[iPho] = reduceFloat(swissCross,nBits_); 
          patPhoton_scR9[iPho] = reduceFloat(e3x3/scRef->rawEnergy(),nBits_); 
          patPhoton_scEMax[iPho] = reduceFloat(eMax,nBits_); 
          patPhoton_scSigmaIEtaIEta[iPho] = reduceFloat(phoSS.sigmaIetaIeta,nBits_);
          patPhoton_scSigmaIEtaIPhi[iPho] = reduceFloat(phoSS.sigmaIetaIphi,nBits_);
          patPhoton_scSigmaIPhiIPhi[iPho] = reduceFloat(phoSS.sigmaIphiIphi,nBits_);   
          patPhoton_full5x5_scR9[iPho] = reduceFloat(full5x5_e3x3/scRef->rawEnergy(),nBits_);  
          patPhoton_full5x5_scEMax[iPho] = reduceFloat(full5x5_eMax,nBits_); 
          patPhoton_full5x5_scSigmaIEtaIEta[iPho] = reduceFloat(full5x5_phoSS.sigmaIetaIeta,nBits_);
          patPhoton_full5x5_scSigmaIEtaIPhi[iPho] = reduceFloat(full5x5_phoSS.sigmaIetaIphi,nBits_);
          patPhoton_full5x5_scSigmaIPhiIPhi[iPho] = reduceFloat(full5x5_phoSS.sigmaIphiIphi,nBits_);  
          patPhoton_egmCutBasedPhotonIDloose[iPho] = iPhoton.photonID(egmCutBasedPhotonIDloose_.c_str());
          patPhoton_egmCutBasedPhotonIDmedium[iPho] = iPhoton.photonID(egmCutBasedPhotonIDmedium_.c_str());
          patPhoton_egmCutBasedPhotonIDtight[iPho] = iPhoton.photonID(egmCutBasedPhotonIDtight_.c_str());
          patPhoton_egmMVAPhotonIDmedium[iPho] = iPhoton.photonID(egmMVAPhotonIDmedium_.c_str());
          patPhoton_egmMVAPhotonIDtight[iPho] = iPhoton.photonID(egmMVAPhotonIDtight_.c_str());    
          patPhoton_passPreselections[iPho] = passPreselections(&iPhoton);
          
          iPho++; 
      }
   }
   
   //fill tree for each event
   tree->Fill();
}

void ElePhoDumper::beginJob()
{

}

void ElePhoDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void ElePhoDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("rho", &rho, "rho/F"); 
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   if(isMC_){ 
      tree->Branch("truePU", &truePU, "truePU/D");
      tree->Branch("obsPU", &obsPU, "obsPU/D");
      tree->Branch("genParticle_size", &genParticle_size, "genParticle_size/I");  
      tree->Branch("genParticle_pdgId","std::vector<int>",&genParticle_pdgId);
      tree->Branch("genParticle_status","std::vector<int>",&genParticle_status); 
      tree->Branch("genParticle_statusFlag","std::vector<int>",&genParticle_statusFlag); 
      tree->Branch("genParticle_statusFlag_isPrompt","std::vector<int>",&genParticle_statusFlag_isPrompt); 
      tree->Branch("genParticle_statusFlag_isLastCopyBeforeFSR","std::vector<int>",&genParticle_statusFlag_isLastCopyBeforeFSR); 
      tree->Branch("genParticle_statusFlag_isLastCopy","std::vector<int>",&genParticle_statusFlag_isLastCopy); 
      tree->Branch("genParticle_statusFlag_isFirstCopy","std::vector<int>",&genParticle_statusFlag_isFirstCopy); 
      tree->Branch("genParticle_statusFlag_fromHardProcessBeforeFSR","std::vector<int>",&genParticle_statusFlag_fromHardProcessBeforeFSR); 
      tree->Branch("genParticle_statusFlag_isDirectHardProcessTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isDirectHardProcessTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_isHardProcessTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isHardProcessTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_fromHardProcess","std::vector<int>",&genParticle_statusFlag_fromHardProcess); 
      tree->Branch("genParticle_statusFlag_isHardProcess","std::vector<int>",&genParticle_statusFlag_isHardProcess); 
      tree->Branch("genParticle_statusFlag_isDirectHadronDecayProduct","std::vector<int>",&genParticle_statusFlag_isDirectHadronDecayProduct); 
      tree->Branch("genParticle_statusFlag_isDirectPromptTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isDirectPromptTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_isDirectTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isDirectTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_isPromptTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isPromptTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_isTauDecayProduct","std::vector<int>",&genParticle_statusFlag_isTauDecayProduct); 
      tree->Branch("genParticle_statusFlag_isDecayedLeptonHadron","std::vector<int>",&genParticle_statusFlag_isDecayedLeptonHadron); 
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
   }
   if(saveElectrons_){
      tree->Branch("patElectron_size", &patElectron_size, "patElectron_size/I");  
      tree->Branch("patElectron_index","std::vector<int>",&patElectron_index); 
      tree->Branch("patElectron_classification","std::vector<int>",&patElectron_classification); 
      tree->Branch("patElectron_scNPFClusters","std::vector<int>",&patElectron_scNPFClusters); 
      tree->Branch("patElectron_isEB","std::vector<bool> ",&patElectron_isEB); 
      tree->Branch("patElectron_isEE","std::vector<bool> ",&patElectron_isEE); 
      tree->Branch("patElectron_isEBEEGap","std::vector<bool> ",&patElectron_isEBEEGap);
      tree->Branch("patElectron_isEBEtaGap","std::vector<bool> ",&patElectron_isEBEtaGap);
      tree->Branch("patElectron_isEBPhiGap","std::vector<bool> ",&patElectron_isEBPhiGap);
      tree->Branch("patElectron_isEEDeeGap","std::vector<bool> ",&patElectron_isEEDeeGap);
      tree->Branch("patElectron_isEERingGap","std::vector<bool> ",&patElectron_isEERingGap);
      tree->Branch("patElectron_isEcalDriven","std::vector<bool> ",&patElectron_isEcalDriven);
      tree->Branch("patElectron_isTrackerDriven","std::vector<bool> ",&patElectron_isTrackerDriven);
      tree->Branch("patElectron_passConversionVeto","std::vector<bool> ",&patElectron_passConversionVeto);
      tree->Branch("patElectron_eta","std::vector<float>",&patElectron_eta); 
      tree->Branch("patElectron_phi","std::vector<float>",&patElectron_phi); 
      tree->Branch("patElectron_p","std::vector<float>",&patElectron_p); 
      tree->Branch("patElectron_pt","std::vector<float>",&patElectron_pt); 
      tree->Branch("patElectron_pIn","std::vector<float>",&patElectron_pIn); 
      tree->Branch("patElectron_pOut","std::vector<float>",&patElectron_pOut); 
      tree->Branch("patElectron_trackFbrem","std::vector<float>",&patElectron_trackFbrem);
      tree->Branch("patElectron_superClusterFbrem","std::vector<float>",&patElectron_superClusterFbrem);
      tree->Branch("patElectron_energy","std::vector<float>",&patElectron_energy);   
      tree->Branch("patElectron_energyErr","std::vector<float>",&patElectron_energyErr);   
      tree->Branch("patElectron_ecalEnergy","std::vector<float>",&patElectron_ecalEnergy); 
      tree->Branch("patElectron_ecalEnergyErr","std::vector<float>",&patElectron_ecalEnergyErr);  
      tree->Branch("patElectron_et","std::vector<float>",&patElectron_et); 
      tree->Branch("patElectron_HoE","std::vector<float>",&patElectron_HoE);    
      tree->Branch("patElectron_scEta","std::vector<float>",&patElectron_scEta); 
      tree->Branch("patElectron_scPhi","std::vector<float>",&patElectron_scPhi);
      tree->Branch("patElectron_scEnergy","std::vector<float>",&patElectron_scEnergy); 
      tree->Branch("patElectron_scRawEnergy","std::vector<float>",&patElectron_scRawEnergy); 
      tree->Branch("patElectron_scRawESEnergy","std::vector<float>",&patElectron_scRawESEnergy); 
      tree->Branch("patElectron_scEt","std::vector<float>",&patElectron_scEt); 
      tree->Branch("patElectron_scPhiWidth","std::vector<float>",&patElectron_scPhiWidth);  
      tree->Branch("patElectron_scEtaWidth","std::vector<float>",&patElectron_scEtaWidth);  
      tree->Branch("patElectron_scEoP","std::vector<float>",&patElectron_scEoP);  
      tree->Branch("patElectron_scSwissCross","std::vector<float>",&patElectron_scSwissCross);  
      tree->Branch("patElectron_scR9","std::vector<float>",&patElectron_scR9); 
      tree->Branch("patElectron_scEMax","std::vector<float>",&patElectron_scEMax); 
      tree->Branch("patElectron_scSigmaIEtaIEta","std::vector<float>",&patElectron_scSigmaIEtaIEta);  
      tree->Branch("patElectron_scSigmaIEtaIPhi","std::vector<float>",&patElectron_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_scSigmaIPhiIPhi","std::vector<float>",&patElectron_scSigmaIPhiIPhi);  
      tree->Branch("patElectron_full5x5_scR9","std::vector<float>",&patElectron_full5x5_scR9);  
      tree->Branch("patElectron_full5x5_scEMax","std::vector<float>",&patElectron_full5x5_scEMax);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIEta","std::vector<float>",&patElectron_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIPhi","std::vector<float>",&patElectron_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_full5x5_scSigmaIPhiIPhi","std::vector<float>",&patElectron_full5x5_scSigmaIPhiIPhi);  
      tree->Branch("patElectron_egmCutBasedElectronIDVeto","std::vector<int>",&patElectron_egmCutBasedElectronIDVeto);  
      tree->Branch("patElectron_egmCutBasedElectronIDloose","std::vector<int>",&patElectron_egmCutBasedElectronIDloose);  
      tree->Branch("patElectron_egmCutBasedElectronIDmedium","std::vector<int>",&patElectron_egmCutBasedElectronIDmedium);
      tree->Branch("patElectron_egmCutBasedElectronIDtight","std::vector<int>",&patElectron_egmCutBasedElectronIDtight);  
      tree->Branch("patElectron_egmMVAElectronIDloose","std::vector<int>",&patElectron_egmMVAElectronIDloose);  
      tree->Branch("patElectron_egmMVAElectronIDmedium","std::vector<int>",&patElectron_egmMVAElectronIDmedium);
      tree->Branch("patElectron_egmMVAElectronIDtight","std::vector<int>",&patElectron_egmMVAElectronIDtight);  
      tree->Branch("patElectron_egmMVAElectronIDlooseNoIso","std::vector<int>",&patElectron_egmMVAElectronIDlooseNoIso);  
      tree->Branch("patElectron_egmMVAElectronIDmediumNoIso","std::vector<int>",&patElectron_egmMVAElectronIDmediumNoIso);
      tree->Branch("patElectron_egmMVAElectronIDtightNoIso","std::vector<int>",&patElectron_egmMVAElectronIDtightNoIso);  
      tree->Branch("patElectron_heepElectronID","std::vector<int>",&patElectron_heepElectronID);  
   }
   if(savePhotons_){
      tree->Branch("patPhoton_size", &patPhoton_size, "patPhoton_size/I");  
      tree->Branch("patPhoton_index","std::vector<int>",&patPhoton_index); 
      tree->Branch("patPhoton_scNPFClusters","std::vector<int>",&patPhoton_scNPFClusters);   
      tree->Branch("patPhoton_passElectronVeto","std::vector<bool> ",&patPhoton_passElectronVeto); 
      tree->Branch("patPhoton_hasPixelSeed","std::vector<bool> ",&patPhoton_hasPixelSeed); 
      tree->Branch("patPhoton_hasConversionTracks","std::vector<bool> ",&patPhoton_hasConversionTracks); 
      tree->Branch("patPhoton_nConversions","std::vector<int>",&patPhoton_nConversions);     
      tree->Branch("patPhoton_nConversionsOneLeg","std::vector<int>",&patPhoton_nConversionsOneLeg);     
      tree->Branch("patPhoton_isEB","std::vector<bool> ",&patPhoton_isEB); 
      tree->Branch("patPhoton_isEE","std::vector<bool> ",&patPhoton_isEE); 
      tree->Branch("patPhoton_isEBEEGap","std::vector<bool> ",&patPhoton_isEBEEGap);
      tree->Branch("patPhoton_isEBEtaGap","std::vector<bool> ",&patPhoton_isEBEtaGap);
      tree->Branch("patPhoton_isEBPhiGap","std::vector<bool> ",&patPhoton_isEBPhiGap);
      tree->Branch("patPhoton_isEEDeeGap","std::vector<bool> ",&patPhoton_isEEDeeGap);
      tree->Branch("patPhoton_isEERingGap","std::vector<bool> ",&patPhoton_isEERingGap);  
      tree->Branch("patPhoton_eta","std::vector<float>",&patPhoton_eta); 
      tree->Branch("patPhoton_phi","std::vector<float>",&patPhoton_phi); 
      tree->Branch("patPhoton_energy","std::vector<float>",&patPhoton_energy);   
      tree->Branch("patPhoton_energyErr","std::vector<float>",&patPhoton_energyErr);   
      tree->Branch("patPhoton_ecalEnergy","std::vector<float>",&patPhoton_ecalEnergy); 
      tree->Branch("patPhoton_ecalEnergyErr","std::vector<float>",&patPhoton_ecalEnergyErr);  
      tree->Branch("patPhoton_et","std::vector<float>",&patPhoton_et); 
      tree->Branch("patPhoton_HoE","std::vector<float>",&patPhoton_HoE);  
      tree->Branch("patPhoton_scEta","std::vector<float>",&patPhoton_scEta); 
      tree->Branch("patPhoton_scPhi","std::vector<float>",&patPhoton_scPhi);  
      tree->Branch("patPhoton_scEnergy","std::vector<float>",&patPhoton_scEnergy); 
      tree->Branch("patPhoton_scRawEnergy","std::vector<float>",&patPhoton_scRawEnergy); 
      tree->Branch("patPhoton_scRawESEnergy","std::vector<float>",&patPhoton_scRawESEnergy); 
      tree->Branch("patPhoton_scEt","std::vector<float>",&patPhoton_scEt); 
      tree->Branch("patPhoton_scPhiWidth","std::vector<float>",&patPhoton_scPhiWidth);  
      tree->Branch("patPhoton_scEtaWidth","std::vector<float>",&patPhoton_scEtaWidth);  
      tree->Branch("patPhoton_scSwissCross","std::vector<float>",&patPhoton_scSwissCross);  
      tree->Branch("patPhoton_scR9","std::vector<float>",&patPhoton_scR9); 
      tree->Branch("patPhoton_scEMax","std::vector<float>",&patPhoton_scEMax); 
      tree->Branch("patPhoton_scSigmaIEtaIEta","std::vector<float>",&patPhoton_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_scSigmaIEtaIPhi","std::vector<float>",&patPhoton_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_scSigmaIPhiIPhi","std::vector<float>",&patPhoton_scSigmaIPhiIPhi);  
      tree->Branch("patPhoton_full5x5_scR9","std::vector<float>",&patPhoton_full5x5_scR9);  
      tree->Branch("patPhoton_full5x5_scEMax","std::vector<float>",&patPhoton_full5x5_scEMax);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIEta","std::vector<float>",&patPhoton_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIPhi","std::vector<float>",&patPhoton_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_full5x5_scSigmaIPhiIPhi","std::vector<float>",&patPhoton_full5x5_scSigmaIPhiIPhi);  
      tree->Branch("patPhoton_egmCutBasedPhotonIDloose","std::vector<int>",&patPhoton_egmCutBasedPhotonIDloose);  
      tree->Branch("patPhoton_egmCutBasedPhotonIDmedium","std::vector<int>",&patPhoton_egmCutBasedPhotonIDmedium);
      tree->Branch("patPhoton_egmCutBasedPhotonIDtight","std::vector<int>",&patPhoton_egmCutBasedPhotonIDtight);   
      tree->Branch("patPhoton_egmMVAPhotonIDmedium","std::vector<int>",&patPhoton_egmMVAPhotonIDmedium);
      tree->Branch("patPhoton_egmMVAPhotonIDtight","std::vector<int>",&patPhoton_egmMVAPhotonIDtight);  
      tree->Branch("patPhoton_passPreselections","std::vector<int>",&patPhoton_passPreselections);
   }
}

void ElePhoDumper::setVectors(int nElectrons, int nPhotons)
{
   
   genParticle_pdgId.clear();
   genParticle_status.clear();
   genParticle_statusFlag.clear();
   genParticle_statusFlag_isLastCopyBeforeFSR.clear();              
   genParticle_statusFlag_isLastCopy.clear();  
   genParticle_statusFlag_isFirstCopy.clear();  
   genParticle_statusFlag_fromHardProcessBeforeFSR.clear();  
   genParticle_statusFlag_isDirectHardProcessTauDecayProduct.clear();    
   genParticle_statusFlag_isHardProcessTauDecayProduct.clear();   
   genParticle_statusFlag_fromHardProcess.clear();
   genParticle_statusFlag_isHardProcess.clear();  
   genParticle_statusFlag_isDirectHadronDecayProduct.clear();   
   genParticle_statusFlag_isDirectPromptTauDecayProduct.clear();  
   genParticle_statusFlag_isDirectTauDecayProduct.clear();   
   genParticle_statusFlag_isPromptTauDecayProduct.clear();   
   genParticle_statusFlag_isTauDecayProduct.clear();   
   genParticle_statusFlag_isDecayedLeptonHadron.clear(); 
   genParticle_statusFlag_isPrompt.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
    
   patElectron_index.clear();
   patElectron_classification.clear();
   patElectron_scNPFClusters.clear();
   patElectron_charge.clear();
   patElectron_isEB.clear();
   patElectron_isEE.clear();  
   patElectron_isEBEEGap.clear();
   patElectron_isEBEtaGap.clear();
   patElectron_isEBPhiGap.clear();
   patElectron_isEEDeeGap.clear();
   patElectron_isEERingGap.clear(); 
   patElectron_isEcalDriven.clear();
   patElectron_isTrackerDriven.clear();
   patElectron_passConversionVeto.clear();
   patElectron_eta.clear();
   patElectron_phi.clear();
   patElectron_p.clear();
   patElectron_pt.clear();
   patElectron_pIn.clear();
   patElectron_pOut.clear();
   patElectron_trackFbrem.clear();
   patElectron_superClusterFbrem.clear();
   patElectron_energy.clear();
   patElectron_energyErr.clear();
   patElectron_ecalEnergy.clear();
   patElectron_ecalEnergyErr.clear();
   patElectron_et.clear();
   patElectron_HoE.clear();
   patElectron_scEta.clear();
   patElectron_scPhi.clear();
   patElectron_scEnergy.clear();
   patElectron_scRawEnergy.clear();
   patElectron_scRawESEnergy.clear();
   patElectron_scEt.clear();
   patElectron_scPhiWidth.clear();
   patElectron_scEtaWidth.clear();
   patElectron_scEoP.clear();
   patElectron_scSwissCross.clear();
   patElectron_scR9.clear();
   patElectron_scEMax.clear();
   patElectron_scSigmaIEtaIEta.clear();
   patElectron_scSigmaIEtaIPhi.clear();
   patElectron_scSigmaIPhiIPhi.clear();
   patElectron_full5x5_scR9.clear();
   patElectron_full5x5_scEMax.clear();
   patElectron_full5x5_scSigmaIEtaIEta.clear();
   patElectron_full5x5_scSigmaIEtaIPhi.clear();
   patElectron_full5x5_scSigmaIPhiIPhi.clear();
   patElectron_egmCutBasedElectronIDVeto.clear();
   patElectron_egmCutBasedElectronIDloose.clear();
   patElectron_egmCutBasedElectronIDmedium.clear();
   patElectron_egmCutBasedElectronIDtight.clear();
   patElectron_egmMVAElectronIDloose.clear();
   patElectron_egmMVAElectronIDmedium.clear();
   patElectron_egmMVAElectronIDtight.clear();
   patElectron_egmMVAElectronIDlooseNoIso.clear();
   patElectron_egmMVAElectronIDmediumNoIso.clear();
   patElectron_egmMVAElectronIDtightNoIso.clear();
   patElectron_heepElectronID.clear();
   if(saveElectrons_){
      patElectron_index.resize(nElectrons);
      patElectron_classification.resize(nElectrons);
      patElectron_scNPFClusters.resize(nElectrons);
      patElectron_charge.resize(nElectrons);
      patElectron_isEB.resize(nElectrons);
      patElectron_isEE.resize(nElectrons);  
      patElectron_isEBEEGap.resize(nElectrons);
      patElectron_isEBEtaGap.resize(nElectrons);
      patElectron_isEBPhiGap.resize(nElectrons);
      patElectron_isEEDeeGap.resize(nElectrons);
      patElectron_isEERingGap.resize(nElectrons); 
      patElectron_isEcalDriven.resize(nElectrons);
      patElectron_isTrackerDriven.resize(nElectrons);
      patElectron_passConversionVeto.resize(nElectrons);
      patElectron_eta.resize(nElectrons);
      patElectron_phi.resize(nElectrons);
      patElectron_p.resize(nElectrons);
      patElectron_pt.resize(nElectrons);
      patElectron_pIn.resize(nElectrons);
      patElectron_pOut.resize(nElectrons);
      patElectron_trackFbrem.resize(nElectrons);
      patElectron_superClusterFbrem.resize(nElectrons);
      patElectron_energy.resize(nElectrons);
      patElectron_energyErr.resize(nElectrons);
      patElectron_ecalEnergy.resize(nElectrons);
      patElectron_ecalEnergyErr.resize(nElectrons);
      patElectron_et.resize(nElectrons);
      patElectron_HoE.resize(nElectrons);
      patElectron_scEta.resize(nElectrons);
      patElectron_scPhi.resize(nElectrons);
      patElectron_scEnergy.resize(nElectrons);
      patElectron_scRawEnergy.resize(nElectrons);
      patElectron_scRawESEnergy.resize(nElectrons);
      patElectron_scEt.resize(nElectrons);
      patElectron_scPhiWidth.resize(nElectrons);
      patElectron_scEtaWidth.resize(nElectrons);
      patElectron_scEoP.resize(nElectrons);
      patElectron_scSwissCross.resize(nElectrons);
      patElectron_scR9.resize(nElectrons);
      patElectron_scEMax.resize(nElectrons);
      patElectron_scSigmaIEtaIEta.resize(nElectrons);
      patElectron_scSigmaIEtaIPhi.resize(nElectrons);
      patElectron_scSigmaIPhiIPhi.resize(nElectrons);
      patElectron_full5x5_scR9.resize(nElectrons);
      patElectron_full5x5_scEMax.resize(nElectrons);
      patElectron_full5x5_scSigmaIEtaIEta.resize(nElectrons);
      patElectron_full5x5_scSigmaIEtaIPhi.resize(nElectrons);
      patElectron_full5x5_scSigmaIPhiIPhi.resize(nElectrons);
      patElectron_egmCutBasedElectronIDVeto.resize(nElectrons);
      patElectron_egmCutBasedElectronIDloose.resize(nElectrons);
      patElectron_egmCutBasedElectronIDmedium.resize(nElectrons);
      patElectron_egmCutBasedElectronIDtight.resize(nElectrons);
      patElectron_egmMVAElectronIDloose.resize(nElectrons);
      patElectron_egmMVAElectronIDmedium.resize(nElectrons);
      patElectron_egmMVAElectronIDtight.resize(nElectrons);
      patElectron_egmMVAElectronIDlooseNoIso.resize(nElectrons);
      patElectron_egmMVAElectronIDmediumNoIso.resize(nElectrons);
      patElectron_egmMVAElectronIDtightNoIso.resize(nElectrons);
      patElectron_heepElectronID.resize(nElectrons);
   }

   patPhoton_index.clear();
   patPhoton_scNPFClusters.clear();
   patPhoton_isEB.clear();
   patPhoton_isEE.clear();
   patPhoton_isEBEEGap.clear(); 
   patPhoton_isEBEtaGap.clear();
   patPhoton_isEBPhiGap.clear();
   patPhoton_isEEDeeGap.clear();
   patPhoton_isEERingGap.clear(); 
   patPhoton_passElectronVeto.clear();
   patPhoton_hasPixelSeed.clear(); 
   patPhoton_hasConversionTracks.clear();
   patPhoton_nConversions.clear();
   patPhoton_nConversionsOneLeg.clear();
   patPhoton_eta.clear();
   patPhoton_phi.clear();
   patPhoton_energy.clear();
   patPhoton_energyErr.clear();
   patPhoton_ecalEnergy.clear();
   patPhoton_ecalEnergyErr.clear();
   patPhoton_et.clear();
   patPhoton_HoE.clear();
   patPhoton_scEta.clear();
   patPhoton_scPhi.clear();
   patPhoton_scEnergy.clear();
   patPhoton_scRawEnergy.clear();
   patPhoton_scRawESEnergy.clear();
   patPhoton_scEt.clear();
   patPhoton_scPhiWidth.clear();
   patPhoton_scEtaWidth.clear();
   patPhoton_scSwissCross.clear();
   patPhoton_scR9.clear();
   patPhoton_scEMax.clear();
   patPhoton_scSigmaIEtaIEta.clear();
   patPhoton_scSigmaIEtaIPhi.clear();
   patPhoton_scSigmaIPhiIPhi.clear(); 
   patPhoton_full5x5_scR9.clear();
   patPhoton_full5x5_scEMax.clear();
   patPhoton_full5x5_scSigmaIEtaIEta.clear();
   patPhoton_full5x5_scSigmaIEtaIPhi.clear();
   patPhoton_full5x5_scSigmaIPhiIPhi.clear();
   patPhoton_egmCutBasedPhotonIDloose.clear();
   patPhoton_egmCutBasedPhotonIDmedium.clear();
   patPhoton_egmCutBasedPhotonIDtight.clear();
   patPhoton_egmMVAPhotonIDmedium.clear();
   patPhoton_egmMVAPhotonIDtight.clear();
   patPhoton_passPreselections.clear();
   if(savePhotons_){
      patPhoton_index.resize(nPhotons);
      patPhoton_scNPFClusters.resize(nPhotons);
      patPhoton_isEB.resize(nPhotons);
      patPhoton_isEE.resize(nPhotons);
      patPhoton_isEBEEGap.resize(nPhotons); 
      patPhoton_isEBEtaGap.resize(nPhotons);
      patPhoton_isEBPhiGap.resize(nPhotons);
      patPhoton_isEEDeeGap.resize(nPhotons);
      patPhoton_isEERingGap.resize(nPhotons); 
      patPhoton_passElectronVeto.resize(nPhotons);
      patPhoton_hasPixelSeed.resize(nPhotons); 
      patPhoton_hasConversionTracks.resize(nPhotons);
      patPhoton_nConversions.resize(nPhotons);
      patPhoton_nConversionsOneLeg.resize(nPhotons);
      patPhoton_eta.resize(nPhotons);
      patPhoton_phi.resize(nPhotons);
      patPhoton_energy.resize(nPhotons);
      patPhoton_energyErr.resize(nPhotons);
      patPhoton_ecalEnergy.resize(nPhotons);
      patPhoton_ecalEnergyErr.resize(nPhotons);
      patPhoton_et.resize(nPhotons);
      patPhoton_HoE.resize(nPhotons);
      patPhoton_scEta.resize(nPhotons);
      patPhoton_scPhi.resize(nPhotons);
      patPhoton_scEnergy.resize(nPhotons);
      patPhoton_scRawEnergy.resize(nPhotons);
      patPhoton_scRawESEnergy.resize(nPhotons);
      patPhoton_scEt.resize(nPhotons);
      patPhoton_scPhiWidth.resize(nPhotons);
      patPhoton_scEtaWidth.resize(nPhotons);
      patPhoton_scSwissCross.resize(nPhotons);
      patPhoton_scR9.resize(nPhotons);
      patPhoton_scEMax.resize(nPhotons);
      patPhoton_scSigmaIEtaIEta.resize(nPhotons);
      patPhoton_scSigmaIEtaIPhi.resize(nPhotons);
      patPhoton_scSigmaIPhiIPhi.resize(nPhotons); 
      patPhoton_full5x5_scR9.resize(nPhotons);
      patPhoton_full5x5_scEMax.resize(nPhotons);
      patPhoton_full5x5_scSigmaIEtaIEta.resize(nPhotons);
      patPhoton_full5x5_scSigmaIEtaIPhi.resize(nPhotons);
      patPhoton_full5x5_scSigmaIPhiIPhi.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDloose.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDmedium.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDtight.resize(nPhotons);
      patPhoton_egmMVAPhotonIDmedium.resize(nPhotons);
      patPhoton_egmMVAPhotonIDtight.resize(nPhotons);
      patPhoton_passPreselections.resize(nPhotons);
   }
   
}

int ElePhoDumper::getGenStatusFlag(const reco::GenParticle* genParticle)
{
   int statusFlag = 
   genParticle->statusFlags().isLastCopyBeforeFSR()                  * 16384 +
   genParticle->statusFlags().isLastCopy()                           * 8192  +
   genParticle->statusFlags().isFirstCopy()                          * 4096  +
   genParticle->statusFlags().fromHardProcessBeforeFSR()             * 2048  +
   genParticle->statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +
   genParticle->statusFlags().isHardProcessTauDecayProduct()         * 512   +
   genParticle->statusFlags().fromHardProcess()                      * 256   +  
   genParticle->statusFlags().isHardProcess()                        * 128   +
   genParticle->statusFlags().isDirectHadronDecayProduct()           * 64    +
   genParticle->statusFlags().isDirectPromptTauDecayProduct()        * 32    +
   genParticle->statusFlags().isDirectTauDecayProduct()              * 16    +
   genParticle->statusFlags().isPromptTauDecayProduct()              * 8     +
   genParticle->statusFlags().isTauDecayProduct()                    * 4     +
   genParticle->statusFlags().isDecayedLeptonHadron()                * 2     +
   genParticle->statusFlags().isPrompt()                             * 1      ;
   
   return statusFlag; 
}

void ElePhoDumper::printGenStatusFlag(const reco::GenParticle* genParticle)
{
   std::cout << " -- Flags: "; 
   if(genParticle->statusFlags().isLastCopyBeforeFSR()) std::cout << "isLastCopyBeforeFSR() && ";
   if(genParticle->statusFlags().isLastCopy()) std::cout << "isLastCopy() && ";     
   if(genParticle->statusFlags().isFirstCopy()) std::cout << "isFirstCopy() && ";
   if(genParticle->statusFlags().fromHardProcessBeforeFSR()) std::cout << "fromHardProcessBeforeFSR() && ";
   if(genParticle->statusFlags().isDirectHardProcessTauDecayProduct()) std::cout << "isDirectHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().isHardProcessTauDecayProduct()) std::cout << "isHardProcessTauDecayProduct() && ";
   if(genParticle->statusFlags().fromHardProcess()) std::cout << "fromHardProcess() && "; 
   if(genParticle->statusFlags().isHardProcess()) std::cout << "isHardProcess() && ";
   if(genParticle->statusFlags().isDirectHadronDecayProduct()) std::cout << "isDirectHadronDecayProduct() && ";
   if(genParticle->statusFlags().isDirectPromptTauDecayProduct()) std::cout << "isDirectPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isDirectTauDecayProduct()) std::cout << "isDirectTauDecayProduct() && ";
   if(genParticle->statusFlags().isPromptTauDecayProduct()) std::cout << "isPromptTauDecayProduct() && ";
   if(genParticle->statusFlags().isTauDecayProduct()) std::cout << "isTauDecayProduct() && ";
   if(genParticle->statusFlags().isDecayedLeptonHadron()) std::cout << "isDecayedLeptonHadron() && ";
   if(genParticle->statusFlags().isPrompt()) std::cout << "isPrompt() ";
   std::cout << "-- " << getGenStatusFlag(genParticle) << " -- ";
}

int ElePhoDumper::passPreselections(const pat::Photon* photon)
{
   int pass = 0;

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
   ) pass = 1;
   
   return pass;
}   


float ElePhoDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElePhoDumper);
