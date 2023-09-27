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


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"


#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TVirtualFitter.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
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
  // ----------member data ---------------------------
  TH1F* Profile;
  TF1* f;
  TTree* tree;

  long long int runId;
  long long int eventId;
  long long int lumiId;
  int genParticle_pdgId;
  int genParticle_status; 
  double genParticle_energy;
  double genParticle_pt;
  double genParticle_eta;
  double genParticle_phi; 
  double xtal_energyMax;
  double xtal_binMaxVal;
  double xtal_tMax;
  double xtal_tMaxError;
  double xtal_tMax_fitPar0;
  double xtal_tMax_fitPar1;
  double xtal_tMax_fitPar2;
  double xtal_tMax_fitParError0;
  double xtal_tMax_fitParError1;
  double xtal_tMax_fitParError2;
  int xtal_tMax_fitStatus;
  int xtal_ieta;
  int xtal_iphi;
  int xtal_iz;

  std::map<long long int,double> energyMapEB;
  std::map<long long int,double> energyMapEE;

  edm::Handle<std::vector<reco::GenParticle> > genParticle;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticleToken;

  edm::Handle<edm::PCaloHitContainer> pCaloHits_EB_Handle;
  edm::EDGetTokenT<edm::PCaloHitContainer> pCaloHits_EB_Token;
  
  edm::Handle<edm::PCaloHitContainer> pCaloHits_EE_Handle;
  edm::EDGetTokenT<edm::PCaloHitContainer> pCaloHits_EE_Token;  

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
	genParticleToken = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
        pCaloHits_EB_Token = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("pcaloHitsEB"));
        pCaloHits_EE_Token = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("pcaloHitsEE"));

        edm::Service<TFileService> fs;
        tree = fs->make<TTree>("caloTree","caloTree"); 
        setTree(tree);
        //Profile = fs->make<TH1F>("Profile","Profile",260,0., 26.); 
        Profile = new TH1F("Profile","Profile",260,0., 26.); 
        Profile->Sumw2(); 
        

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
void simAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  runId = iEvent.id().run();
  eventId = iEvent.id().event();
  lumiId = iEvent.id().luminosityBlock();
  
  iEvent.getByToken(genParticleToken, genParticle);         //GenParticles
  iEvent.getByToken(pCaloHits_EB_Token, pCaloHits_EB_Handle); //PCaloHits EB
  iEvent.getByToken(pCaloHits_EE_Token, pCaloHits_EE_Handle); //PCaloHits EE

  energyMapEB.clear();
  energyMapEE.clear();
  Profile->Reset();
  
  EBDetId* DidEB;
  EEDetId* DidEE;

  edm::PCaloHitContainer::const_iterator caloHitsItr;

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

  std::map<long long int,double>::iterator maxValEB = std::max_element(energyMapEB.begin(),energyMapEB.end(),[] (const std::pair<long long int,double>& a, const std::pair<long long int,double>& b)->bool{ return a.second < b.second; } );
  std::map<long long int,double>::iterator maxValEE = std::max_element(energyMapEE.begin(),energyMapEE.end(),[] (const std::pair<long long int,double>& a, const std::pair<long long int,double>& b)->bool{ return a.second < b.second; } );
  
  bool isEB = false;
  if(maxValEB->second > maxValEE->second) isEB = true;

  long long int rawId = 0;
  double energyMax = -999.;
  if(isEB){ 
     rawId = maxValEB->first;
     energyMax = maxValEB->second;
  }else{ 
     rawId = maxValEE->first;
     energyMax = maxValEE->second; 
  }
  xtal_energyMax = energyMax;

  if(isEB){
     for(caloHitsItr = pCaloHits_EB_Handle->begin(); caloHitsItr != pCaloHits_EB_Handle->end(); caloHitsItr++){
         DidEB = new EBDetId(caloHitsItr->id());
         if(DidEB->rawId() != rawId) continue;

         xtal_ieta = DidEB->ieta();
         xtal_iphi = DidEB->iphi();
         xtal_iz = 0;
 
         int temp = caloHitsItr->depth(); 
         if((temp == 1) || (temp == 2)) continue;
         double depth = (caloHitsItr->depth()>>3)/100.;
         Profile->Fill(depth, caloHitsItr->energy());
     }
  }else{
     for(caloHitsItr = pCaloHits_EE_Handle->begin(); caloHitsItr != pCaloHits_EE_Handle->end(); caloHitsItr++){
         DidEE = new EEDetId(caloHitsItr->id());
         if(DidEE->rawId() != rawId) continue;

         xtal_ieta = DidEE->ix();
         xtal_iphi = DidEE->iy();
         xtal_iz = DidEE->zside();;

         int temp = caloHitsItr->depth(); 
         if((temp == 1) || (temp == 2)) continue;
         double depth = (caloHitsItr->depth()>>3)/100.;
         Profile->Fill(depth, caloHitsItr->energy());
     } 
  }

  int binmax = Profile->GetMaximumBin(); 
  double x = Profile->GetXaxis()->GetBinCenter(binmax);
  double x_min = x-1*Profile->GetRMS();
  double x_max = x+1*Profile->GetRMS();

  f = new TF1("f","[0]*TMath::Power(x,[1]*[2])*TMath::Exp(-1*[2]*x)",x_min,x_max);
  f->SetLineColor(kRed+1);
  f->SetLineWidth(2);
  f->SetParameter(0,0.001);
  f->SetParLimits(0,0.,100.); 
  f->SetParameter(1,x);
  f->SetParLimits(1,0.,26.);  
  f->SetParameter(2,0.5); 
  f->SetParLimits(2,0.,100.);

  //TVirtualFitter::SetDefaultFitter("Minuit2");
  int status = 0;
  TFitResultPtr rp;  
  rp = Profile->Fit("f", "QMERS+");
  while(status != 0)
  {
      rp = Profile->Fit("f", "QMERS+");
      status = rp;
      if(status == 0) break;
  }
  
  x_min = x-2*Profile->GetRMS();
  x_max = x+2*Profile->GetRMS();
  f->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));
  f->SetParLimits(0,f->GetParameter(0)-0.5*f->GetParError(0),f->GetParameter(0)+0.5*f->GetParError(0)); 
  f->SetParLimits(1,f->GetParameter(1)-0.1*f->GetParError(1),f->GetParameter(1)+0.1*f->GetParError(1)); 
  f->SetParLimits(2,f->GetParameter(2)-0.5*f->GetParError(2),f->GetParameter(2)+0.5*f->GetParError(2)); 
  status = 0;
  rp = Profile->Fit("f", "QMERS+");
  while(status != 0)
  {
      rp = Profile->Fit("f", "QMERS+");
      status = rp;
      if(status == 0) break;
  }

  xtal_tMax = f->GetParameter(1);  
  xtal_tMaxError = f->GetParError(1);  
  //std::cout << "tMax:" << xtal_tMax << " - " << status << " - " << xtal_energyMax << std::endl;

  xtal_tMax_fitPar0 = f->GetParameter(0);
  xtal_tMax_fitPar1 = f->GetParameter(1);
  xtal_tMax_fitPar2 = f->GetParameter(2);
  xtal_tMax_fitParError0 = f->GetParError(0);
  xtal_tMax_fitParError1 = f->GetParError(1);
  xtal_tMax_fitParError2 = f->GetParError(2);
  xtal_tMax_fitStatus = status;

  TH1F* h = (TH1F*)Profile->Clone();
  h->Rebin(2); 
  xtal_binMaxVal = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()); 

  //Gen Stuff
  for(auto &p : *genParticle){
      genParticle_pdgId = p.pdgId();
      genParticle_status = p.status();
      genParticle_energy = p.energy();
      genParticle_pt = p.pt();
      genParticle_eta = p.eta();
      genParticle_phi = p.phi();
  }

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
  //Profile = new TH1F("Profile","Profile",260,0., 26.);
  //Profile = fs->make<TH1F>("Profile","Profile",260,0., 26.);  
  //Profile->Sumw2(); 
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
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

void simAnalysis::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/L");
   tree->Branch("runId", &runId, "runId/L");
   tree->Branch("genParticle_pdgId", &genParticle_pdgId, "genParticle_pdgId/I");
   tree->Branch("genParticle_status", &genParticle_status, "genParticle_status/D"); 
   tree->Branch("genParticle_energy", &genParticle_energy, "genParticle_energy/D");
   tree->Branch("genParticle_pt", &genParticle_pt, "genParticle_pt/D");
   tree->Branch("genParticle_eta", &genParticle_eta, "genParticle_eta/D");
   tree->Branch("genParticle_phi", &genParticle_phi, "genParticle_phi/D");
   tree->Branch("xtal_energyMax", &xtal_energyMax, "xtal_energyMax/D");
   tree->Branch("xtal_binMaxVal", &xtal_binMaxVal, "xtal_binMaxVal/D");
   tree->Branch("xtal_tMax", &xtal_tMax, "xtal_tMax/D");
   tree->Branch("xtal_tMaxError", &xtal_tMaxError, "xtal_tMaxError/D");
   tree->Branch("xtal_tMax_fitPar0", &xtal_tMax_fitPar0, "xtal_tMax_fitPar0/D"); 
   tree->Branch("xtal_tMax_fitPar1", &xtal_tMax_fitPar1, "xtal_tMax_fitPar1/D"); 
   tree->Branch("xtal_tMax_fitPar2", &xtal_tMax_fitPar2, "xtal_tMax_fitPar2/D"); 
   tree->Branch("xtal_tMax_fitParError0", &xtal_tMax_fitParError0, "xtal_tMax_fitParError0/D"); 
   tree->Branch("xtal_tMax_fitParError1", &xtal_tMax_fitParError1, "xtal_tMax_fitParError1/D"); 
   tree->Branch("xtal_tMax_fitParError2", &xtal_tMax_fitParError2, "xtal_tMax_fitParError2/D");  
   tree->Branch("xtal_tMax_fitStatus", &xtal_tMax_fitStatus, "xtal_tMax_fitStatus/I");  
   tree->Branch("xtal_ieta", &xtal_ieta, "xtal_ieta/I");
   tree->Branch("xtal_iphi", &xtal_iphi, "xtal_iphi/I");
   tree->Branch("xtal_iz", &xtal_iz, "xtal_iz/I");
}  

//define this as a plug-in
DEFINE_FWK_MODULE(simAnalysis);
