// -*- C++ -*-
//
// Package:    Test
// Class:      Test
// 
/**\class Test Test.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/Test.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Aniello Spiezia,21 1-007,+41227676459,
//         Created:  Mon Sep  9 13:14:05 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//new inclusion
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TMath.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"


//
// class declaration
//

class Test : public edm::EDAnalyzer {
public:
  explicit Test(const edm::ParameterSet&);
  ~Test();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  TH1D* jetPt_mh; TH1D* jetEta_mh; TH1D* jetMass_mh; 
  TH1D* jetMuonEnergyFraction_mh;
  TH1D* jetPhotonEnergyFraction_mh;
  TH1D* jetChargedEmEnergyFraction_mh;
  TH1D* jetNeutralHadronEnergyFraction_mh;
  TH1D* jetChargedHadronEnergyFraction_mh;

  TH1D* jetPt_hh; TH1D* jetEta_hh; TH1D* jetMass_hh; 
  TH1D* jetMuonEnergyFraction_hh;
  TH1D* jetPhotonEnergyFraction_hh;
  TH1D* jetChargedEmEnergyFraction_hh;
  TH1D* jetNeutralHadronEnergyFraction_hh;
  TH1D* jetChargedHadronEnergyFraction_hh;

  TH1D* jetPt_mh_NoMu; TH1D* jetEta_mh_NoMu; TH1D* jetMass_mh_NoMu; 
  TH1D* jetMuonEnergyFraction_mh_NoMu;
  TH1D* jetPhotonEnergyFraction_mh_NoMu;
  TH1D* jetChargedEmEnergyFraction_mh_NoMu;
  TH1D* jetNeutralHadronEnergyFraction_mh_NoMu;
  TH1D* jetChargedHadronEnergyFraction_mh_NoMu;

  TH1D* jetPt_hh_NoMu; TH1D* jetEta_hh_NoMu; TH1D* jetMass_hh_NoMu; 
  TH1D* jetMuonEnergyFraction_hh_NoMu;
  TH1D* jetPhotonEnergyFraction_hh_NoMu;
  TH1D* jetChargedEmEnergyFraction_hh_NoMu;
  TH1D* jetNeutralHadronEnergyFraction_hh_NoMu;
  TH1D* jetChargedHadronEnergyFraction_hh_NoMu;

  TH1D* jetCA8Mass; TH1D* jetCA8Pt; TH1D* jetCA8Tau21;

  // ----------member data ---------------------------
};

using namespace edm;
using namespace std;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Test::Test(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed

}


Test::~Test()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Test::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  
  //Handle<reco::PFJetCollection> PFJets;
  //iEvent.getByLabel("ak5PFJetsNoMu", PFJets);
  //Handle<reco::PFJetCollection> PFJetsSiMu;
  //iEvent.getByLabel("ak5PFJets", PFJetsSiMu);

  //Handle<vector<reco::GenParticle> > genParts;
  //iEvent.getByLabel("genParticles", genParts);

    //Handle<reco::PFTauCollection> tauHandle;
  //iEvent.getByLabel("hpsPFTauProducer",tauHandle);
  //Handle<reco::PFTauDiscriminator> decayModeFinding;
  //iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",decayModeFinding);

  edm::Handle<reco::PFJetCollection> CA8JetswithQjets;
  iEvent.getByLabel("ca8PFJetsCHSwithNsub", CA8JetswithQjets);
  edm::Handle<reco::BasicJetCollection> CA8JetsPruned;
  iEvent.getByLabel("ca8PFJetsCHSpruned", CA8JetsPruned);
  vector<reco::PFJetCollection::const_iterator> SelectedJet;
  vector<float> SelectedPrunedMass;
  int Njet; float massZ=-9999; float ptZ=-999; bool foundZ=false; float tau21Z=-9999;
  SelectedJet.clear(); SelectedPrunedMass.clear();
  for(reco::PFJetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(reco::BasicJetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.7 ){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass=jetPruned->mass();
      }
    }
    if(jet->muonEnergyFraction()>=0.99) continue;
    if(jet->photonEnergyFraction()>=0.99) continue;
    if(jet->chargedEmEnergyFraction()>=0.99) continue;
    if(jet->neutralHadronEnergyFraction()>=0.99) continue;
    if(jet->chargedHadronEnergyFraction()<=0.00) continue;
    if(jet->pt()<80) continue;
    if(abs(jet->eta())>2.4) continue;
    if(!(mass>70 && mass<110)) continue;
    if(jet->mass()>0.75) continue;
    //if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    //if(fabs(mass-91)<fabs(massZ-91)){
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->mass();
      foundZ=true;
    }
    Njet = Njet + 1;
    SelectedJet.push_back(jet);
    SelectedPrunedMass.push_back(mass);
   }
  if(foundZ){
    jetCA8Mass->Fill(massZ);
    jetCA8Pt->Fill(ptZ);
    jetCA8Tau21->Fill(tau21Z);
  }
		
  /*
  int ele = 0; int muo = 0;
  vector<const reco::Candidate *> GEN;
  math::PtEtaPhiELorentzVector gentau_prov;
  math::PtEtaPhiELorentzVector gentauHad_prov;
  math::PtEtaPhiELorentzVector gentauWithNeutrino_prov;
  math::PtEtaPhiELorentzVector gentauHadWithNeutrino_prov;
  vector<math::PtEtaPhiELorentzVector> gentau;
  vector<math::PtEtaPhiELorentzVector> gentauHad;
  vector<math::PtEtaPhiELorentzVector> gentauWithNeutrino;
  vector<math::PtEtaPhiELorentzVector> gentauHadWithNeutrino;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    const reco::Candidate * mom = genPart.mother();
    if(abs(genPart.pdgId())==15 && genPart.status()!=3 && (abs(mom->pdgId())==23 || abs(mom->pdgId())==15)){
      if((!genPart.pt()>10)) continue;
      gentau_prov=genPart.p4();
      gentauWithNeutrino_prov=genPart.p4();
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=13){
	  gentauHad_prov=genPart.p4();
	  gentauHadWithNeutrino_prov=genPart.p4();
	}
	if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	if(abs(daughter->pdgId())==13 && daughter->status()==1) {
	  gentauHad_prov = gentauHad_prov - daughter->p4();
	  muo = muo + 1;
	}
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  gentau_prov    = gentau_prov    - daughter->p4();
	  gentauHad_prov = gentauHad_prov - daughter->p4();
	}
      }
      if(gentau_prov.pt()>10){
	gentau.push_back(gentau_prov);
	gentauHad.push_back(gentauHad_prov);
	gentauWithNeutrino.push_back(gentauWithNeutrino_prov);
	gentauHadWithNeutrino.push_back(gentauWithNeutrino_prov);
      }
    }
  }

  if(ele==0 && muo==2){
    for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet){
      if(iJet->pt()<20) continue;
      jetPt_mh_NoMu->Fill(iJet->pt());
      jetEta_mh_NoMu->Fill(iJet->eta());
      jetMass_mh_NoMu->Fill(iJet->mass());
      jetMuonEnergyFraction_mh_NoMu->Fill(iJet->muonEnergyFraction());
      jetPhotonEnergyFraction_mh_NoMu->Fill(iJet->photonEnergyFraction());
      jetChargedEmEnergyFraction_mh_NoMu->Fill(iJet->chargedEmEnergyFraction());
      jetNeutralHadronEnergyFraction_mh_NoMu->Fill(iJet->neutralHadronEnergyFraction());
      jetChargedHadronEnergyFraction_mh_NoMu->Fill(iJet->chargedHadronEnergyFraction());
    }
    for(reco::PFJetCollection::const_iterator iJet = PFJetsSiMu->begin(); iJet != PFJetsSiMu->end(); ++iJet){
      if(iJet->pt()<20) continue;
      jetPt_mh->Fill(iJet->pt());
      jetEta_mh->Fill(iJet->eta());
      jetMass_mh->Fill(iJet->mass());
      jetMuonEnergyFraction_mh->Fill(iJet->muonEnergyFraction());
      jetPhotonEnergyFraction_mh->Fill(iJet->photonEnergyFraction());
      jetChargedEmEnergyFraction_mh->Fill(iJet->chargedEmEnergyFraction());
      jetNeutralHadronEnergyFraction_mh->Fill(iJet->neutralHadronEnergyFraction());
      jetChargedHadronEnergyFraction_mh->Fill(iJet->chargedHadronEnergyFraction());
    }
    
    float deltaR=99.;
    int tauIndex = -1;
    for(unsigned int i=0;i<tauHandle->size();++i) {
      reco::PFTauRef pfTauRef(tauHandle,i);
      if(!(pfTauRef->pt()>10)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<10 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<deltaR){
	deltaR = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4());
	tauIndex = i;
      }
    }
    if(tauIndex!=-1) {
      reco::PFTauRef pfTauRefSelected(tauHandle,tauIndex);
    }
  }

  if(ele==0 && muo==0){
    for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet){
      if(iJet->pt()<20) continue;
      jetPt_hh_NoMu->Fill(iJet->pt());
      jetEta_hh_NoMu->Fill(iJet->eta());
      jetMass_hh_NoMu->Fill(iJet->mass());
      jetMuonEnergyFraction_hh_NoMu->Fill(iJet->muonEnergyFraction());
      jetPhotonEnergyFraction_hh_NoMu->Fill(iJet->photonEnergyFraction());
      jetChargedEmEnergyFraction_hh_NoMu->Fill(iJet->chargedEmEnergyFraction());
      jetNeutralHadronEnergyFraction_hh_NoMu->Fill(iJet->neutralHadronEnergyFraction());
      jetChargedHadronEnergyFraction_hh_NoMu->Fill(iJet->chargedHadronEnergyFraction());
    }
    for(reco::PFJetCollection::const_iterator iJet = PFJetsSiMu->begin(); iJet != PFJetsSiMu->end(); ++iJet){
      if(iJet->pt()<20) continue;
      jetPt_hh->Fill(iJet->pt());
      jetEta_hh->Fill(iJet->eta());
      jetMass_hh->Fill(iJet->mass());
      jetMuonEnergyFraction_hh->Fill(iJet->muonEnergyFraction());
      jetPhotonEnergyFraction_hh->Fill(iJet->photonEnergyFraction());
      jetChargedEmEnergyFraction_hh->Fill(iJet->chargedEmEnergyFraction());
      jetNeutralHadronEnergyFraction_hh->Fill(iJet->neutralHadronEnergyFraction());
      jetChargedHadronEnergyFraction_hh->Fill(iJet->chargedHadronEnergyFraction());
    }
    
    float deltaR1=99.; float deltaR2=99.;
    int tauIndex1 = -1; int tauIndex2 = -1;
    for(unsigned int i=0;i<tauHandle->size();++i) {
      reco::PFTauRef pfTauRef(tauHandle,i);
      if(!(pfTauRef->pt()>10)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<1.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<deltaR1){
	deltaR1 = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4());
	tauIndex1 = i;
      }
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[1],pfTauRef->p4())<1.5 && ROOT::Math::VectorUtil::DeltaR(gentauHad[1],pfTauRef->p4())<deltaR2){
	deltaR2 = ROOT::Math::VectorUtil::DeltaR(gentauHad[1],pfTauRef->p4());
	tauIndex2 = i;
      }
    }
  }

  */

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Test::beginJob()
{
  Service<TFileService> fs;
  /*
  jetPt_mh   = fs->make<TH1D>("jetPt_mh",   "jetPt_mh",   1000, 0, 1000);
  jetEta_mh  = fs->make<TH1D>("jetEta_mh",  "jetEta_mh",  600, -3, 3   );
  jetMass_mh = fs->make<TH1D>("jetMass_mh", "jetMass_mh", 300,  0, 300 );
  jetMuonEnergyFraction_mh          = fs->make<TH1D>("jetMuonEnergyFraction_mh",          "jetMuonEnergyFraction_mh",          100, 0, 1);
  jetPhotonEnergyFraction_mh        = fs->make<TH1D>("jetPhotonEnergyFraction_mh",        "jetPhotonEnergyFraction_mh",        100, 0, 1);
  jetChargedEmEnergyFraction_mh     = fs->make<TH1D>("jetChargedEmEnergyFraction_mh",     "jetChargedEmEnergyFraction_mh",     100, 0, 1);
  jetNeutralHadronEnergyFraction_mh = fs->make<TH1D>("jetNeutralHadronEnergyFraction_mh", "jetNeutralHadronEnergyFraction_mh", 100, 0, 1);
  jetChargedHadronEnergyFraction_mh = fs->make<TH1D>("jetChargedHadronEnergyFraction_mh", "jetChargedHadronEnergyFraction_mh", 100, 0, 1);
  jetPt_hh   = fs->make<TH1D>("jetPt_hh",   "jetPt_hh",   1000, 0, 1000);
  jetEta_hh  = fs->make<TH1D>("jetEta_hh",  "jetEta_hh",  600, -3, 3   );
  jetMass_hh = fs->make<TH1D>("jetMass_hh", "jetMass_hh", 300,  0, 300 );
  jetMuonEnergyFraction_hh          = fs->make<TH1D>("jetMuonEnergyFraction_hh",          "jetMuonEnergyFraction_hh",          100, 0, 1);
  jetPhotonEnergyFraction_hh        = fs->make<TH1D>("jetPhotonEnergyFraction_hh",        "jetPhotonEnergyFraction_hh",        100, 0, 1);
  jetChargedEmEnergyFraction_hh     = fs->make<TH1D>("jetChargedEmEnergyFraction_hh",     "jetChargedEmEnergyFraction_hh",     100, 0, 1);
  jetNeutralHadronEnergyFraction_hh = fs->make<TH1D>("jetNeutralHadronEnergyFraction_hh", "jetNeutralHadronEnergyFraction_hh", 100, 0, 1);
  jetChargedHadronEnergyFraction_hh = fs->make<TH1D>("jetChargedHadronEnergyFraction_hh", "jetChargedHadronEnergyFraction_hh", 100, 0, 1);

  jetPt_mh_NoMu   = fs->make<TH1D>("jetPt_mh_NoMu",   "jetPt_mh_NoMu",   1000, 0, 1000);
  jetEta_mh_NoMu  = fs->make<TH1D>("jetEta_mh_NoMu",  "jetEta_mh_NoMu",  600, -3, 3   );
  jetMass_mh_NoMu = fs->make<TH1D>("jetMass_mh_NoMu", "jetMass_mh_NoMu", 300,  0, 300 );
  jetMuonEnergyFraction_mh_NoMu          = fs->make<TH1D>("jetMuonEnergyFraction_mh_NoMu",          "jetMuonEnergyFraction_mh_NoMu",          100, 0, 1);
  jetPhotonEnergyFraction_mh_NoMu        = fs->make<TH1D>("jetPhotonEnergyFraction_mh_NoMu",        "jetPhotonEnergyFraction_mh_NoMu",        100, 0, 1);
  jetChargedEmEnergyFraction_mh_NoMu     = fs->make<TH1D>("jetChargedEmEnergyFraction_mh_NoMu",     "jetChargedEmEnergyFraction_mh_NoMu",     100, 0, 1);
  jetNeutralHadronEnergyFraction_mh_NoMu = fs->make<TH1D>("jetNeutralHadronEnergyFraction_mh_NoMu", "jetNeutralHadronEnergyFraction_mh_NoMu", 100, 0, 1);
  jetChargedHadronEnergyFraction_mh_NoMu = fs->make<TH1D>("jetChargedHadronEnergyFraction_mh_NoMu", "jetChargedHadronEnergyFraction_mh_NoMu", 100, 0, 1);
  jetPt_hh_NoMu   = fs->make<TH1D>("jetPt_hh_NoMu",   "jetPt_hh_NoMu",   1000, 0, 1000);
  jetEta_hh_NoMu  = fs->make<TH1D>("jetEta_hh_NoMu",  "jetEta_hh_NoMu",  600, -3, 3   );
  jetMass_hh_NoMu = fs->make<TH1D>("jetMass_hh_NoMu", "jetMass_hh_NoMu", 300,  0, 300 );
  jetMuonEnergyFraction_hh_NoMu          = fs->make<TH1D>("jetMuonEnergyFraction_hh_NoMu",          "jetMuonEnergyFraction_hh_NoMu",          100, 0, 1);
  jetPhotonEnergyFraction_hh_NoMu        = fs->make<TH1D>("jetPhotonEnergyFraction_hh_NoMu",        "jetPhotonEnergyFraction_hh_NoMu",        100, 0, 1);
  jetChargedEmEnergyFraction_hh_NoMu     = fs->make<TH1D>("jetChargedEmEnergyFraction_hh_NoMu",     "jetChargedEmEnergyFraction_hh_NoMu",     100, 0, 1);
  jetNeutralHadronEnergyFraction_hh_NoMu = fs->make<TH1D>("jetNeutralHadronEnergyFraction_hh_NoMu", "jetNeutralHadronEnergyFraction_hh_NoMu", 100, 0, 1);
  jetChargedHadronEnergyFraction_hh_NoMu = fs->make<TH1D>("jetChargedHadronEnergyFraction_hh_NoMu", "jetChargedHadronEnergyFraction_hh_NoMu", 100, 0, 1);
*/

  jetCA8Mass = fs->make<TH1D>("jetCA8Mass", "jetCA8Mass", 1000, 0, 200);
  jetCA8Pt = fs->make<TH1D>("jetCA8Pt", "jetCA8Pt", 2000, 0, 2000);
  jetCA8Tau21 = fs->make<TH1D>("jetCA8Tau21", "jetCA8Tau21", 300, -1, 2);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Test::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Test::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Test::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Test::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Test::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Test::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(Test);
