// -*- C++ -*-
//
// Package:    JetIdStudy
// Class:      JetIdStudy
// 
/**\class JetIdStudy JetIdStudy.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/JetIdStudy.cc

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

//new inclusion
#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"

//
// class declaration
//

class JetIdStudy : public edm::EDAnalyzer {
public:
  explicit JetIdStudy(const edm::ParameterSet&);
  ~JetIdStudy();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //FUNCTION
  void SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, std::vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  bool  ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho);

  EGammaMvaEleEstimator* myMVANonTrig;
 
  //HISTOGRAMS
  TH1D* eleTight;

  TH1D* jetPt;
  TH1D* jetMass;
  TH1D* jetEta;
  TH1D* jetSubjettiness;
  TH1D* jetMuonEF;
  TH1D* jetPhotonEF;
  TH1D* jetChargedEmEF;
  TH1D* jetNeutralHEF;
  TH1D* jetChargedHEF;
  TH1D* cutPt;
  TH1D* cutMass;
  TH1D* cutEta;
  TH1D* cutSubjettiness;
  TH1D* cutMuonEF;
  TH1D* cutPhotonEF;
  TH1D* cutChargedEmEF;
  TH1D* cutNeutralHEF;
  TH1D* cutChargedHEF;
  TH1D* cutMatching;

  // ----------member data ---------------------------
};

using namespace edm;
using namespace std;
using NSVfitStandalone::Vector;
using NSVfitStandalone::LorentzVector;
using NSVfitStandalone::MeasuredTauLepton;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetIdStudy::JetIdStudy(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
}


JetIdStudy::~JetIdStudy()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetIdStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  edm::Handle<pat::JetCollection> CA8JetswithQjets;
  iEvent.getByLabel("selectedPatJetsCA8CHSwithQjets", CA8JetswithQjets);
  edm::Handle<pat::JetCollection> CA8JetsPruned;
  iEvent.getByLabel("selectedPatJetsCA8CHSpruned", CA8JetsPruned);

  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  vector<math::PtEtaPhiELorentzVector> genquark;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())>0 && abs(genPart.pdgId())<7){
      const reco::Candidate * mom = genPart.mother();
      if(mom->pdgId()==23){
	math::PtEtaPhiELorentzVector gen_prov; gen_prov=genPart.p4();
	genquark.push_back(gen_prov);
      }
    }
  }

  bool matching=false;
  float dRmatch = 9999.; pat::JetCollection::const_iterator SelectedJet;
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    for(unsigned int i=0; i<genquark.size(); i++){
      if(ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i])<dRmatch && ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i])<0.8){
	dRmatch = ROOT::Math::VectorUtil::DeltaR(jet->p4(),genquark[i]);
	SelectedJet = jet;
	matching=true;
      }
    }
  }

  if(matching){
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(SelectedJet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.7 ){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass=jetPruned->mass();
      }
    }
    
    jetPt->Fill(SelectedJet->pt());
    jetMass->Fill(mass);
    jetEta->Fill(SelectedJet->eta());
    jetSubjettiness->Fill(SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1"));
    jetMuonEF->Fill(SelectedJet->muonEnergyFraction());
    jetPhotonEF->Fill(SelectedJet->photonEnergyFraction());
    jetChargedEmEF->Fill(SelectedJet->chargedEmEnergyFraction());
    jetNeutralHEF->Fill(SelectedJet->neutralHadronEnergyFraction());
    jetChargedHEF->Fill(SelectedJet->chargedHadronEnergyFraction());
    cutPt->Fill(SelectedJet->pt()>350);
    cutMass->Fill((mass>70 && mass<110));
    cutEta->Fill(abs(SelectedJet->eta())<2.4);
    cutSubjettiness->Fill((SelectedJet->userFloat("tau2")/SelectedJet->userFloat("tau1"))<0.75);
    cutMuonEF->Fill(SelectedJet->muonEnergyFraction()<0.99);
    cutPhotonEF->Fill(SelectedJet->photonEnergyFraction()<0.99);
    cutChargedEmEF->Fill(SelectedJet->chargedEmEnergyFraction()<0.99);
    cutNeutralHEF->Fill(SelectedJet->neutralHadronEnergyFraction()<0.99);
    cutChargedHEF->Fill(SelectedJet->chargedHadronEnergyFraction()>0.0);
  }
  cutMatching->Fill(matching);

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
JetIdStudy::beginJob()
{
  Service<TFileService> fs;
  jetPt           = fs->make<TH1D>("jetPt",           "jetPt",           3000, 0, 3000);
  jetMass         = fs->make<TH1D>("jetMass",         "jetMass",         1000, 0, 1000);
  jetEta          = fs->make<TH1D>("jetEta",          "jetEta",          800, -4, 4   );
  jetSubjettiness = fs->make<TH1D>("jetSubjettiness", "jetSubjettiness", 1000, 0, 1   );
  jetMuonEF       = fs->make<TH1D>("jetMuonEF",       "jetMuonEF",       100,  0, 1   );
  jetPhotonEF     = fs->make<TH1D>("jetPhotonEF",     "jetPhotonEF",     100,  0, 1   );
  jetChargedEmEF  = fs->make<TH1D>("jetChargedEmEF",  "jetChargedEmEF",  100,  0, 1   );
  jetNeutralHEF   = fs->make<TH1D>("jetNeutralHEF",   "jetNeutralHEF",   100,  0, 1   );
  jetChargedHEF   = fs->make<TH1D>("jetChargedHEF",   "jetChargedHEF",   100,  0, 1   );
  cutPt           = fs->make<TH1D>("cutPt",           "cutPt",           2, -0.5, 1.5 );
  cutMass         = fs->make<TH1D>("cutMass",         "cutMass",         2, -0.5, 1.5 );
  cutEta          = fs->make<TH1D>("cutEta",          "cutEta",          2, -0.5, 1.5 );
  cutSubjettiness = fs->make<TH1D>("cutSubjettiness", "cutSubjettiness", 2, -0.5, 1.5 );
  cutMuonEF       = fs->make<TH1D>("cutMuonEF",       "cutMuonEF",       2, -0.5, 1.5 );
  cutPhotonEF     = fs->make<TH1D>("cutPhotonEF",     "cutPhotonEF",     2, -0.5, 1.5 );
  cutChargedEmEF  = fs->make<TH1D>("cutChargedEmEF",  "cutChargedEmEF",  2, -0.5, 1.5 );
  cutNeutralHEF   = fs->make<TH1D>("cutNeutralHEF",   "cutNeutralHEF",   2, -0.5, 1.5 );
  cutChargedHEF   = fs->make<TH1D>("cutChargedHEF",   "cutChargedHEF",   2, -0.5, 1.5 );
  cutMatching     = fs->make<TH1D>("cutMatching",     "cutMatching",     2, -0.5, 1.5 );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetIdStudy::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
JetIdStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetIdStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetIdStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetIdStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetIdStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(JetIdStudy);
