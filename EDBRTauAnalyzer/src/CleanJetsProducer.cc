// -*- C++ -*-
//
// Package:    CleanJetsProducer
// Class:      CleanJetsProducer
// 
/**\class CleanJetsProducer CleanJetsProducer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/CleanJetsProducer.cc

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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 

//new inclusion
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDProducer.h"


using namespace edm;
using namespace std;
using namespace reco;


//
// class declaration
//

class CleanJetsProducer : public edm::EDProducer {
public:
  explicit CleanJetsProducer(const edm::ParameterSet&);
  ~CleanJetsProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

 

  //HISTOGRAMS
  //TH1D* metPt;

  // ----------member data ---------------------------
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
CleanJetsProducer::CleanJetsProducer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //produces<PFJetCollection>( "ak5PFJetsNoMu" );
  produces<PFJetCollection>( );
}


CleanJetsProducer::~CleanJetsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CleanJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<reco::PFJetCollection> PFJets;
  iEvent.getByLabel("ak5PFJets", PFJets);
  auto_ptr<PFJetCollection> SetOfJets( new PFJetCollection );

  Handle<MuonRefVector> muons;
  iEvent.getByLabel("muons", muons);

  Handle<reco::VertexCollection> vertices;
  //iEvent.getByLabel("offlinePrimaryVertices", vertices);
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);
  

  for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet){
    
    math::XYZTLorentzVector pfmomentum;  
    std::vector<edm::Ptr<Candidate> > jetConstituents;
    jetConstituents.clear();
    std::vector<reco::PFCandidatePtr> JetPFCands = iJet->getPFConstituents();
    reco::PFJet::Specific specs = iJet->getSpecific();
    
    for(std::vector<edm::Ptr<reco::PFCandidate> >::iterator i = JetPFCands.begin(); i != JetPFCands.end(); ++i){
      reco::PFCandidate pfCand = *i;
      
      if (pfCand.particleId() == 3){ //the PFCandidate is a muon
	reco::MuonRef theRecoMuon = pfCand.muonRef();
	if(muon::isSoftMuon(*theRecoMuon, primaryVertex) && theRecoMuon->pt()>5){
	  //cout<<iJet->pt()<<" "<<theRecoMuon->pt()<<endl;
	  specs.mMuonEnergy -= pfCand.p4().e();
	  specs.mMuonMultiplicity -= 1;
	  specs.mChargedMuEnergy -= pfCand.p4().e();
	  specs.mChargedMultiplicity -= 1;
	}
	else{
	  pfmomentum += pfCand.p4();
	  jetConstituents.push_back((*i));
	}
	
      }//the PFCandidate is a muon  
      else {
	pfmomentum += pfCand.p4();
	jetConstituents.push_back((*i));
      }
    }//loop over PFCandidate

    PFJet muonfreePFJet(pfmomentum, specs, jetConstituents);
    SetOfJets->push_back( muonfreePFJet );
  
  }//loop over PFJet

  //const OrphanHandle<PFJetCollection> cleanedJetsRefProd = iEvent.put( SetOfJets, "ak5PFJetsNoMu" ); 
  iEvent.put( SetOfJets );


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
CleanJetsProducer::beginJob()
{
  Service<TFileService> fs;
  //metPt             = fs->make<TH1D>("metPt",            "metPt",            5000, 0, 5000);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CleanJetsProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
CleanJetsProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CleanJetsProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CleanJetsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CleanJetsProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CleanJetsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(CleanJetsProducer);
