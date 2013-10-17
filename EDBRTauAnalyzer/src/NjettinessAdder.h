#ifndef ExoDiBosonResonances_EDBRTauAnalyzer_NjettinessAdder_h
#define ExoDiBosonResonances_EDBRTauAnalyzer_NjettinessAdder_h

#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class NjettinessAdder : public edm::EDProducer { 
public:
  explicit NjettinessAdder(const edm::ParameterSet& iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    cone_(iConfig.getParameter<double>("cone"))
    {
    produces<std::vector<reco::PFJet> >();
  }
  
  virtual ~NjettinessAdder() {}
  
  void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) ;
  float getTau(int num,edm::Ptr<reco::PFJet> object) const;

private:	
  edm::InputTag src_ ;
  double cone_ ;
};


#endif
