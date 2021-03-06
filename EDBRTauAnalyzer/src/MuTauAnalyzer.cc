// -*- C++ -*-
//
// Package:    MuTauAnalyzer
// Class:      MuTauAnalyzer
// 
/**\class MuTauAnalyzer MuTauAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/MuTauAnalyzer.cc

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
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "TMath.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"


//
// class declaration
//

class MuTauAnalyzer : public edm::EDAnalyzer {
public:
  explicit MuTauAnalyzer(const edm::ParameterSet&);
  ~MuTauAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  float MuonPFIso(reco::MuonCollection::const_iterator muon, bool highpt);
  bool signal_;

  TH1D* efficiency;
  TH1D* jetMass; TH1D* jetPt; TH1D* jetTau21;
  TH1D* tauPt; TH1D* muonPt; TH1D* metPt;
  TH1D* muonDRMet; TH1D* muonDPhiMet; TH1D* muonPFIso; 
  TH1D* jetDRMet; TH1D* jetDPhiMet; TH1D* tauDRMet; TH1D* tauDPhiMet;
  TH1D* jetDRTau; TH1D* jetDRMuo; TH1D* tauDRMuo;
  TH1D* ZDRZ; TH1D* MassSvfitTauMuo; TH1D* XMassSVFit;

  TH1D* jetPtNO;
  TH1D* jetMassNO;
  TH1D* jetTau21NO;
  TH1D* tauPtNO;
  TH1D* tauDRMuoNO;
  TH1D* muonPtNO;
  TH1D* metPtNO;
  TH1D* jetDPhiMetNO;
  TH1D* muonDPhiMetNO;
  TH1D* tauDPhiMetNO;
  TH1D* jetDRTauNO;
  TH1D* jetDRMuoNO;

  TTree *TreeVariable;
  float m_jetMass;
  float m_jetPt;
  float m_jetTau21;
  float m_tauPt;
  float m_muonPt;
  float m_metPt;
  float m_jetDPhiMet;
  float m_jetDRMet;
  float m_jetDRMuo;
  float m_jetDRTau;
  float m_tauDPhiMet;
  float m_tauDRMet;
  float m_tauDRMuo;
  float m_muonDPhiMet;
  float m_muonDRMet;
  float m_muonPFIso;
  float m_ZDRZ;
  float m_MassSvfitTauMuo;
  float m_XMassSVFit;

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
MuTauAnalyzer::MuTauAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  signal_ = iConfig.getUntrackedParameter<bool>("signal");
}


MuTauAnalyzer::~MuTauAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //GEN INFORMATION
  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);
  int ele = 0; int muo = 0; int neu = 0;
  vector<const reco::Candidate *> GEN;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==15 && genPart.status()!=3){
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	if(abs(daughter->pdgId())==13 && daughter->status()==1) muo = muo + 1;
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 || abs(daughter->pdgId())==16) && daughter->status()==1) neu = neu + 1;
      }
    }
  }
  
  //JET SELECTION
  edm::Handle<reco::PFJetCollection> CA8JetswithQjets;
  iEvent.getByLabel("ca8PFJetsCHSwithNsub", CA8JetswithQjets);
  edm::Handle<reco::BasicJetCollection> CA8JetsPruned;
  iEvent.getByLabel("ca8PFJetsCHSpruned", CA8JetsPruned);
  reco::PFJetCollection::const_iterator SelectedJet;
  float massZ=-9999; float ptZ=-999; bool foundJet=false; float tau21Z=-9999;
  for(reco::PFJetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(reco::BasicJetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
      float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
      if(dRtmp<dRmin && dRtmp<0.8){//matching failed if greater than jet radius
	dRmin=dRtmp;
	mass=jetPruned->mass();
      }
    }
    if(jet->muonEnergyFraction()>=0.99) continue;
    if(jet->photonEnergyFraction()>=0.99) continue;
    if(jet->chargedEmEnergyFraction()>=0.99) continue;
    if(jet->neutralHadronEnergyFraction()>=0.99) continue;
    if(jet->chargedHadronEnergyFraction()<=0.00) continue;
    if(jet->pt()<400) continue;
    if(abs(jet->eta())>2.4) continue;
    if(!(mass>70 && mass<110)) continue;
    if(jet->mass()>0.75) continue;
    foundJet=true;
    if(jet->pt()>ptZ){
      massZ=mass;
      ptZ=jet->pt();
      tau21Z=jet->mass();
      SelectedJet=jet;
    }
  }

  //TAU SELECTION
  Handle<reco::PFTauCollection> tauHandle;
  iEvent.getByLabel("hpsPFTauProducer",tauHandle);
  Handle<reco::PFTauDiscriminator> decayModeFinding;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",decayModeFinding);
  Handle<reco::PFTauDiscriminator> ByVLooseCombinedIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr",ByVLooseCombinedIsolationDBSumPtCorr);
  Handle<reco::PFTauDiscriminator> againstElectronLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection",againstElectronLoose);
  Handle<reco::PFTauDiscriminator> againstMuonLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection",againstMuonLoose);
  Handle<reco::PFTauDiscriminator> againstMuonTight3;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightMuonRejection3",againstMuonTight3);
  float ptTau=-99; bool foundTau=false;
  reco::PFTauRef SelectedTau;
  for(unsigned int i=0;i<tauHandle->size();++i) {  
    reco::PFTauRef PFTau(tauHandle,i);
    if(PFTau->pt()<20) continue;
    if(abs(PFTau->eta())>2.4) continue;
    if((*decayModeFinding)[PFTau]<0.5) continue;
    if((*ByVLooseCombinedIsolationDBSumPtCorr)[PFTau]<0.5) continue;
    if((*againstElectronLoose)[PFTau]<0.5) continue;
    if((*againstMuonLoose)[PFTau]<0.5) continue;
    foundTau=true;
    if(PFTau->pt()>ptTau){
      SelectedTau=PFTau;
      ptTau=PFTau->pt();
    }
  }

  //MUON SELECTION
  edm::Handle<reco::MuonCollection> muoH;
  iEvent.getByLabel("muons", muoH);
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);
  float ptMuon=-99; bool foundMuon=false;
  reco::MuonCollection::const_iterator SelectedMuon;
  for(reco::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(abs(cktTrack->eta())>2.4) continue;
    if(abs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatchedStations()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    foundMuon=true;
    if(cktTrack->pt()>ptMuon){
      SelectedMuon=muon;
      ptMuon=cktTrack->pt();
    }
  }

  edm::Handle<reco::PFMETCollection> met;
  iEvent.getByLabel("pfMet", met);

  //PLOT
  efficiency->Fill(0);
  if(signal_){
    if(ele==0 && muo==1){
      efficiency->Fill(1);
      if(foundJet && foundTau && foundMuon){
	
	TMatrixD covMET(2, 2); // PFMET significance matrix
	covMET[0][0] = (met->front() ).getSignificanceMatrix()(0,0);
	covMET[1][0] = (met->front() ).getSignificanceMatrix()(1,0);
	covMET[0][1] = (met->front() ).getSignificanceMatrix()(0,1);
	covMET[1][1] = (met->front() ).getSignificanceMatrix()(1,1);
	std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedTau->p4()));
	measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedMuon->p4()));
	NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
	algo.addLogM(false);
	algo.integrateMarkovChain();
	if(algo.pt()>0){
	  
	  TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
	  math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
	  TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E()); 
	  m_jetMass=massZ;
	  m_jetPt=SelectedJet->pt();
	  m_jetTau21=tau21Z;
	  m_tauPt=SelectedTau->pt();
	  m_muonPt=SelectedMuon->pt();
	  m_metPt=met->begin()->pt();
	  m_jetDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
	  m_jetDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
	  m_jetDRMuo=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
	  m_jetDRTau=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4());
	  m_tauDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4());
	  m_tauDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTau->p4());
	  m_tauDRMuo=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4());
	  m_muonDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4());
	  m_muonDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4());
	  m_muonPFIso=MuonPFIso(SelectedMuon,true);
	  m_ZDRZ=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
	  m_MassSvfitTauMuo=algo.getMass();
	  m_XMassSVFit=(SVFitTauTau+PrunedJet).M();
	  TreeVariable->Fill();
	  
	  if(met->begin()->pt()>50
	     && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>2.5 
	     && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>2.5
	     && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
	     && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
	     && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
	     && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()))>2 ){
	    if(algo.getMass()>20){
	      
	      efficiency->Fill(2);
	      
	      jetPt->Fill(SelectedJet->pt());
	      jetMass->Fill(massZ);
	      jetTau21->Fill(tau21Z);
	      tauPt->Fill(SelectedTau->pt());
	      muonPt->Fill(SelectedMuon->pt());
	      metPt->Fill(met->begin()->pt());
	      jetDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()));
	      jetDRMuo->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4()));
	      jetDRTau->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4()));
	      tauDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()));
	      muonDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()));
	      tauDRMuo->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4()));
	      jetDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4()));
	      tauDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTau->p4()));
	      muonDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4()));
	      muonPFIso->Fill(MuonPFIso(SelectedMuon,true));
	      ZDRZ->Fill(ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4()));
	      MassSvfitTauMuo->Fill(algo.getMass());
	      XMassSVFit->Fill((SVFitTauTau+PrunedJet).M());
	    }
	  }
	}
      }
    }
  }
  else{
    if(foundJet && foundTau && foundMuon){
	
      TMatrixD covMET(2, 2); // PFMET significance matrix
      covMET[0][0] = (met->front() ).getSignificanceMatrix()(0,0);
      covMET[1][0] = (met->front() ).getSignificanceMatrix()(1,0);
      covMET[0][1] = (met->front() ).getSignificanceMatrix()(0,1);
      covMET[1][1] = (met->front() ).getSignificanceMatrix()(1,1);
      std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
      measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedTau->p4()));
      measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, SelectedMuon->p4()));
      NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
      algo.addLogM(false);
      algo.integrateMarkovChain();
      if(algo.pt()>0){
	
	TLorentzVector SVFitTauTau; SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
	math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet->pt(),SelectedJet->eta(),SelectedJet->phi(),massZ);
	TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E()); 
	m_jetMass=massZ;
	m_jetPt=SelectedJet->pt();
	m_jetTau21=tau21Z;
	m_tauPt=SelectedTau->pt();
	m_muonPt=SelectedMuon->pt();
	m_metPt=met->begin()->pt();
	m_jetDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4());
	m_jetDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4());
	m_jetDRMuo=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4());
	m_jetDRTau=ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4());
	m_tauDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4());
	m_tauDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTau->p4());
	m_tauDRMuo=ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4());
	m_muonDPhiMet=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4());
	m_muonDRMet=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4());
	m_muonPFIso=MuonPFIso(SelectedMuon,true);
	m_ZDRZ=ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4());
	m_MassSvfitTauMuo=algo.getMass();
	m_XMassSVFit=(SVFitTauTau+PrunedJet).M();
	TreeVariable->Fill();

	if(met->begin()->pt()>50
	   && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>2.5 
	   && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>2.5
	   && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
	   && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
	   && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
	   && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()))>2 ){
	  if(algo.getMass()>20){
	  
	    efficiency->Fill(1);
	    jetMass->Fill(massZ);
	    jetPt->Fill(SelectedJet->pt());
	    jetTau21->Fill(tau21Z);
	    tauPt->Fill(SelectedTau->pt());
	    muonPt->Fill(SelectedMuon->pt());
	    metPt->Fill(met->begin()->pt());
	    jetDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()));
	    jetDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4()));
	    jetDRMuo->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4()));
	    jetDRTau->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4()));
	    tauDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()));
	    tauDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTau->p4()));
	    tauDRMuo->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4()));
	    muonDPhiMet->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()));
	    muonDRMet->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4()));
	    muonPFIso->Fill(MuonPFIso(SelectedMuon,true));
	    ZDRZ->Fill(ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4()));
	    MassSvfitTauMuo->Fill(algo.getMass());
	    XMassSVFit->Fill((SVFitTauTau+PrunedJet).M());
	  }
	}
      }
    }
  }


  //CONTROL PLOT
  if(foundJet && foundTau && foundMuon){
    metPtNO->Fill(met->begin()->pt());
    if(met->begin()->pt()>50){
      jetDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()));
      muonDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()));
      tauDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()));
      jetDRTauNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4()));
      jetDRMuoNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4()));
      tauDRMuoNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4()));
    }
  }

  if(foundTau && foundMuon){
    reco::PFJetCollection::const_iterator SelectedJet1;
    reco::PFJetCollection::const_iterator SelectedJet2;
    reco::PFJetCollection::const_iterator SelectedJet3;
    float ptMax1=-999; float ptMax2=-999; float ptMax3=-999;
    float massSelected = 0.;
    for(reco::PFJetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
      float dRmin = 9999.; float mass = 0.;
      for(reco::BasicJetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
	float dRtmp = ROOT::Math::VectorUtil::DeltaR(jet->p4(),jetPruned->p4());
	if(dRtmp<dRmin && dRtmp<0.8){
	  dRmin=dRtmp;
	  mass=jetPruned->mass();
	}
      }
      if(jet->muonEnergyFraction()>=0.99) continue;
      if(jet->photonEnergyFraction()>=0.99) continue;
      if(jet->chargedEmEnergyFraction()>=0.99) continue;
      if(jet->neutralHadronEnergyFraction()>=0.99) continue;
      if(jet->chargedHadronEnergyFraction()<=0.00) continue;
      if(abs(jet->eta())>2.4) continue;
      if(jet->pt()>400 && (mass>20 && mass<70)){
	if(jet->pt()>ptMax1){
	  ptMax1=jet->pt();
	  SelectedJet1=jet;
	}
      }
      if(jet->mass()<0.75 && (mass>20 && mass<70)){
	if(jet->pt()>ptMax2){
	  ptMax2=jet->pt();
	  SelectedJet2=jet;
	}
      }
      if(jet->pt()>400 && jet->mass()<0.75){
	if(jet->pt()>ptMax3){
	  ptMax3=jet->pt();
	  massSelected=mass;
	  SelectedJet3=jet;
	}
      }
    }
    if(ptMax1>0 && met->begin()->pt()>50
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet1->p4())>2.5 
       && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet1->p4())>2.5
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet1->p4()))>2 ) jetTau21NO->Fill(SelectedJet1->mass());
    if(ptMax2>0 && met->begin()->pt()>50
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet2->p4())>2.5 
       && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet2->p4())>2.5
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet2->p4()))>2 ) jetPtNO->Fill(SelectedJet2->pt());
    if(ptMax3>0 && met->begin()->pt()>50
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet3->p4())>2.5 
       && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet3->p4())>2.5
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet3->p4()))>2 ) jetMassNO->Fill(massSelected);
  }

  if(foundJet && foundMuon){
    ptTau=-99;
    for(unsigned int i=0;i<tauHandle->size();++i) {  
      reco::PFTauRef PFTau(tauHandle,i);
      if(abs(PFTau->eta())>2.4) continue;
      if((*decayModeFinding)[PFTau]<0.5) continue;
      if((*ByVLooseCombinedIsolationDBSumPtCorr)[PFTau]<0.5) continue;
      if((*againstElectronLoose)[PFTau]<0.5) continue;
      if((*againstMuonLoose)[PFTau]<0.5) continue;
      if(PFTau->pt()>ptTau){
	SelectedTau=PFTau;
	ptTau=PFTau->pt();
      }
    }
    if(ptTau>0 && met->begin()->pt()>50
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>2.5 
       && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>2.5
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()))>2 ) tauPtNO->Fill(SelectedTau->pt());
  }

  if(foundJet && foundTau){
    ptMuon=-99;
    for(reco::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
      if(!(muon->isGlobalMuon())) continue;
      reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
      if(abs(cktTrack->eta())>2.4) continue;
      if(abs(cktTrack->phi())>3.2) continue;
      if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
      if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
      if(muon->numberOfMatchedStations()<=1) continue;
      if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
      if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
      if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
      if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
      foundMuon=true;
      if(cktTrack->pt()>ptMuon){
	SelectedMuon=muon;
	ptMuon=cktTrack->pt();
      }
    }
    if(ptMuon>0 && met->begin()->pt()>50
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>2.5 
       && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>2.5
       && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
       && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()))>2 ) muonPtNO->Fill(ptMuon);
  }



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
MuTauAnalyzer::beginJob()
{
  Service<TFileService> fs;
  TreeVariable = fs->make<TTree>("TreeVariable", "TreeVariable");
  TreeVariable->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeVariable->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeVariable->Branch("jetTau21", &m_jetTau21, "jetTau21/f");
  TreeVariable->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeVariable->Branch("muonPt", &m_muonPt, "muonPt/f");
  TreeVariable->Branch("metPt", &m_metPt, "metPt/f");
  TreeVariable->Branch("jetDPhiMet", &m_jetDPhiMet, "jetDPhiMet/f");
  TreeVariable->Branch("jetDRMet", &m_jetDRMet, "jetDRMet/f");
  TreeVariable->Branch("jetDRMuo", &m_jetDRMuo, "jetDRMuo/f");
  TreeVariable->Branch("jetDRTau", &m_jetDRTau, "jetDRTau/f");
  TreeVariable->Branch("tauDPhiMet", &m_tauDPhiMet, "tauDPhiMet/f");
  TreeVariable->Branch("tauDRMet", &m_tauDRMet, "tauDRMet/f");
  TreeVariable->Branch("tauDRMuo", &m_tauDRMuo, "tauDRMuo/f");
  TreeVariable->Branch("muonDPhiMet", &m_muonDPhiMet, "muonDPhiMet/f");
  TreeVariable->Branch("muonDRMet", &m_muonDRMet, "muonDRMet/f");
  TreeVariable->Branch("muonPFIso", &m_muonPFIso, "muonPFIso/f");
  TreeVariable->Branch("ZDRZ", &m_ZDRZ, "ZDRZ/f");
  TreeVariable->Branch("MassSvfitTauMuo", &m_MassSvfitTauMuo, "MassSvfitTauMuo/f");
  TreeVariable->Branch("XMassSVFit", &m_XMassSVFit, "XMassSVFit/f");

  efficiency      = fs->make<TH1D>("efficiency",      "efficiency",      3, -0.5, 2.5);
  jetMass         = fs->make<TH1D>("jetMass",         "jetMass",         400, 0, 200);
  jetPt           = fs->make<TH1D>("jetPt",           "jetPt",           2000, 0, 2000);
  jetTau21        = fs->make<TH1D>("jetTau21",        "jetTau21",        100, 0, 1);
  jetDRMet        = fs->make<TH1D>("jetDRMet",        "jetDRMet",        100, 0, 5);
  jetDRTau        = fs->make<TH1D>("jetDRTau",        "jetDRTau",        100, 0, 5);
  jetDRMuo        = fs->make<TH1D>("jetDRMuo",        "jetDRMuo",        100, 0, 5);
  jetDPhiMet      = fs->make<TH1D>("jetDPhiMet",      "jetDPhiMet",      100, -5, 5);
  tauPt           = fs->make<TH1D>("tauPt",           "tauPt",           2000, 0, 2000);
  tauDRMet        = fs->make<TH1D>("tauDRMet",        "tauDRMet",        100, 0, 5);
  tauDRMuo        = fs->make<TH1D>("tauDRMuo",        "tauDRMuo",        100, 0, 5);
  tauDPhiMet      = fs->make<TH1D>("tauDPhiMet",      "tauDPhiMet",      100, -5, 5);
  muonPt          = fs->make<TH1D>("muonPt",          "muonPt",          2000, 0, 2000);
  muonDRMet       = fs->make<TH1D>("muonDRMet",       "muonDRMet",       100, 0, 5);
  muonDPhiMet     = fs->make<TH1D>("muonDPhiMet",     "muonDPhiMet",     100, -5, 5);
  muonPFIso       = fs->make<TH1D>("muonPFIso",       "muonPFIso",       500, 0, 5);
  metPt           = fs->make<TH1D>("metPt",           "metPt",           2000, 0, 2000);
  ZDRZ            = fs->make<TH1D>("ZDRZ",            "ZDRZ",            100, 0, 5);
  MassSvfitTauMuo = fs->make<TH1D>("MassSvfitTauMuo", "MassSvfitTauMuo", 500, 0, 500);
  XMassSVFit      = fs->make<TH1D>("XMassSVFit",      "XMassSVFit",      300, 0, 3000);

  jetPtNO           = fs->make<TH1D>("jetPtNO",           "jetPtNO",           2000, 0, 2000);
  jetMassNO         = fs->make<TH1D>("jetMassNO",         "jetMassNO",         400, 0, 200);
  jetTau21NO        = fs->make<TH1D>("jetTau21NO",        "jetTau21NO",        100, 0, 1);
  tauPtNO           = fs->make<TH1D>("tauPtNO",           "tauPtNO",           2000, 0, 2000);
  muonPtNO          = fs->make<TH1D>("muonPtNO",          "muonPtNO",          2000, 0, 2000);
  metPtNO           = fs->make<TH1D>("metPtNO",           "metPtNO",           2000, 0, 2000);
  jetDPhiMetNO      = fs->make<TH1D>("jetDPhiMetNO",      "jetDPhiMetNO",      100, -5, 5);
  muonDPhiMetNO     = fs->make<TH1D>("muonDPhiMetNO",     "muonDPhiMetNO",     100, -5, 5);
  tauDPhiMetNO      = fs->make<TH1D>("tauDPhiMetNO",      "tauDPhiMetNO",      100, -5, 5);
  tauDRMuoNO        = fs->make<TH1D>("tauDRMuoNO",        "tauDRMuoNO",        100, 0, 5);
  jetDRTauNO        = fs->make<TH1D>("jetDRTauNO",        "jetDRTauNO",        100, 0, 5);
  jetDRMuoNO        = fs->make<TH1D>("jetDRMuoNO",        "jetDRMuoNO",        100, 0, 5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuTauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MuTauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MuTauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuTauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuTauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float MuTauAnalyzer::MuonPFIso(reco::MuonCollection::const_iterator muon, bool highpt){
  float sumChargedHadronPt = muon->pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon->pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon->pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon->pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt){
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuTauAnalyzer);
