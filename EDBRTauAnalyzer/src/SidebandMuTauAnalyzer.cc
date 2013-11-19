// -*- C++ -*-
//
// Package:    SidebandMuTauAnalyzer
// Class:      SidebandMuTauAnalyzer
// 
/**\class SidebandMuTauAnalyzer SidebandMuTauAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/SidebandMuTauAnalyzer.cc

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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TMath.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"


//
// class declaration
//

class SidebandMuTauAnalyzer : public edm::EDAnalyzer {
public:
  explicit SidebandMuTauAnalyzer(const edm::ParameterSet&);
  ~SidebandMuTauAnalyzer();
  
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

  TH1D* NVertices;
  TH1D* jetMass; TH1D* jetPt; TH1D* jetTau21;
  TH1D* tauPt; TH1D* muonPt; TH1D* metPt;
  TH1D* muonDRMet; TH1D* muonDPhiMet; TH1D* muonPFIso; 
  TH1D* jetDRMet; TH1D* jetDPhiMet; TH1D* tauDRMet; TH1D* tauDPhiMet;
  TH1D* jetDRTau; TH1D* jetDRMuo; TH1D* tauDRMuo;
  TH1D* ZDRZ; TH1D* MassSvfitTauMuo; TH1D* XMassSVFit;
  TH1D* metPx;  TH1D* metPy;  TH1D* metEta;  TH1D* metPhi;

  TH1D* NVerticesNO;
  TH1D* jetMassNO; TH1D* jetPtNO; TH1D* jetTau21NO;
  TH1D* tauPtNO; TH1D* muonPtNO; TH1D* metPtNO;
  TH1D* muonDRMetNO; TH1D* muonDPhiMetNO; TH1D* muonPFIsoNO; 
  TH1D* jetDRMetNO; TH1D* jetDPhiMetNO; TH1D* tauDRMetNO; TH1D* tauDPhiMetNO;
  TH1D* jetDRTauNO; TH1D* jetDRMuoNO; TH1D* tauDRMuoNO;
  TH1D* ZDRZNO; TH1D* MassSvfitTauMuoNO; TH1D* XMassSVFitNO;
  TH1D* metPxNO;  TH1D* metPyNO;  TH1D* metEtaNO;  TH1D* metPhiNO;

  TTree *TreeVariable;
  int   m_NVertices;
  float m_jetMass;
  float m_jetPt;
  float m_jetTau21;
  float m_tauPt;
  float m_muonPt;
  float m_metPt;
  float m_metPx;
  float m_metPy;
  float m_metEta;
  float m_metPhi;
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
  float m_PUWeight;
  
  edm::LumiReWeighting LumiWeights_;
  bool isData; 

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
SidebandMuTauAnalyzer::SidebandMuTauAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  isData = iConfig.getUntrackedParameter<bool>("isData_");


  // True number of interaction for data produced as in: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
  TFile *da_=new TFile ("/afs/cern.ch/work/a/aspiezia/EDBRTauAnalyzer/CMSSW_5_3_11_patch6/src/ExoDiBosonResonances/EDBRTauAnalyzer/data/MyDataPileupHistogram_True.root");
  TH1F *da = (TH1F*) da_->Get("pileup");
  
  // MC distribution of true number of interactions as in: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
  Double_t dat[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06 };
  
  //PileUp weights calculation
  double d,m;
  std::vector< float > mcNum; 
  std::vector< float > dataNum;
  for (Int_t i=1; i< 50; i++){
    m=dat[i-1];
    d=da->GetBinContent(i);
    mcNum.push_back(m);
    dataNum.push_back(d); 
  }
  LumiWeights_=edm::LumiReWeighting(mcNum, dataNum);
}


SidebandMuTauAnalyzer::~SidebandMuTauAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SidebandMuTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);
  NVertices->Fill(vertices->size());

  //Trigget paths  
  bool isFired_HLT_HT650 = false;
  bool isFired_HLT_PFJet320 = false;
  edm::Handle<edm::TriggerResults> trigResults;
  edm::InputTag trigResultsTag("TriggerResults","","HLT");  
  iEvent.getByLabel(trigResultsTag,trigResults);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
  unsigned int TrggIndex_PFHT650_v5( trigNames.triggerIndex("HLT_PFHT650_v5"));
  unsigned int TrggIndex_PFHT650_v6( trigNames.triggerIndex("HLT_PFHT650_v6"));
  unsigned int TrggIndex_PFHT650_v7( trigNames.triggerIndex("HLT_PFHT650_v7"));
  unsigned int TrggIndex_PFHT650_v8( trigNames.triggerIndex("HLT_PFHT650_v8"));
  unsigned int TrggIndex_PFHT650_v9( trigNames.triggerIndex("HLT_PFHT650_v9"));
  unsigned int TrggIndex_PFNoPUHT650_v1( trigNames.triggerIndex("HLT_PFNoPUHT650_v1"));
  unsigned int TrggIndex_PFNoPUHT650_v3( trigNames.triggerIndex("HLT_PFNoPUHT650_v3"));
  unsigned int TrggIndex_PFNoPUHT650_v4( trigNames.triggerIndex("HLT_PFNoPUHT650_v4"));
  unsigned int TrggIndex_PFJet320_v3( trigNames.triggerIndex("HLT_PFJet320_v3"));
  unsigned int TrggIndex_PFJet320_v4( trigNames.triggerIndex("HLT_PFJet320_v4"));
  unsigned int TrggIndex_PFJet320_v5( trigNames.triggerIndex("HLT_PFJet320_v5"));
  unsigned int TrggIndex_PFJet320_v6( trigNames.triggerIndex("HLT_PFJet320_v6"));
  unsigned int TrggIndex_PFJet320_v8( trigNames.triggerIndex("HLT_PFJet320_v8"));
  unsigned int TrggIndex_PFJet320_v9( trigNames.triggerIndex("HLT_PFJet320_v9"));
  if(TrggIndex_PFHT650_v5 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v5);
  if(TrggIndex_PFHT650_v6 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v6);
  if(TrggIndex_PFHT650_v7 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v7);
  if(TrggIndex_PFHT650_v8 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v8);
  if(TrggIndex_PFHT650_v9 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFHT650_v9);
  if(TrggIndex_PFNoPUHT650_v1 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v1);
  if(TrggIndex_PFNoPUHT650_v3 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v3);
  if(TrggIndex_PFNoPUHT650_v4 < trigResults->size()) isFired_HLT_HT650 = trigResults->accept(TrggIndex_PFNoPUHT650_v4);
  if(TrggIndex_PFJet320_v3 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v3);
  if(TrggIndex_PFJet320_v4 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v4);
  if(TrggIndex_PFJet320_v5 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v5);
  if(TrggIndex_PFJet320_v6 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v6);
  if(TrggIndex_PFJet320_v8 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v8);
  if(TrggIndex_PFJet320_v9 < trigResults->size()) isFired_HLT_PFJet320 = trigResults->accept(TrggIndex_PFJet320_v9);

  //PILEUP WEIGHT
  double MyWeight = 1;
  if(!isData){
    //edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float Tnpv = -1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
	Tnpv = PVI->getTrueNumInteractions();
	continue;
      }
    }
    MyWeight = LumiWeights_.weight( Tnpv );
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
    if(!(mass>20 && mass<70)) continue;
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
    //if(PFTau->pt()<20) continue;
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
  float ptMuon=-99; bool foundMuon=false;
  reco::MuonCollection::const_iterator SelectedMuon;
  for(reco::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    //if(cktTrack->pt()<10) continue;
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
  
  //MET
  edm::Handle<reco::PFMETCollection> met;
  iEvent.getByLabel("pfMet", met);

  
  //PLOT
  if(foundJet && foundTau && foundMuon && (isFired_HLT_PFJet320==true || isFired_HLT_HT650==true)){
    
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
      m_NVertices=vertices->size();
      m_PUWeight=MyWeight;
      m_jetMass=massZ;
      m_jetPt=SelectedJet->pt();
      m_jetTau21=tau21Z;
      m_tauPt=SelectedTau->pt();
      m_muonPt=SelectedMuon->pt();
      m_metPt=met->begin()->pt();
      m_metPx=met->begin()->px();
      m_metPy=met->begin()->py();
      m_metEta=met->begin()->eta();
      m_metPhi=met->begin()->phi();
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
      
      NVerticesNO->Fill(vertices->size());
      jetMassNO->Fill(massZ);
      jetPtNO->Fill(SelectedJet->pt());
      jetTau21NO->Fill(tau21Z);
      tauPtNO->Fill(SelectedTau->pt());
      muonPtNO->Fill(SelectedMuon->pt());
      metPtNO->Fill(met->begin()->pt());
      metPxNO->Fill(met->begin()->px());
      metPyNO->Fill(met->begin()->py());
      metEtaNO->Fill(met->begin()->eta());
      metPhiNO->Fill(met->begin()->phi());
      jetDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()));
      jetDRMetNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet->p4()));
      jetDRMuoNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4()));
      jetDRTauNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(),SelectedJet->p4()));
      tauDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()));
      tauDRMetNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedTau->p4()));
      tauDRMuoNO->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4()));
      muonDPhiMetNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()));
      muonDRMetNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuon->p4()));
      muonPFIsoNO->Fill(MuonPFIso(SelectedMuon,true));
      ZDRZNO->Fill(ROOT::Math::VectorUtil::DeltaR(SVFitTauTau,SelectedJet->p4()));
      MassSvfitTauMuoNO->Fill(algo.getMass());
      XMassSVFitNO->Fill((SVFitTauTau+PrunedJet).M());
      
      if(met->begin()->pt()>50
	 && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedJet->p4())>2.5 
	 && ROOT::Math::VectorUtil::DeltaR(SelectedTau->p4(), SelectedJet->p4())>2.5
	 && ROOT::Math::VectorUtil::DeltaR(SelectedMuon->p4(),SelectedTau->p4())>0.05
	 && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTau->p4()))<1 
	 && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuon->p4()))<1
	 && fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet->p4()))>2 
	 && SelectedTau->pt()>20 && SelectedMuon->pt()>10 && algo.getMass()>20){
	  
	jetMass->Fill(massZ);
	jetPt->Fill(SelectedJet->pt());
	jetTau21->Fill(tau21Z);
	tauPt->Fill(SelectedTau->pt());
	muonPt->Fill(SelectedMuon->pt());
	metPt->Fill(met->begin()->pt());
	metPx->Fill(met->begin()->px());
	metPy->Fill(met->begin()->py());
	metEta->Fill(met->begin()->eta());
	metPhi->Fill(met->begin()->phi());
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
SidebandMuTauAnalyzer::beginJob()
{
  Service<TFileService> fs;
  TreeVariable = fs->make<TTree>("TreeVariable", "TreeVariable");
  TreeVariable->Branch("NVertices", &m_NVertices, "NVertices/i");
  TreeVariable->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  TreeVariable->Branch("jetMass", &m_jetMass, "jetMass/f");
  TreeVariable->Branch("jetPt", &m_jetPt, "jetPt/f");
  TreeVariable->Branch("jetTau21", &m_jetTau21, "jetTau21/f");
  TreeVariable->Branch("tauPt", &m_tauPt, "tauPt/f");
  TreeVariable->Branch("muonPt", &m_muonPt, "muonPt/f");
  TreeVariable->Branch("metPt", &m_metPt, "metPt/f");
  TreeVariable->Branch("metPx", &m_metPx, "metPx/f");
  TreeVariable->Branch("metPy", &m_metPy, "metPy/f");
  TreeVariable->Branch("metEta", &m_metEta, "metEta/f");
  TreeVariable->Branch("metPhi", &m_metPhi, "metPhi/f");
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
  
  NVertices       = fs->make<TH1D>("NVertices",       "NVertices",       200, 0, 200);
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
  metPx           = fs->make<TH1D>("metPx",           "metPx",           2000, 0, 2000);
  metPy           = fs->make<TH1D>("metPy",           "metPy",           2000, 0, 2000);
  metEta          = fs->make<TH1D>("metEta",          "metEta",          100, -5, 5);
  metPhi          = fs->make<TH1D>("metPhi",          "metPhi",          100, -5, 5);
  ZDRZ            = fs->make<TH1D>("ZDRZ",            "ZDRZ",            100, 0, 5);
  MassSvfitTauMuo = fs->make<TH1D>("MassSvfitTauMuo", "MassSvfitTauMuo", 500, 0, 500);
  XMassSVFit      = fs->make<TH1D>("XMassSVFit",      "XMassSVFit",      300, 0, 3000);
  
  NVerticesNO       = fs->make<TH1D>("NVerticesNO",       "NVerticesNO",       200, 0, 200);
  jetMassNO         = fs->make<TH1D>("jetMassNO",         "jetMassNO",         400, 0, 200);
  jetPtNO           = fs->make<TH1D>("jetPtNO",           "jetPtNO",           2000, 0, 2000);
  jetTau21NO        = fs->make<TH1D>("jetTau21NO",        "jetTau21NO",        100, 0, 1);
  jetDRMetNO        = fs->make<TH1D>("jetDRMetNO",        "jetDRMetNO",        100, 0, 5);
  jetDRTauNO        = fs->make<TH1D>("jetDRTauNO",        "jetDRTauNO",        100, 0, 5);
  jetDRMuoNO        = fs->make<TH1D>("jetDRMuoNO",        "jetDRMuoNO",        100, 0, 5);
  jetDPhiMetNO      = fs->make<TH1D>("jetDPhiMetNO",      "jetDPhiMetNO",      100, -5, 5);
  tauPtNO           = fs->make<TH1D>("tauPtNO",           "tauPtNO",           2000, 0, 2000);
  tauDRMetNO        = fs->make<TH1D>("tauDRMetNO",        "tauDRMetNO",        100, 0, 5);
  tauDRMuoNO        = fs->make<TH1D>("tauDRMuoNO",        "tauDRMuoNO",        100, 0, 5);
  tauDPhiMetNO      = fs->make<TH1D>("tauDPhiMetNO",      "tauDPhiMetNO",      100, -5, 5);
  muonPtNO          = fs->make<TH1D>("muonPtNO",          "muonPtNO",          2000, 0, 2000);
  muonDRMetNO       = fs->make<TH1D>("muonDRMetNO",       "muonDRMetNO",       100, 0, 5);
  muonDPhiMetNO     = fs->make<TH1D>("muonDPhiMetNO",     "muonDPhiMetNO",     100, -5, 5);
  muonPFIsoNO       = fs->make<TH1D>("muonPFIsoNO",       "muonPFIsoNO",       500, 0, 5);
  metPtNO           = fs->make<TH1D>("metPtNO",           "metPtNO",           2000, 0, 2000);
  metPxNO           = fs->make<TH1D>("metPxNO",           "metPxNO",           2000, 0, 2000);
  metPyNO           = fs->make<TH1D>("metPyNO",           "metPyNO",           2000, 0, 2000);
  metEtaNO          = fs->make<TH1D>("metEtaNO",          "metEtaNO",          100, -5, 5);
  metPhiNO          = fs->make<TH1D>("metPhiNO",          "metPhiNO",          100, -5, 5);
  ZDRZNO            = fs->make<TH1D>("ZDRZNO",            "ZDRZNO",            100, 0, 5);
  MassSvfitTauMuoNO = fs->make<TH1D>("MassSvfitTauMuoNO", "MassSvfitTauMuoNO", 500, 0, 500);
  XMassSVFitNO      = fs->make<TH1D>("XMassSVFitNO",      "XMassSVFitNO",      300, 0, 3000);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SidebandMuTauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SidebandMuTauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SidebandMuTauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SidebandMuTauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SidebandMuTauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SidebandMuTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

float SidebandMuTauAnalyzer::MuonPFIso(reco::MuonCollection::const_iterator muon, bool highpt){
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
DEFINE_FWK_MODULE(SidebandMuTauAnalyzer);
