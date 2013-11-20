// -*- C++ -*-
//
// Package:    EleIdStudy
// Class:      EleIdStudy
// 
/**\class EleIdStudy EleIdStudy.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/EleIdStudy.cc

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

//
// class declaration
//

class EleIdStudy : public edm::EDAnalyzer {
public:
  explicit EleIdStudy(const edm::ParameterSet&);
  ~EleIdStudy();
  
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

 
  //HISTOGRAMS
  TH1D* eleTight; TH1D* eleMedium; TH1D* eleLoose;
  TH1D* eleHEEPID;
  TH1D* isoFar1;
  TH1D* isoFar2;
  TH1D* chargedHadronIsoFar1;
  TH1D* chargedHadronIsoFar2;
  TH1D* neutralHadronIsoFar1;
  TH1D* neutralHadronIsoFar2;
  TH1D* photonIsoFar1;
  TH1D* photonIsoFar2;
  TH1D* isoNear1;
  TH1D* isoNear2;
  TH1D* chargedHadronIsoNear1;
  TH1D* chargedHadronIsoNear2;
  TH1D* neutralHadronIsoNear1;
  TH1D* neutralHadronIsoNear2;
  TH1D* photonIsoNear1;
  TH1D* photonIsoNear2;

  TH2D* isoFar1VSpt1;
  TH2D* isoFar2VSpt1;
  TH2D* chargedHadronIsoFar1VSpt1;
  TH2D* chargedHadronIsoFar2VSpt1;
  TH2D* neutralHadronIsoFar1VSpt1;
  TH2D* neutralHadronIsoFar2VSpt1;
  TH2D* photonIsoFar1VSpt1;
  TH2D* photonIsoFar2VSpt1;
  TH2D* isoNear1VSpt1;
  TH2D* isoNear2VSpt1;
  TH2D* chargedHadronIsoNear1VSpt1;
  TH2D* chargedHadronIsoNear2VSpt1;
  TH2D* neutralHadronIsoNear1VSpt1;
  TH2D* neutralHadronIsoNear2VSpt1;
  TH2D* photonIsoNear1VSpt1;
  TH2D* photonIsoNear2VSpt1;
        
  TH2D* isoFar1VSpt2;
  TH2D* isoFar2VSpt2;
  TH2D* chargedHadronIsoFar1VSpt2;
  TH2D* chargedHadronIsoFar2VSpt2;
  TH2D* neutralHadronIsoFar1VSpt2;
  TH2D* neutralHadronIsoFar2VSpt2;
  TH2D* photonIsoFar1VSpt2;
  TH2D* photonIsoFar2VSpt2;
  TH2D* isoNear1VSpt2;
  TH2D* isoNear2VSpt2;
  TH2D* chargedHadronIsoNear1VSpt2;
  TH2D* chargedHadronIsoNear2VSpt2;
  TH2D* neutralHadronIsoNear1VSpt2;
  TH2D* neutralHadronIsoNear2VSpt2;
  TH2D* photonIsoNear1VSpt2;
  TH2D* photonIsoNear2VSpt2;

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
EleIdStudy::EleIdStudy(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  /*
   //MVA IDENDTIFICATION
   // NOTE: It is safer and crab-compliant to get the files locally, i.e in EGamma/EGammaAnalysisTools/data
   // (see the downloard.url file in that directory)
   // Alternatively (for tests), they can be read from AFS:
   std::vector<std::string> myManualCatWeigths;
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
   myManualCatWeigths.push_back("/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");
   
   Bool_t manualCat = true;
   myMVANonTrig = new EGammaMvaEleEstimator();
   myMVANonTrig->initialize("BDT",
			    EGammaMvaEleEstimator::kNonTrig,
			    manualCat, 
			    myManualCatWeigths);
  */
}


EleIdStudy::~EleIdStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EleIdStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel("primaryVertexFilter", vertices);
   reco::Vertex primaryVertex;
   primaryVertex = vertices->at(0);

   edm::Handle<pat::JetCollection> CA8JetswithQjets;
   iEvent.getByLabel("selectedPatJetsCA8CHSwithQjets", CA8JetswithQjets);
   edm::Handle<pat::JetCollection> CA8JetsPruned;
   iEvent.getByLabel("selectedPatJetsCA8CHSpruned", CA8JetsPruned);

   edm::Handle<pat::METCollection> met;
   iEvent.getByLabel("patMETs", met);

   edm::Handle<pat::ElectronCollection> eleH;
   iEvent.getByLabel("patElectronsWithTrigger", eleH);

   edm::Handle<pat::MuonCollection> muoH;
   iEvent.getByLabel("patMuonsWithTrigger", muoH);

   edm::Handle<pat::METCollection> metRaw;
   iEvent.getByLabel("patMETsRaw", metRaw);

   edm::Handle<double> rhoHandle;
   iEvent.getByLabel("kt6PFJets", "rho", rhoHandle);
   float rho = *(rhoHandle.product());
   
   Handle<vector<reco::GenParticle> > genParts;
   iEvent.getByLabel("genParticles", genParts);

   edm::Handle<pat::JetCollection> ak5jetCands;
   iEvent.getByLabel("patJetsWithVarCHS",ak5jetCands);	
  
   int ele = 0; int muo = 0; int neu = 0;
   vector<math::PtEtaPhiELorentzVector> genele;
   for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
     const reco::GenParticle & genPart = (*genParts)[ngenPart];
     if(abs(genPart.pdgId())==15 && genPart.status()!=3){
       for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	 const reco::Candidate * daughter = genPart.daughter(ndaugh);
	 if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	 if(abs(daughter->pdgId())==13 && daughter->status()==1) muo = muo + 1;
	 if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 || abs(daughter->pdgId())==16) && daughter->status()==1) neu = neu + 1;
	 math::PtEtaPhiELorentzVector gen_prov; gen_prov=daughter->p4();
	 if(abs(daughter->pdgId())==11 && daughter->status()==1) genele.push_back(gen_prov);
       }
     }
   }

   int Nele = 0;
   vector<pat::ElectronCollection::const_iterator> SelectedInitialEle;
   SelectElectronCutBased(eleH, SelectedInitialEle, Nele, rho, primaryVertex);

   float DRele1=99.; int eleInd1=-1; float DRele2=99.; int eleInd2=-1;
   for(unsigned int i=0; i<genele.size(); i++){
     for(unsigned int j=0; j<SelectedInitialEle.size(); j++){
       if(i==0 && ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4())<0.2 && ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4())<DRele1){
	 DRele1=ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4());
	 eleInd1=j;
       }
       if(i==1 && ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4())<0.2 && ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4())<DRele2){
	 DRele2=ROOT::Math::VectorUtil::DeltaR(genele[i],SelectedInitialEle[j]->p4());
	 eleInd2=j;
       }
     }
   }
   if(eleInd1!=-1 && eleInd2!=-1){
     if(SelectedInitialEle[eleInd2]->pt()>SelectedInitialEle[eleInd1]->pt()){
       int eleInd3=eleInd1;
       eleInd1=eleInd2;
       eleInd2=eleInd3;
     }
   }



   //LOOSE SELECTION
   if(eleInd1!=-1 && eleInd2!=-1){
     bool firstEle=false; bool secondEle=false;
     eleLoose->Fill(0);
     if(SelectedInitialEle[eleInd1]->pt()>20) eleLoose->Fill(1);
     if(SelectedInitialEle[eleInd1]->pt()>20 && fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=2.5){
       eleLoose->Fill(2);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleLoose->Fill(3);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dz(primaryVertex.position()))<0.2){
	 eleLoose->Fill(4);
	 if(fabs(1/SelectedInitialEle[eleInd1]->ecalEnergy()-SelectedInitialEle[eleInd1]->eSuperClusterOverP()/SelectedInitialEle[eleInd1]->ecalEnergy())<0.05){
	   eleLoose->Fill(5);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0) eleLoose->Fill(6);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0 && SelectedInitialEle[eleInd1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<2){
	     eleLoose->Fill(7);
	     if(fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.007) eleLoose->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.007 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.015){
		 eleLoose->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01) eleLoose->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.12) {
		   eleLoose->Fill(11);
		   firstEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.009) eleLoose->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.009 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.010){
		 eleLoose->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03) eleLoose->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.10) {
		   eleLoose->Fill(11);
		   firstEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20) eleLoose->Fill(12);
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20 && fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=2.5){
       eleLoose->Fill(13);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleLoose->Fill(14);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dz(primaryVertex.position()))<0.2){
	 eleLoose->Fill(15);
	 if(fabs(1/SelectedInitialEle[eleInd2]->ecalEnergy()-SelectedInitialEle[eleInd2]->eSuperClusterOverP()/SelectedInitialEle[eleInd2]->ecalEnergy())<0.05){
	   eleLoose->Fill(16);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0) eleLoose->Fill(17);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0 && SelectedInitialEle[eleInd2]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<2){
	     eleLoose->Fill(18);
	     if(fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.007) eleLoose->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.007 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.015){
		 eleLoose->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01) eleLoose->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.12) {
		   eleLoose->Fill(22);
		   secondEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.009) eleLoose->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.009 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.010){
		 eleLoose->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03) eleLoose->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.10) {
		   eleLoose->Fill(22);
		   secondEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(secondEle){
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.15) eleLoose->Fill(23);
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.15 && ElectronPFIso(SelectedInitialEle[eleInd2],rho)<0.15) eleLoose->Fill(24);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho)) eleLoose->Fill(25);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho) && ElectronDETIso(SelectedInitialEle[eleInd2],rho)) eleLoose->Fill(26);
     }
   }



   //MEDIUM SELECTION
   if(eleInd1!=-1 && eleInd2!=-1){
     bool firstEle=false; bool secondEle=false;
     eleMedium->Fill(0);
     if(SelectedInitialEle[eleInd1]->pt()>20) eleMedium->Fill(1);
     if(SelectedInitialEle[eleInd1]->pt()>20 && fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=2.5){
       eleMedium->Fill(2);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleMedium->Fill(3);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dz(primaryVertex.position()))<0.1){
	 eleMedium->Fill(4);
	 if(fabs(1/SelectedInitialEle[eleInd1]->ecalEnergy()-SelectedInitialEle[eleInd1]->eSuperClusterOverP()/SelectedInitialEle[eleInd1]->ecalEnergy())<0.05){
	   eleMedium->Fill(5);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0) eleMedium->Fill(6);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0 && SelectedInitialEle[eleInd1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<2){
	     eleMedium->Fill(7);
	     if(fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.004) eleMedium->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.004 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.060){
		 eleMedium->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01) eleMedium->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.12) {
		   eleMedium->Fill(11);
		   firstEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.007) eleMedium->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.007 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.030){
		 eleMedium->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03) eleMedium->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.10) {
		   eleMedium->Fill(11);
		   firstEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20) eleMedium->Fill(12);
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20 && fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=2.5){
       eleMedium->Fill(13);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleMedium->Fill(14);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dz(primaryVertex.position()))<0.1){
	 eleMedium->Fill(15);
	 if(fabs(1/SelectedInitialEle[eleInd2]->ecalEnergy()-SelectedInitialEle[eleInd2]->eSuperClusterOverP()/SelectedInitialEle[eleInd2]->ecalEnergy())<0.05){
	   eleMedium->Fill(16);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0) eleMedium->Fill(17);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0 && SelectedInitialEle[eleInd2]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<2){
	     eleMedium->Fill(18);
	     if(fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.004) eleMedium->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.004 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.060){
		 eleMedium->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01) eleMedium->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.12) {
		   eleMedium->Fill(22);
		   secondEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.007) eleMedium->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.007 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.030){
		 eleMedium->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03) eleMedium->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.10) {
		   eleMedium->Fill(22);
		   secondEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(secondEle){
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.15) eleMedium->Fill(23);
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.15 && ElectronPFIso(SelectedInitialEle[eleInd2],rho)<0.15) eleMedium->Fill(24);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho)) eleMedium->Fill(25);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho) && ElectronDETIso(SelectedInitialEle[eleInd2],rho)) eleMedium->Fill(26);
     }
   }




   //TIGHT SELECTION
   if(eleInd1!=-1 && eleInd2!=-1){
     bool firstEle=false; bool secondEle=false;
     eleTight->Fill(0);
     if(SelectedInitialEle[eleInd1]->pt()>20) eleTight->Fill(1);
     if(SelectedInitialEle[eleInd1]->pt()>20 && fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=2.5){
       eleTight->Fill(2);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleTight->Fill(3);
       if(fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd1]->gsfTrack()->dz(primaryVertex.position()))<0.1){
	 eleTight->Fill(4);
	 if(fabs(1/SelectedInitialEle[eleInd1]->ecalEnergy()-SelectedInitialEle[eleInd1]->eSuperClusterOverP()/SelectedInitialEle[eleInd1]->ecalEnergy())<0.05){
	   eleTight->Fill(5);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0) eleTight->Fill(6);
	   if(SelectedInitialEle[eleInd1]->passConversionVeto()!=0 && SelectedInitialEle[eleInd1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0){
	     eleTight->Fill(7);
	     if(fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.004) eleTight->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.004 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.030){
		 eleTight->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01) eleTight->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.12) {
		   eleTight->Fill(11);
		   firstEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.005) eleTight->Fill(8);
	       if(SelectedInitialEle[eleInd1]->deltaEtaSuperClusterTrackAtVtx()<0.005 && SelectedInitialEle[eleInd1]->deltaPhiSuperClusterTrackAtVtx()<0.020){
		 eleTight->Fill(9);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03) eleTight->Fill(10);
		 if(SelectedInitialEle[eleInd1]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd1]->hadronicOverEm()<0.10) {
		   eleTight->Fill(11);
		   firstEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20) eleTight->Fill(12);
     if(firstEle && SelectedInitialEle[eleInd2]->pt()>20 && fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=2.5){
       eleTight->Fill(13);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02) eleTight->Fill(14);
       if(fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dxy(primaryVertex.position()))<0.02 && fabs(SelectedInitialEle[eleInd2]->gsfTrack()->dz(primaryVertex.position()))<0.1){
	 eleTight->Fill(15);
	 if(fabs(1/SelectedInitialEle[eleInd2]->ecalEnergy()-SelectedInitialEle[eleInd2]->eSuperClusterOverP()/SelectedInitialEle[eleInd2]->ecalEnergy())<0.05){
	   eleTight->Fill(16);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0) eleTight->Fill(17);
	   if(SelectedInitialEle[eleInd2]->passConversionVeto()!=0 && SelectedInitialEle[eleInd2]->gsfTrack()->trackerExpectedHitsInner().numberOfHits()==0){
	     eleTight->Fill(18);
	     if(fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<=1.479){
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.004) eleTight->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.004 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.030){
		 eleTight->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01) eleTight->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.01 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.12) {
		   eleTight->Fill(22);
		   secondEle=true;
		 }
	       }
	     }else {
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.005) eleTight->Fill(19);
	       if(SelectedInitialEle[eleInd2]->deltaEtaSuperClusterTrackAtVtx()<0.005 && SelectedInitialEle[eleInd2]->deltaPhiSuperClusterTrackAtVtx()<0.020){
		 eleTight->Fill(20);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03) eleTight->Fill(21);
		 if(SelectedInitialEle[eleInd2]->sigmaIetaIeta()<0.03 && SelectedInitialEle[eleInd2]->hadronicOverEm()<0.10) {
		   eleTight->Fill(22);
		   secondEle=true;
		 }
	       }
	     }
	   }
	 }
       }
     }
     if(secondEle){
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.1) eleTight->Fill(23);
       if(ElectronPFIso(SelectedInitialEle[eleInd1],rho)<0.1 && ElectronPFIso(SelectedInitialEle[eleInd2],rho)<0.1) eleTight->Fill(24);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho)) eleTight->Fill(25);
       if(ElectronDETIso(SelectedInitialEle[eleInd1],rho) && ElectronDETIso(SelectedInitialEle[eleInd2],rho)) eleTight->Fill(26);
     }
   }




   //HEEPID SELECTION
   if(eleInd1!=-1 && eleInd2!=-1){
     bool firstEle=false; bool secondEle=false;
     eleHEEPID->Fill(0);
     if(SelectedInitialEle[eleInd1]->pt()>40) {
       eleHEEPID->Fill(1);
       if(fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<1.4442 || (fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())>1.5666 && fabs(SelectedInitialEle[eleInd1]->superCluster()->eta())<2.5)) {
	 eleHEEPID->Fill(2);
	 if(SelectedInitialEle[eleInd1]->userInt("HEEPId")==0) eleHEEPID->Fill(3);
	 if(SelectedInitialEle[eleInd1]->userInt("HEEPId")==0 && SelectedInitialEle[eleInd2]->pt()>40) {
	   eleHEEPID->Fill(4);
	   if(fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<1.4442 || (fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())>1.5666 && fabs(SelectedInitialEle[eleInd2]->superCluster()->eta())<2.5)) {
	     eleHEEPID->Fill(5);
	     if(SelectedInitialEle[eleInd2]->userInt("HEEPId")==0) eleHEEPID->Fill(6);
	     if(SelectedInitialEle[eleInd2]->userInt("HEEPId")==0 && ElectronDETIso(SelectedInitialEle[eleInd1],rho)) eleHEEPID->Fill(7);
	     if(SelectedInitialEle[eleInd2]->userInt("HEEPId")==0 && ElectronDETIso(SelectedInitialEle[eleInd1],rho) && ElectronDETIso(SelectedInitialEle[eleInd2],rho)) 
	       eleHEEPID->Fill(8);
	   }
	 }
       }
     }
   }

   //ISOLATION STUDY
   if(eleInd1!=-1 && eleInd2!=-1){
     float chargedHadronIso1 = SelectedInitialEle[eleInd1]->chargedHadronIso();
     float neutralHadronIso1 = SelectedInitialEle[eleInd1]->neutralHadronIso();
     float photonIso1 = SelectedInitialEle[eleInd1]->photonIso();
     float thiseta1 = SelectedInitialEle[eleInd1]->superCluster()->eta();
     float Aeff1=0.;
     if(thiseta1<1.0)                    Aeff1=0.13;
     if(thiseta1>=1.0 && thiseta1<1.479) Aeff1=0.14;
     if(thiseta1>=1.479 && thiseta1<2.0) Aeff1=0.07;
     if(thiseta1>=2.0 && thiseta1<2.2)   Aeff1=0.09;
     if(thiseta1>=2.2 && thiseta1<2.3)   Aeff1=0.11;
     if(thiseta1>=2.3 && thiseta1<2.4)   Aeff1=0.11;
     if(thiseta1>=2.4)                   Aeff1=0.14;
     float chargedHadronIso2 = SelectedInitialEle[eleInd2]->chargedHadronIso();
     float neutralHadronIso2 = SelectedInitialEle[eleInd2]->neutralHadronIso();
     float photonIso2 = SelectedInitialEle[eleInd2]->photonIso();
     float thiseta2 = SelectedInitialEle[eleInd2]->superCluster()->eta();
     float Aeff2=0.;
     if(thiseta2<1.0)                    Aeff2=0.13;
     if(thiseta2>=1.0 && thiseta2<1.479) Aeff2=0.14;
     if(thiseta2>=1.479 && thiseta2<2.0) Aeff2=0.07;
     if(thiseta2>=2.0 && thiseta2<2.2)   Aeff2=0.09;
     if(thiseta2>=2.2 && thiseta2<2.3)   Aeff2=0.11;
     if(thiseta2>=2.3 && thiseta2<2.4)   Aeff2=0.11;
     if(thiseta2>=2.4)                   Aeff2=0.14;
     float zero = 0.;
     float iso1 = (chargedHadronIso1 + max(zero, neutralHadronIso1 + photonIso1 - rho*Aeff1))/SelectedInitialEle[eleInd1]->pt();
     float iso2 = (chargedHadronIso2 + max(zero, neutralHadronIso2 + photonIso2 - rho*Aeff2))/SelectedInitialEle[eleInd2]->pt();
     if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[eleInd1]->p4(),SelectedInitialEle[eleInd2]->p4())<0.25 && SelectedInitialEle[eleInd1]->pt()>20 && SelectedInitialEle[eleInd2]->pt()>20){
       isoNear1->Fill(iso1);
       isoNear2->Fill(iso2);
       chargedHadronIsoNear1->Fill(chargedHadronIso1);
       chargedHadronIsoNear2->Fill(chargedHadronIso2);
       neutralHadronIsoNear1->Fill(neutralHadronIso1);
       neutralHadronIsoNear2->Fill(neutralHadronIso2);
       photonIsoNear1->Fill(photonIso1);
       photonIsoNear2->Fill(photonIso2);
       isoNear1VSpt1->Fill(iso1,SelectedInitialEle[eleInd1]->pt());
       isoNear2VSpt1->Fill(iso2,SelectedInitialEle[eleInd1]->pt());
       chargedHadronIsoNear1VSpt1->Fill(chargedHadronIso1,SelectedInitialEle[eleInd1]->pt());
       chargedHadronIsoNear2VSpt1->Fill(chargedHadronIso2,SelectedInitialEle[eleInd1]->pt());
       neutralHadronIsoNear1VSpt1->Fill(neutralHadronIso1,SelectedInitialEle[eleInd1]->pt());
       neutralHadronIsoNear2VSpt1->Fill(neutralHadronIso2,SelectedInitialEle[eleInd1]->pt());
       photonIsoNear1VSpt1->Fill(photonIso1,SelectedInitialEle[eleInd1]->pt());
       photonIsoNear2VSpt1->Fill(photonIso2,SelectedInitialEle[eleInd1]->pt());
       isoNear1VSpt2->Fill(iso1,SelectedInitialEle[eleInd2]->pt());
       isoNear2VSpt2->Fill(iso2,SelectedInitialEle[eleInd2]->pt());
       chargedHadronIsoNear1VSpt2->Fill(chargedHadronIso1,SelectedInitialEle[eleInd2]->pt());
       chargedHadronIsoNear2VSpt2->Fill(chargedHadronIso2,SelectedInitialEle[eleInd2]->pt());
       neutralHadronIsoNear1VSpt2->Fill(neutralHadronIso1,SelectedInitialEle[eleInd2]->pt());
       neutralHadronIsoNear2VSpt2->Fill(neutralHadronIso2,SelectedInitialEle[eleInd2]->pt());
       photonIsoNear1VSpt2->Fill(photonIso1,SelectedInitialEle[eleInd2]->pt());
       photonIsoNear2VSpt2->Fill(photonIso2,SelectedInitialEle[eleInd2]->pt());
     }
     if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[eleInd1]->p4(),SelectedInitialEle[eleInd2]->p4())>0.50 && SelectedInitialEle[eleInd1]->pt()>20 && SelectedInitialEle[eleInd2]->pt()>20){
       isoFar1->Fill(iso1);
       isoFar2->Fill(iso2);
       chargedHadronIsoFar1->Fill(chargedHadronIso1);
       chargedHadronIsoFar2->Fill(chargedHadronIso2);
       neutralHadronIsoFar1->Fill(neutralHadronIso1);
       neutralHadronIsoFar2->Fill(neutralHadronIso2);
       photonIsoFar1->Fill(photonIso1);
       photonIsoFar2->Fill(photonIso2);
       isoFar1VSpt1->Fill(iso1,SelectedInitialEle[eleInd1]->pt());
       isoFar2VSpt1->Fill(iso2,SelectedInitialEle[eleInd1]->pt());
       chargedHadronIsoFar1VSpt1->Fill(chargedHadronIso1,SelectedInitialEle[eleInd1]->pt());
       chargedHadronIsoFar2VSpt1->Fill(chargedHadronIso2,SelectedInitialEle[eleInd1]->pt());
       neutralHadronIsoFar1VSpt1->Fill(neutralHadronIso1,SelectedInitialEle[eleInd1]->pt());
       neutralHadronIsoFar2VSpt1->Fill(neutralHadronIso2,SelectedInitialEle[eleInd1]->pt());
       photonIsoFar1VSpt1->Fill(photonIso1,SelectedInitialEle[eleInd1]->pt());
       photonIsoFar2VSpt1->Fill(photonIso2,SelectedInitialEle[eleInd1]->pt());
       isoFar1VSpt2->Fill(iso1,SelectedInitialEle[eleInd2]->pt());
       isoFar2VSpt2->Fill(iso2,SelectedInitialEle[eleInd2]->pt());
       chargedHadronIsoFar1VSpt2->Fill(chargedHadronIso1,SelectedInitialEle[eleInd2]->pt());
       chargedHadronIsoFar2VSpt2->Fill(chargedHadronIso2,SelectedInitialEle[eleInd2]->pt());
       neutralHadronIsoFar1VSpt2->Fill(neutralHadronIso1,SelectedInitialEle[eleInd2]->pt());
       neutralHadronIsoFar2VSpt2->Fill(neutralHadronIso2,SelectedInitialEle[eleInd2]->pt());
       photonIsoFar1VSpt2->Fill(photonIso1,SelectedInitialEle[eleInd2]->pt());
       photonIsoFar2VSpt2->Fill(photonIso2,SelectedInitialEle[eleInd2]->pt());
     }
   }
   
   //MVA IDENTIFICATION
   //double myMVANonTrigMethod1 = myMVANonTrig->mvaValue(electron,*primaryVertex,thebuilder,lazyTools,debugMVAclass);

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
EleIdStudy::beginJob()
{
  Service<TFileService> fs;

  eleLoose  = fs->make<TH1D>("eleLoose",  "eleLoose",  30, -0.5, 29.5);
  eleMedium = fs->make<TH1D>("eleMedium", "eleMedium", 30, -0.5, 29.5);
  eleTight  = fs->make<TH1D>("eleTight",  "eleTight",  30, -0.5, 29.5);
  eleHEEPID = fs->make<TH1D>("eleHEEPID", "eleHEEPID", 30, -0.5, 29.5);

  isoFar1  = fs->make<TH1D>("isoFar1",  "isoFar1",  1000, 0, 10);
  isoFar2  = fs->make<TH1D>("isoFar2",  "isoFar2",  1000, 0, 10);
  chargedHadronIsoFar1  = fs->make<TH1D>("chargedHadronIsoFar1",  "chargedHadronIsoFar1",  10000, 0, 1000);
  chargedHadronIsoFar2  = fs->make<TH1D>("chargedHadronIsoFar2",  "chargedHadronIsoFar2",  10000, 0, 1000);
  neutralHadronIsoFar1  = fs->make<TH1D>("neutralHadronIsoFar1",  "neutralHadronIsoFar1",  10000, 0, 1000);
  neutralHadronIsoFar2  = fs->make<TH1D>("neutralHadronIsoFar2",  "neutralHadronIsoFar2",  10000, 0, 1000);
  photonIsoFar1  = fs->make<TH1D>("photonIsoFar1",  "photonIsoFar1", 10000, 0, 1000);
  photonIsoFar2  = fs->make<TH1D>("photonIsoFar2",  "photonIsoFar2", 10000, 0, 1000);
  isoNear1  = fs->make<TH1D>("isoNear1",  "isoNear1",  1000, 0, 10);
  isoNear2  = fs->make<TH1D>("isoNear2",  "isoNear2",  1000, 0, 10);
  chargedHadronIsoNear1  = fs->make<TH1D>("chargedHadronIsoNear1",  "chargedHadronIsoNear1",  10000, 0, 1000);
  chargedHadronIsoNear2  = fs->make<TH1D>("chargedHadronIsoNear2",  "chargedHadronIsoNear2",  10000, 0, 1000);
  neutralHadronIsoNear1  = fs->make<TH1D>("neutralHadronIsoNear1",  "neutralHadronIsoNear1",  10000, 0, 1000);
  neutralHadronIsoNear2  = fs->make<TH1D>("neutralHadronIsoNear2",  "neutralHadronIsoNear2",  10000, 0, 1000);
  photonIsoNear1  = fs->make<TH1D>("photonIsoNear1",  "photonIsoNear1", 10000, 0, 1000);
  photonIsoNear2  = fs->make<TH1D>("photonIsoNear2",  "photonIsoNear2", 10000, 0, 1000);

  isoFar1VSpt1  = fs->make<TH2D>("isoFar1VSpt1",  "isoFar1VSpt1",  50, 0, 5,  500, 0, 500);
  isoFar2VSpt1  = fs->make<TH2D>("isoFar2VSpt1",  "isoFar2VSpt1",  50, 0, 5,  500, 0, 500);
  chargedHadronIsoFar1VSpt1  = fs->make<TH2D>("chargedHadronIsoFar1VSpt1",  "chargedHadronIsoFar1VSpt1",  500, 0, 500,  500, 0, 500);
  chargedHadronIsoFar2VSpt1  = fs->make<TH2D>("chargedHadronIsoFar2VSpt1",  "chargedHadronIsoFar2VSpt1",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoFar1VSpt1  = fs->make<TH2D>("neutralHadronIsoFar1VSpt1",  "neutralHadronIsoFar1VSpt1",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoFar2VSpt1  = fs->make<TH2D>("neutralHadronIsoFar2VSpt1",  "neutralHadronIsoFar2VSpt1",  500, 0, 500,  500, 0, 500);
  photonIsoFar1VSpt1  = fs->make<TH2D>("photonIsoFar1VSpt1",  "photonIsoFar1VSpt1", 500, 0, 500,  500, 0, 500);
  photonIsoFar2VSpt1  = fs->make<TH2D>("photonIsoFar2VSpt1",  "photonIsoFar2VSpt1", 500, 0, 500,  500, 0, 500);
  isoNear1VSpt1  = fs->make<TH2D>("isoNear1VSpt1",  "isoNear1VSpt1",  50, 0, 5,  500, 0, 500);
  isoNear2VSpt1  = fs->make<TH2D>("isoNear2VSpt1",  "isoNear2VSpt1",  50, 0, 5,  500, 0, 500);
  chargedHadronIsoNear1VSpt1  = fs->make<TH2D>("chargedHadronIsoNear1VSpt1",  "chargedHadronIsoNear1VSpt1",  500, 0, 500,  500, 0, 500);
  chargedHadronIsoNear2VSpt1  = fs->make<TH2D>("chargedHadronIsoNear2VSpt1",  "chargedHadronIsoNear2VSpt1",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoNear1VSpt1  = fs->make<TH2D>("neutralHadronIsoNear1VSpt1",  "neutralHadronIsoNear1VSpt1",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoNear2VSpt1  = fs->make<TH2D>("neutralHadronIsoNear2VSpt1",  "neutralHadronIsoNear2VSpt1",  500, 0, 500,  500, 0, 500);
  photonIsoNear1VSpt1  = fs->make<TH2D>("photonIsoNear1VSpt1",  "photonIsoNear1VSpt1", 500, 0, 500,  500, 0, 500);
  photonIsoNear2VSpt1  = fs->make<TH2D>("photonIsoNear2VSpt1",  "photonIsoNear2VSpt1", 500, 0, 500,  500, 0, 500);

  isoFar1VSpt2  = fs->make<TH2D>("isoFar1VSpt2",  "isoFar1VSpt2",  50, 0, 5,  500, 0, 500);
  isoFar2VSpt2  = fs->make<TH2D>("isoFar2VSpt2",  "isoFar2VSpt2",  50, 0, 5,  500, 0, 500);
  chargedHadronIsoFar1VSpt2  = fs->make<TH2D>("chargedHadronIsoFar1VSpt2",  "chargedHadronIsoFar1VSpt2",  500, 0, 500,  500, 0, 500);
  chargedHadronIsoFar2VSpt2  = fs->make<TH2D>("chargedHadronIsoFar2VSpt2",  "chargedHadronIsoFar2VSpt2",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoFar1VSpt2  = fs->make<TH2D>("neutralHadronIsoFar1VSpt2",  "neutralHadronIsoFar1VSpt2",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoFar2VSpt2  = fs->make<TH2D>("neutralHadronIsoFar2VSpt2",  "neutralHadronIsoFar2VSpt2",  500, 0, 500,  500, 0, 500);
  photonIsoFar1VSpt2  = fs->make<TH2D>("photonIsoFar1VSpt2",  "photonIsoFar1VSpt2", 500, 0, 500,  500, 0, 500);
  photonIsoFar2VSpt2  = fs->make<TH2D>("photonIsoFar2VSpt2",  "photonIsoFar2VSpt2", 500, 0, 500,  500, 0, 500);
  isoNear1VSpt2  = fs->make<TH2D>("isoNear1VSpt2",  "isoNear1VSpt2",  50, 0, 5,  500, 0, 500);
  isoNear2VSpt2  = fs->make<TH2D>("isoNear2VSpt2",  "isoNear2VSpt2",  50, 0, 5,  500, 0, 500);
  chargedHadronIsoNear1VSpt2  = fs->make<TH2D>("chargedHadronIsoNear1VSpt2",  "chargedHadronIsoNear1VSpt2",  500, 0, 500,  500, 0, 500);
  chargedHadronIsoNear2VSpt2  = fs->make<TH2D>("chargedHadronIsoNear2VSpt2",  "chargedHadronIsoNear2VSpt2",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoNear1VSpt2  = fs->make<TH2D>("neutralHadronIsoNear1VSpt2",  "neutralHadronIsoNear1VSpt2",  500, 0, 500,  500, 0, 500);
  neutralHadronIsoNear2VSpt2  = fs->make<TH2D>("neutralHadronIsoNear2VSpt2",  "neutralHadronIsoNear2VSpt2",  500, 0, 500,  500, 0, 500);
  photonIsoNear1VSpt2  = fs->make<TH2D>("photonIsoNear1VSpt2",  "photonIsoNear1VSpt2", 500, 0, 500,  500, 0, 500);
  photonIsoNear2VSpt2  = fs->make<TH2D>("photonIsoNear2VSpt2",  "photonIsoNear2VSpt2", 500, 0, 500,  500, 0, 500);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EleIdStudy::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
EleIdStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EleIdStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EleIdStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EleIdStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EleIdStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void EleIdStudy::SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex){
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    SelectedEle.push_back(electron);
  }
}


float EleIdStudy::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
  float chargedHadronIso = electron->chargedHadronIso();
  float neutralHadronIso = electron->neutralHadronIso();
  float photonIso = electron->photonIso();
  float thiseta = electron->superCluster()->eta();
  float Aeff=0.;
  if(thiseta<1.0)                   Aeff=0.13;
  if(thiseta>=1.0 && thiseta<1.479) Aeff=0.14;
  if(thiseta>=1.479 && thiseta<2.0) Aeff=0.07;
  if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.09;
  if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.11;
  if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.11;
  if(thiseta>=2.4)                  Aeff=0.14;
  float zero = 0.;
  float iso = (chargedHadronIso + max(zero, neutralHadronIso + photonIso - rho*Aeff))/electron->pt();
  return iso;
}

bool EleIdStudy::ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho){
  bool iso = false;
  if(electron->userIso(0) < 5.0){
    bool inBarrel = electron->isEE();
    double ECALIsol = electron->userIso(1);
    double HCALIsol = electron->userIso(2);
    double sumCaloEt = ECALIsol + HCALIsol;
    double sumCaloEtLimit = -1;
    double et = electron->et();
    if(inBarrel) sumCaloEtLimit = 2.0 + 0.03*et + 0.28*rho;
    if(!inBarrel){
      if(et<50)  sumCaloEtLimit = 2.5 + 0.28*rho;
      else       sumCaloEtLimit = 2.5 + 0.03*(et-50) + 0.28*rho;
    }
    if(sumCaloEt<sumCaloEtLimit) iso = true;
  }
  return iso;
}


//define this as a plug-in
DEFINE_FWK_MODULE(EleIdStudy);
