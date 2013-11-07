// -*- C++ -*-
//
// Package:    SidebandFullyLeptonicAnalyzer
// Class:      SidebandFullyLeptonicAnalyzer
// 
/**\class SidebandFullyLeptonicAnalyzer SidebandFullyLeptonicAnalyzer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/SidebandFullyLeptonicAnalyzer.cc

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

class SidebandFullyLeptonicAnalyzer : public edm::EDAnalyzer {
public:
  explicit SidebandFullyLeptonicAnalyzer(const edm::ParameterSet&);
  ~SidebandFullyLeptonicAnalyzer();
  
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
  void SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets, edm::Handle<pat::JetCollection> CA8JetsPruned, 
		 std::vector<pat::JetCollection::const_iterator> & SelectedJet, std::vector<float> & SelectedPrunedMass, int & Njet);
  void SelectElectron(edm::Handle<pat::ElectronCollection> eleH, std::vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho);
  void SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, std::vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex);
  void SelectTightMuon(  edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectHighptMuon( edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, std::vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex);
  void SelectTrackerGlobalID(std::vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo);
  float ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho);
  bool  ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho);
  float MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt);
  float MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, std::vector<pat::MuonCollection::const_iterator> SelectedHighptMuo);
  void BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, pat::JetCollection::const_iterator SelectedJet);
  float SVFitMass(edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> metRaw, TLorentzVector & SVFitTauTau, LorentzVector lep1, LorentzVector lep2);
  void CollinearApproximation(edm::Handle<pat::METCollection> met, TLorentzVector & CATauTau, TLorentzVector PrunedJet, LorentzVector lep1, LorentzVector lep2, bool & CA);
  void jetPlot(edm::Handle<pat::METCollection> met, std::vector<pat::JetCollection::const_iterator> SelectedJet, std::vector<float> SelectedPrunedMass, 
	       std::vector<pat::MuonCollection::const_iterator> SelectedMuo, std::vector<pat::ElectronCollection::const_iterator> SelectedEle, int jetInd);
  void muonPlot(edm::Handle<pat::METCollection> met, std::vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex);
  void electronPlot(edm::Handle<pat::METCollection> met, std::vector<pat::ElectronCollection::const_iterator> SelectedEle, float rho);
  void massTauPlot(edm::Handle<pat::METCollection> met, std::vector<pat::ElectronCollection::const_iterator> SelectedEle, 
		   std::vector<pat::MuonCollection::const_iterator> SelectedMuo, float MassSVFit, float MassCA);
  void massXPlot(edm::Handle<pat::METCollection> met, std::vector<pat::ElectronCollection::const_iterator> SelectedEle, 
		 std::vector<pat::MuonCollection::const_iterator> SelectedMuo, math::PtEtaPhiMLorentzVector PrunedJet_prov, float XmassCA, float XMassSVFit);
  void otherPlot(edm::Handle<pat::METCollection> met, int Njet, int Nele, int Nmuo, float dRJetZ, int nbtagsL, int nbtagsM, int nbtagsT, bool selected);

 

  //HISTOGRAMS
  TH1D* metPt; TH1D* NJet; TH1D* NEle; TH1D* NMuo; TH1D* NLep; 
  TH1D* metPtSelected; TH1D* NJetSelected; TH1D* NEleSelected; TH1D* NMuoSelected; TH1D* NLepSelected;
  TH1D* NbtagsLSelected; TH1D* NbtagsMSelected; TH1D* NbtagsTSelected;
  TH1D* JetZdRSelected; TH1D* tauDecay; TH1D* tauDecaySelected; TH1D* EleMuoDRSelected;

  TH1D* MassVisEleEle; TH1D* MassEffEleEle; TH1D* MassVisMuoMuo; 
  TH1D* MassEffMuoMuo; TH1D* MassVisEleMuo; TH1D* MassEffEleMuo;
  TH1D* MassSvfitEleEle; TH1D* MassSvfitMuoMuo; TH1D* MassSvfitEleMuo; 
  TH1D* MassCAEleEle; TH1D* MassCAMuoMuo;TH1D* MassCAEleMuo; 
  TH1D* XMassVisEleEle; TH1D* XMassVisMuoMuo; TH1D* XMassVisEleMuo;
  TH1D* XMassEffEleEle; TH1D* XMassEffMuoMuo; TH1D* XMassEffEleMuo;
  TH1D* XMassSVFitEleEle; TH1D* XMassSVFitMuoMuo; TH1D* XMassSVFitEleMuo;
  TH1D* XMassCAEleEle; TH1D* XMassCAMuoMuo; TH1D* XMassCAEleMuo;

  TH1D* jetPtSelected; TH1D* jetEtaSelected; TH1D* jetMassSelected; TH1D* jetSubjettinessSelected; 
  TH1D* jetDPhiMetSelected; TH1D* jetDRMetSelected; TH1D* jetDRLepSelected; TH1D* jetDRVisZSelected; 

  TH1D* electronPtLead; TH1D* electronEtaLead; TH1D* electronDRMetLead; 
  TH1D* electronPFIsoLead; TH1D* electronDetIsoLead; TH1D* electronDPhiMetLead;
  TH1D* electronPtSublead; TH1D* electronEtaSublead; TH1D* electronDRMetSublead; 
  TH1D* electronPFIsoSublead; TH1D* electronDetIsoSublead; TH1D* electronDPhiMetSublead;
  TH1D* electronPtEMU; TH1D* electronEtaEMU; TH1D* electronDRMetEMU; TH1D* electronDR; 
  TH1D* electronPFIsoEMU; TH1D* electronDetIsoEMU; TH1D* electronDPhiMetEMU;

  TH1D* muonPtLead; TH1D* muonEtaLead; TH1D* muonDRMetLead; TH1D* muonDPhiMetLead; TH1D* muonPFIsoLead; TH1D* muonDetIsoLead; TH1D* muonDR; 
  TH1D* muonPtSublead; TH1D* muonEtaSublead; TH1D* muonDRMetSublead; TH1D* muonDPhiMetSublead; TH1D* muonPFIsoSublead; TH1D* muonDetIsoSublead; 
  TH1D* muonPtEMU; TH1D* muonEtaEMU; TH1D* muonDRMetEMU; TH1D* muonDPhiMetEMU; TH1D* muonPFIsoEMU; TH1D* muonDetIsoEMU;

  TH1D* metPtSelectedNO;
  TH1D* muonDRMetLeadNO;   TH1D* muonDPhiMetLeadNO;   TH1D* muonDRMetSubleadNO;   TH1D* muonDPhiMetSubleadNO;
  TH1D* electronDRMetLeadNO;   TH1D* electronDPhiMetLeadNO;   TH1D* electronDRMetSubleadNO;   TH1D* electronDPhiMetSubleadNO; 
  TH1D* muonDRMetEMUNO;   TH1D* muonDPhiMetEMUNO;   TH1D* electronDRMetEMUNO;   TH1D* electronDPhiMetEMUNO;
  TH1D* jetPtNO; TH1D* jetSubjettinessNO; TH1D* jetMassNO;

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
SidebandFullyLeptonicAnalyzer::SidebandFullyLeptonicAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


SidebandFullyLeptonicAnalyzer::~SidebandFullyLeptonicAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SidebandFullyLeptonicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   edm::Handle<pat::JetCollection> ak5jetCands;
   iEvent.getByLabel("patJetsWithVarCHS",ak5jetCands);	

   int Njet = 0;
   vector<pat::JetCollection::const_iterator> SelectedJet;
   vector<float> SelectedPrunedMass;
   SelectJet(CA8JetswithQjets, CA8JetsPruned, SelectedJet, SelectedPrunedMass, Njet);
   float jetpt=-99; int jetInd = -1;
   for(unsigned int j=0; j<SelectedJet.size(); j++){
     if(SelectedJet[j]->pt()>jetpt) {
       jetInd = j;
       jetpt = SelectedJet[j]->pt();
     }
   }

   int Nele = 0;
   vector<pat::ElectronCollection::const_iterator> SelectedInitialEle;
   vector<pat::ElectronCollection::const_iterator> SelectedEle;
   SelectElectronCutBased(eleH, SelectedInitialEle, Nele, rho, primaryVertex);
   for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
     if(jetInd==-1) continue;
     if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
     if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
     if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
     SelectedEle.push_back(SelectedInitialEle[i]);
   }
   
   int Nmuo = 0;
   vector<pat::MuonCollection::const_iterator> SelectedTrackerMuo;
   vector<pat::MuonCollection::const_iterator> SelectedHighptMuo;
   vector<pat::MuonCollection::const_iterator> SelectedMuo;
   SelectHighptMuon( muoH, SelectedHighptMuo,  Nmuo, primaryVertex);
   SelectTrackerMuon(muoH, SelectedTrackerMuo, Nmuo, primaryVertex);
   for(unsigned int i=0; i<SelectedTrackerMuo.size(); i++){
     if(jetInd==-1) continue;
     if(ROOT::Math::VectorUtil::DeltaR(SelectedTrackerMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
     if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedTrackerMuo[i]->p4())>1) continue;
     if((MuonDETIso(SelectedTrackerMuo[i],SelectedTrackerMuo))>0.2) continue;
     SelectedMuo.push_back(SelectedTrackerMuo[i]);
   }
   bool hasAtLeastOneHighPtMuo = false;
   SelectTrackerGlobalID(SelectedMuo, primaryVertex, hasAtLeastOneHighPtMuo);
   
   int nbtagsL=0; int nbtagsM=0; int nbtagsT=0;
   if(jetInd!=-1) BtagVeto(ak5jetCands, nbtagsL, nbtagsM, nbtagsT, SelectedJet[jetInd]);
   otherPlot(met, Njet, Nele, Nmuo, 0, nbtagsL, nbtagsM, nbtagsT, false);

   //DI-ELECTRON SELECTION
   if(SelectedEle.size()==2 && SelectedJet.size()>0){
     math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
     TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

     //SVFIT
     TLorentzVector SVFitTauTau; 
     float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedEle[1]->p4());
     float XMassSVFit = (SVFitTauTau+PrunedJet).M();
 
     //COLLINEAR APPROXIMATION
     TLorentzVector CATauTau; bool CA = false;
     CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedEle[1]->p4(), CA);
     float MassCA = 0;
     float XmassCA = 0.;
     float dRJetZ  = 0.;
     if(CA){
       MassCA = CATauTau.M();
       XmassCA = (CATauTau+PrunedJet).M();
       dRJetZ = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
     }
     float dRJetLep1 = ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet[jetInd]->p4());
     float dRJetLep2 = ROOT::Math::VectorUtil::DeltaR(SelectedEle[1]->p4(),SelectedJet[jetInd]->p4());
     float dRLepLep  = ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedEle[1]->p4());
     float dPhiMetLep1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
     float dPhiMetLep2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()));
     math::PtEtaPhiELorentzVector dilep; dilep = SelectedEle[0]->p4()+SelectedEle[1]->p4();

     //PLOT
     if(met->begin()->pt()>100 && dRJetLep1>0.8 && dRJetLep2>0.8 && dRLepLep>0.02 && dPhiMetLep1<1 && dPhiMetLep2<1 && dilep.M()<70 && nbtagsM==0
	&& SelectedEle[0]->pt()>10 && SelectedEle[1]->pt()>10 && SelectedJet[jetInd]->pt()>400){
       electronPlot(met, SelectedEle, rho);
       jetPlot(met, SelectedJet, SelectedPrunedMass, SelectedMuo, SelectedEle, jetInd);
       massTauPlot(met,SelectedEle, SelectedMuo, MassSVFit, MassCA);
       massXPlot(met, SelectedEle, SelectedMuo, PrunedJet_prov, XmassCA, XMassSVFit);
       otherPlot(met, Njet, Nele, Nmuo, dRJetZ, nbtagsL, nbtagsM, nbtagsT, true);
     }
     //CONTROL PLOTS
     metPtSelectedNO->Fill(met->begin()->pt());
     if(met->begin()->pt()>100){
       electronDRMetLeadNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[0]->p4()));
       electronDPhiMetLeadNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
       electronDRMetSubleadNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[1]->p4()));
       electronDPhiMetSubleadNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()));
     }
   }
   //CONTROL PLOTS
   if(SelectedEle.size()==2){
     pat::JetCollection::const_iterator SelectedJet1;
     pat::JetCollection::const_iterator SelectedJet2;
     pat::JetCollection::const_iterator SelectedJet3;
     float ptMax1=-999; float ptMax2=-999; float ptMax3=-999;
     float massSelected = 0.;
     for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
       float dRmin = 9999.; float mass = 0.;
       for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
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
       if(abs(jet->eta())>2.4) continue;
       if(jet->pt()>400 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax1){
	   ptMax1=jet->pt();
	   SelectedJet1=jet;
	 }
       }
       if((jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax2){
	   ptMax2=jet->pt();
	   SelectedJet2=jet;
	 }
       }
       if(jet->pt()>400 && (jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75){
	 if(jet->pt()>ptMax3){
	   ptMax3=jet->pt();
	   massSelected=mass;
	   SelectedJet3=jet;
	 }
       }
     }
     if(ptMax1>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[1]->p4(), SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedEle[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet1->p4()))>2 ) jetSubjettinessNO->Fill(SelectedJet1->userFloat("tau2")/SelectedJet1->userFloat("tau1")>0.75);
     if(ptMax2>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet2->p4())>0.8 
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[1]->p4(), SelectedJet2->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedEle[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet2->p4()))>2 ) jetPtNO->Fill(SelectedJet2->pt());
     if(ptMax3>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[1]->p4(), SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedEle[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet3->p4()))>2 ) jetMassNO->Fill(massSelected);
   }





   //DI-MUON SELECTION   
   if(SelectedMuo.size()==2  && hasAtLeastOneHighPtMuo == true && SelectedJet.size()>0){
     math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
     TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

     //SVFIT
     TLorentzVector SVFitTauTau; 
     float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedMuo[0]->p4(), SelectedMuo[1]->p4());
     float XMassSVFit = (SVFitTauTau+PrunedJet).M();

     //COLLINEAR APPROXIMATION
     TLorentzVector CATauTau; bool CA = false;
     CollinearApproximation(met, CATauTau, PrunedJet, SelectedMuo[0]->p4(), SelectedMuo[1]->p4(), CA);
     float MassCA = 0;
     float XmassCA = 0.;
     float dRJetZ  = 0.;
     if(CA){
       MassCA = CATauTau.M();
       XmassCA = (CATauTau+PrunedJet).M();
       dRJetZ = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
     }
     float dRJetLep1 = ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet[jetInd]->p4());
     float dRJetLep2 = ROOT::Math::VectorUtil::DeltaR(SelectedMuo[1]->p4(),SelectedJet[jetInd]->p4());
     float dRLepLep  = ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedMuo[1]->p4());
     float dPhiMetLep1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
     float dPhiMetLep2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()));
     math::PtEtaPhiELorentzVector dilep; dilep = SelectedMuo[0]->p4()+SelectedMuo[1]->p4();

     //PLOT
     if(met->begin()->pt()>100 && dRJetLep1>0.8 && dRJetLep2>0.8 && dRLepLep>0.02 && dPhiMetLep1<1 && dPhiMetLep2<1 && dilep.M()<70 && nbtagsM==0
	&& SelectedMuo[0]->pt()>10 && SelectedMuo[1]->pt()>10 && SelectedJet[jetInd]->pt()>400){
       muonPlot(met, SelectedMuo, primaryVertex);
       jetPlot(met, SelectedJet, SelectedPrunedMass, SelectedMuo, SelectedEle, jetInd);
       massTauPlot(met, SelectedEle, SelectedMuo, MassSVFit, MassCA);
       massXPlot(met, SelectedEle, SelectedMuo, PrunedJet_prov, XmassCA, XMassSVFit);
       otherPlot(met, Njet, Nele, Nmuo, dRJetZ, nbtagsL, nbtagsM, nbtagsT, true);
     }
     //CONTROL PLOTS
     metPtSelectedNO->Fill(met->begin()->pt());
     if(met->begin()->pt()>100){
       muonDRMetLeadNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[0]->p4()));
       muonDPhiMetLeadNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
       muonDRMetSubleadNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[1]->p4()));
       muonDPhiMetSubleadNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()));
     }
   }
   //CONTROL PLOTS
   if(SelectedMuo.size()==2  && hasAtLeastOneHighPtMuo == true){
     pat::JetCollection::const_iterator SelectedJet1;
     pat::JetCollection::const_iterator SelectedJet2;
     pat::JetCollection::const_iterator SelectedJet3;
     float ptMax1=-999; float ptMax2=-999; float ptMax3=-999;
     float massSelected = 0.;
     for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
       float dRmin = 9999.; float mass = 0.;
       for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
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
       if(abs(jet->eta())>2.4) continue;
       if(jet->pt()>400 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax1){
	   ptMax1=jet->pt();
	   SelectedJet1=jet;
	 }
       }
       if((jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax2){
	   ptMax2=jet->pt();
	   SelectedJet2=jet;
	 }
       }
       if(jet->pt()>400 && (jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75){
	 if(jet->pt()>ptMax3){
	   ptMax3=jet->pt();
	   massSelected=mass;
	   SelectedJet3=jet;
	 }
       }
     }
     if(ptMax1>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[1]->p4(), SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedMuo[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet1->p4()))>2 ) jetSubjettinessNO->Fill(SelectedJet1->userFloat("tau2")/SelectedJet1->userFloat("tau1")>0.75);
     if(ptMax2>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet2->p4())>0.8 
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[1]->p4(), SelectedJet2->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedMuo[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet2->p4()))>2 ) jetPtNO->Fill(SelectedJet2->pt());
     if(ptMax3>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[1]->p4(), SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedMuo[1]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet3->p4()))>2 ) jetMassNO->Fill(massSelected);
   }



   //ELECTRON-MUON SELECTION
   SelectedEle.clear(); SelectedMuo.clear();
   for(unsigned int i=0; i<SelectedInitialEle.size(); i++){
     if(jetInd==-1) continue;
     if(ROOT::Math::VectorUtil::DeltaR(SelectedInitialEle[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
     if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedInitialEle[i]->p4())>1) continue;
     if(ElectronPFIso(SelectedInitialEle[i],rho)>0.1) continue;
     SelectedEle.push_back(SelectedInitialEle[i]);
   }
   for(unsigned int i=0; i<SelectedHighptMuo.size(); i++){
     if(jetInd==-1) continue;
     if(ROOT::Math::VectorUtil::DeltaR(SelectedHighptMuo[i]->p4(),SelectedJet[jetInd]->p4())<0.8) continue;
     if(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedHighptMuo[i]->p4())>1) continue;
     if(MuonPFIso(SelectedHighptMuo[i], true)>0.2) continue;
     SelectedMuo.push_back(SelectedHighptMuo[i]);
   }

   if(SelectedEle.size()==1 && SelectedMuo.size()==1 && SelectedJet.size()>0){
     math::PtEtaPhiMLorentzVector PrunedJet_prov(SelectedJet[jetInd]->pt(),SelectedJet[jetInd]->eta(),SelectedJet[jetInd]->phi(),SelectedPrunedMass[jetInd]);
     TLorentzVector PrunedJet; PrunedJet.SetPxPyPzE(PrunedJet_prov.px(),PrunedJet_prov.py(),PrunedJet_prov.pz(),PrunedJet_prov.E());

     //SVFIT
     TLorentzVector SVFitTauTau; 
     float MassSVFit = SVFitMass(met, metRaw, SVFitTauTau, SelectedEle[0]->p4(), SelectedMuo[0]->p4()); 
     float XMassSVFit = (SVFitTauTau+PrunedJet).M();

     //COLLINEAR APPROXIMATION
     TLorentzVector CATauTau; bool CA = false;
     CollinearApproximation(met, CATauTau, PrunedJet, SelectedEle[0]->p4(), SelectedMuo[0]->p4(), CA);
     float MassCA = 0;
     float XmassCA = 0.;
     float dRJetZ  = 0.;
     if(CA){
       MassCA = CATauTau.M();
       XmassCA = (CATauTau+PrunedJet).M();
       dRJetZ = ROOT::Math::VectorUtil::DeltaR(CATauTau,PrunedJet);
     }
     float dRJetLep1 = ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet[jetInd]->p4());
     float dRJetLep2 = ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet[jetInd]->p4());
     float dRLepLep  = ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedMuo[0]->p4());
     float dPhiMetLep1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
     float dPhiMetLep2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
     math::PtEtaPhiELorentzVector dilep; dilep = SelectedEle[0]->p4()+SelectedMuo[0]->p4();

     //PLOT
     if(met->begin()->pt()>100 && dRJetLep1>0.8 && dRJetLep2>0.8 && dRLepLep>0.02 && dPhiMetLep1<1 && dPhiMetLep2<1 && dilep.M()<70 && nbtagsM==0
	&& SelectedEle[0]->pt()>10 && SelectedMuo[0]->pt()>10 && SelectedJet[jetInd]->pt()>400){
       EleMuoDRSelected->Fill(dRLepLep);
       electronPlot(met, SelectedEle, rho);
       muonPlot(met, SelectedMuo, primaryVertex);
       jetPlot(met, SelectedJet, SelectedPrunedMass, SelectedMuo, SelectedEle, jetInd);
       massTauPlot(met, SelectedEle, SelectedMuo, MassSVFit, MassCA);
       massXPlot(met, SelectedEle, SelectedMuo, PrunedJet_prov, XmassCA, XMassSVFit);
       otherPlot(met, Njet, Nele, Nmuo, dRJetZ, nbtagsL, nbtagsM, nbtagsT, true);
     }
     //CONTROL PLOTS
     metPtSelectedNO->Fill(met->begin()->pt());
     if(met->begin()->pt()>100){
       electronDRMetEMUNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[0]->p4()));
       electronDPhiMetEMUNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
       muonDRMetEMUNO->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[0]->p4()));
       muonDPhiMetEMUNO->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
     }
   }
   //CONTROL PLOTS
   if(SelectedEle.size()==1 && SelectedMuo.size()==1){
     pat::JetCollection::const_iterator SelectedJet1;
     pat::JetCollection::const_iterator SelectedJet2;
     pat::JetCollection::const_iterator SelectedJet3;
     float ptMax1=-999; float ptMax2=-999; float ptMax3=-999;
     float massSelected = 0.;
     for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
       float dRmin = 9999.; float mass = 0.;
       for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
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
       if(abs(jet->eta())>2.4) continue;
       if(jet->pt()>400 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax1){
	   ptMax1=jet->pt();
	   SelectedJet1=jet;
	 }
       }
       if((jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75 && (mass>20 && mass<70)){
	 if(jet->pt()>ptMax2){
	   ptMax2=jet->pt();
	   SelectedJet2=jet;
	 }
       }
       if(jet->pt()>400 && (jet->userFloat("tau2")/jet->userFloat("tau1")>0.75)<0.75){
	 if(jet->pt()>ptMax3){
	   ptMax3=jet->pt();
	   massSelected=mass;
	   SelectedJet3=jet;
	 }
       }
     }
     if(ptMax1>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(), SelectedJet1->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedEle[0]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet1->p4()))>2 ) jetSubjettinessNO->Fill(SelectedJet1->userFloat("tau2")/SelectedJet1->userFloat("tau1")>0.75);
     if(ptMax2>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet2->p4())>0.8 
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(), SelectedJet2->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedEle[0]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet2->p4()))>2 ) jetPtNO->Fill(SelectedJet2->pt());
     if(ptMax3>0 && met->begin()->pt()>100
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(), SelectedJet3->p4())>0.8
	&& ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedEle[0]->p4())>0.05
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()))<1 
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()))<1
	&& fabs(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet3->p4()))>2 ) jetMassNO->Fill(massSelected);
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
SidebandFullyLeptonicAnalyzer::beginJob()
{
  Service<TFileService> fs;
  
  metPt             = fs->make<TH1D>("metPt",            "metPt",            5000, 0, 5000);
  NJet              = fs->make<TH1D>("NJet",             "NJet",             5, -0.5, 4.5 );
  NEle              = fs->make<TH1D>("NEle",             "NEle",             5, -0.5, 4.5 );
  NMuo              = fs->make<TH1D>("NMuo",             "NMuo",             5, -0.5, 4.5 );
  NLep              = fs->make<TH1D>("NLep",             "NLep",             5, -0.5, 4.5 );
  metPtSelected     = fs->make<TH1D>("metPtSelected",    "metPtSelected",    1000, 0, 1000);
  NJetSelected      = fs->make<TH1D>("NJetSelected",     "NJetSelected",     5, -0.5, 4.5 );
  NEleSelected      = fs->make<TH1D>("NEleSelected",     "NEleSelected",     5, -0.5, 4.5 );
  NMuoSelected      = fs->make<TH1D>("NMuoSelected",     "NMuoSelected",     5, -0.5, 4.5 );
  NLepSelected      = fs->make<TH1D>("NLepSelected",     "NLepSelected",     5, -0.5, 4.5 );
  NbtagsLSelected   = fs->make<TH1D>("NbtagsLSelected",  "NbtagsLSelected", 15, -0.5, 14.5);
  NbtagsMSelected   = fs->make<TH1D>("NbtagsMSelected",  "NbtagsMSelected", 15, -0.5, 14.5);
  NbtagsTSelected   = fs->make<TH1D>("NbtagsTSelected",  "NbtagsTSelected", 15, -0.5, 14.5);
  JetZdRSelected    = fs->make<TH1D>("JetZdRSelected",   "JetZdRSelected",   100,  0, 5   );
  EleMuoDRSelected  = fs->make<TH1D>("EleMuoDRSelected", "EleMuoDRSelected", 500,  0, 5   );
  tauDecay          = fs->make<TH1D>("tauDecay",         "tauDecay",         10,   0, 10  );
  tauDecaySelected  = fs->make<TH1D>("tauDecaySelected", "tauDecaySelected", 10,   0, 10  );

  MassVisEleEle     = fs->make<TH1D>("MassVisEleEle",    "MassVisEleEle",    500,  0, 500 );
  MassVisMuoMuo     = fs->make<TH1D>("MassVisMuoMuo",    "MassVisMuoMuo",    500,  0, 500 );
  MassVisEleMuo     = fs->make<TH1D>("MassVisEleMuo",    "MassVisEleMuo",    500,  0, 500 );
  MassEffEleEle     = fs->make<TH1D>("MassEffEleEle",    "MassEffEleEle",    500,  0, 500 );
  MassEffMuoMuo     = fs->make<TH1D>("MassEffMuoMuo",    "MassEffMuoMuo",    500,  0, 500 );
  MassEffEleMuo     = fs->make<TH1D>("MassEffEleMuo",    "MassEffEleMuo",    500,  0, 500 );
  MassSvfitEleEle   = fs->make<TH1D>("MassSvfitEleEle",  "MassSvfitEleEle",  500,  0, 500 );
  MassSvfitMuoMuo   = fs->make<TH1D>("MassSvfitMuoMuo",  "MassSvfitMuoMuo",  500,  0, 500 );
  MassSvfitEleMuo   = fs->make<TH1D>("MassSvfitEleMuo",  "MassSvfitEleMuo",  500,  0, 500 );
  MassCAEleEle      = fs->make<TH1D>("MassCAEleEle",     "MassCAEleEle",     500,  0, 500 );
  MassCAMuoMuo      = fs->make<TH1D>("MassCAMuoMuo",     "MassCAMuoMuo",     500,  0, 500 );
  MassCAEleMuo      = fs->make<TH1D>("MassCAEleMuo",     "MassCAEleMuo",     500,  0, 500 );
  XMassVisEleEle    = fs->make<TH1D>("XMassVisEleEle",   "XMassVisEleEle",   3000, 0, 3000);
  XMassVisMuoMuo    = fs->make<TH1D>("XMassVisMuoMuo",   "XMassVisMuoMuo",   3000, 0, 3000);
  XMassVisEleMuo    = fs->make<TH1D>("XMassVisEleMuo",   "XMassVisEleMuo",   3000, 0, 3000);
  XMassEffEleEle    = fs->make<TH1D>("XMassEffEleEle",   "XMassEffEleEle",   3000, 0, 3000);
  XMassEffMuoMuo    = fs->make<TH1D>("XMassEffMuoMuo",   "XMassEffMuoMuo",   3000, 0, 3000);
  XMassEffEleMuo    = fs->make<TH1D>("XMassEffEleMuo",   "XMassEffEleMuo",   3000, 0, 3000);
  XMassSVFitEleEle  = fs->make<TH1D>("XMassSVFitEleEle", "XMassSVFitEleEle", 3000, 0, 3000);
  XMassSVFitMuoMuo  = fs->make<TH1D>("XMassSVFitMuoMuo", "XMassSVFitMuoMuo", 3000, 0, 3000);
  XMassSVFitEleMuo  = fs->make<TH1D>("XMassSVFitEleMuo", "XMassSVFitEleMuo", 3000, 0, 3000);
  XMassCAEleEle     = fs->make<TH1D>("XMassCAEleEle",    "XMassCAEleEle",    3000, 0, 3000);
  XMassCAMuoMuo     = fs->make<TH1D>("XMassCAMuoMuo",    "XMassCAMuoMuo",    3000, 0, 3000);
  XMassCAEleMuo     = fs->make<TH1D>("XMassCAEleMuo",    "XMassCAEleMuo",    3000, 0, 3000);

  jetPtSelected           = fs->make<TH1D>("jetPtSelected",           "jetPtSelected",           1000,    0, 1000 );
  jetEtaSelected          = fs->make<TH1D>("jetEtaSelected",          "jetEtaSelected",          600,    -3, 3    );
  jetMassSelected         = fs->make<TH1D>("jetMassSelected",         "jetMassSelected",         300,     0, 300  );
  jetSubjettinessSelected = fs->make<TH1D>("jetSubjettinessSelected", "jetSubjettinessSelected", 100,     0, 1    );
  jetDRLepSelected        = fs->make<TH1D>("jetDRLepSelected",        "jetDRLepSelected",        100,     0, 5    );
  jetDRMetSelected        = fs->make<TH1D>("jetDRMetSelected",        "jetDRMetSelected",        100,     0, 5    );
  jetDRVisZSelected       = fs->make<TH1D>("jetDRVisZSelected",       "jetDRVisZSelected",       100,     0, 5    );
  jetDPhiMetSelected      = fs->make<TH1D>("jetDPhiMetSelected",      "jetDPhiMetSelected",      100,    -5, 5    );

  muonDR           = fs->make<TH1D>("muonDR",     "muonDR",     500,      0, 5    );
  muonPtLead          = fs->make<TH1D>("muonPtLead",          "muonPtLead",          1000,    0, 1000 );
  muonEtaLead         = fs->make<TH1D>("muonEtaLead",         "muonEtaLead",         600,    -3, 3    );
  muonPFIsoLead       = fs->make<TH1D>("muonPFIsoLead",       "muonPFIsoLead",       1500,    0, 1.5  );
  muonDRMetLead       = fs->make<TH1D>("muonDRMetLead",       "muonDRMetLead",       100,     0, 5    );
  muonDPhiMetLead     = fs->make<TH1D>("muonDPhiMetLead",     "muonDPhiMetLead",     100,    -5, 5    );
  muonDetIsoLead      = fs->make<TH1D>("muonDetIsoLead",      "muonDetIsoLead",     1500,     0, 1.5  );
  muonPtSublead          = fs->make<TH1D>("muonPtSublead",          "muonPtSublead",          1000,    0, 1000 );
  muonEtaSublead         = fs->make<TH1D>("muonEtaSublead",         "muonEtaSublead",         600,    -3, 3    );
  muonPFIsoSublead       = fs->make<TH1D>("muonPFIsoSublead",       "muonPFIsoSublead",       1500,    0, 1.5  );
  muonDRMetSublead       = fs->make<TH1D>("muonDRMetSublead",       "muonDRMetSublead",       100,     0, 5    );
  muonDPhiMetSublead     = fs->make<TH1D>("muonDPhiMetSublead",     "muonDPhiMetSublead",     100,    -5, 5    );
  muonDetIsoSublead      = fs->make<TH1D>("muonDetIsoSublead",      "muonDetIsoSublead",     1500,     0, 1.5  );
  muonPtEMU          = fs->make<TH1D>("muonPtEMU",          "muonPtEMU",          1000,    0, 1000 );
  muonEtaEMU         = fs->make<TH1D>("muonEtaEMU",         "muonEtaEMU",         600,    -3, 3    );
  muonPFIsoEMU       = fs->make<TH1D>("muonPFIsoEMU",       "muonPFIsoEMU",       1500,    0, 1.5  );
  muonDRMetEMU       = fs->make<TH1D>("muonDRMetEMU",       "muonDRMetEMU",       100,     0, 5    );
  muonDPhiMetEMU     = fs->make<TH1D>("muonDPhiMetEMU",     "muonDPhiMetEMU",     100,    -5, 5    );
  muonDetIsoEMU      = fs->make<TH1D>("muonDetIsoEMU",      "muonDetIsoEMU",     1500,     0, 1.5  );

  electronDR          = fs->make<TH1D>("electronDR",          "electronDR",          500,     0, 5    );
  electronPtLead      = fs->make<TH1D>("electronPtLead",      "electronPtLead",      1000,    0, 1000 );
  electronEtaLead     = fs->make<TH1D>("electronEtaLead",     "electronEtaLead",     600,    -3, 3    );
  electronPFIsoLead   = fs->make<TH1D>("electronPFIsoLead",   "electronPFIsoLead",   1500,    0, 1.5  );
  electronDRMetLead   = fs->make<TH1D>("electronDRMetLead",   "electronDRMetLead",   100,     0, 5    );
  electronDPhiMetLead = fs->make<TH1D>("electronDPhiMetLead", "electronDPhiMetLead", 100,    -5, 5    );
  electronDetIsoLead  = fs->make<TH1D>("electronDetIsoLead",  "electronDetIsoLead",  1500,    0, 1.5  );
  electronPtSublead      = fs->make<TH1D>("electronPtSublead",      "electronPtSublead",      1000,    0, 1000 );
  electronEtaSublead     = fs->make<TH1D>("electronEtaSublead",     "electronEtaSublead",     600,    -3, 3    );
  electronPFIsoSublead   = fs->make<TH1D>("electronPFIsoSublead",   "electronPFIsoSublead",   1500,    0, 1.5  );
  electronDRMetSublead   = fs->make<TH1D>("electronDRMetSublead",   "electronDRMetSublead",   100,     0, 5    );
  electronDPhiMetSublead = fs->make<TH1D>("electronDPhiMetSublead", "electronDPhiMetSublead", 100,    -5, 5    );
  electronDetIsoSublead  = fs->make<TH1D>("electronDetIsoSublead",  "electronDetIsoSublead",  1500,    0, 1.5  );  
  electronPtEMU      = fs->make<TH1D>("electronPtEMU",      "electronPtEMU",      1000,    0, 1000 );
  electronEtaEMU     = fs->make<TH1D>("electronEtaEMU",     "electronEtaEMU",     600,    -3, 3    );
  electronPFIsoEMU   = fs->make<TH1D>("electronPFIsoEMU",   "electronPFIsoEMU",   1500,    0, 1.5  );
  electronDRMetEMU   = fs->make<TH1D>("electronDRMetEMU",   "electronDRMetEMU",   100,     0, 5    );
  electronDPhiMetEMU = fs->make<TH1D>("electronDPhiMetEMU", "electronDPhiMetEMU", 100,    -5, 5    );
  electronDetIsoEMU  = fs->make<TH1D>("electronDetIsoEMU",  "electronDetIsoEMU",  1500,    0, 1.5  ); 

  metPtSelectedNO           = fs->make<TH1D>("metPtSelectedNO",           "metPtSelectedNO",        1000, 0, 1000);
  muonDRMetLeadNO           = fs->make<TH1D>("muonDRMetLeadNO",           "muonDRMetLeadNO",           100,  0, 5);
  muonDPhiMetLeadNO         = fs->make<TH1D>("muonDPhiMetLeadNO",         "muonDPhiMetLeadNO",         100, -5, 5);
  muonDRMetSubleadNO        = fs->make<TH1D>("muonDRMetSubleadNO",        "muonDRMetSubleadNO",        100,  0, 5);
  muonDPhiMetSubleadNO      = fs->make<TH1D>("muonDPhiMetSubleadNO",      "muonDPhiMetSubleadNO",      100, -5, 5);
  electronDRMetLeadNO       = fs->make<TH1D>("electronDRMetLeadNO",       "electronDRMetLeadNO",       100,  0, 5);
  electronDPhiMetLeadNO     = fs->make<TH1D>("electronDPhiMetLeadNO",     "electronDPhiMetLeadNO",     100, -5, 5);
  electronDRMetSubleadNO    = fs->make<TH1D>("electronDRMetSubleadNO",    "electronDRMetSubleadNO",    100,  0, 5);
  electronDPhiMetSubleadNO  = fs->make<TH1D>("electronDPhiMetSubleadNO",  "electronDPhiMetSubleadNO",  100, -5, 5);
  muonDRMetEMUNO            = fs->make<TH1D>("muonDRMetEMUNO",            "muonDRMetEMUNO",            100,  0, 5);
  muonDPhiMetEMUNO          = fs->make<TH1D>("muonDPhiMetEMUNO",          "muonDPhiMetEMUNO",          100, -5, 5);
  electronDRMetEMUNO        = fs->make<TH1D>("electronDRMetEMUNO",        "electronDRMetEMUNO",        100,  0, 5);
  electronDPhiMetEMUNO      = fs->make<TH1D>("electronDPhiMetEMUNO",      "electronDPhiMetEMUNO",      100, -5, 5);
  jetPtNO           = fs->make<TH1D>("jetPtNO",           "jetPtNO",           1000,    0, 1000 );
  jetMassNO         = fs->make<TH1D>("jetMassNO",         "jetMassNO",         300,     0, 300  );
  jetSubjettinessNO = fs->make<TH1D>("jetSubjettinessNO", "jetSubjettinessNO", 100,     0, 1    );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SidebandFullyLeptonicAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SidebandFullyLeptonicAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SidebandFullyLeptonicAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SidebandFullyLeptonicAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SidebandFullyLeptonicAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SidebandFullyLeptonicAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void SidebandFullyLeptonicAnalyzer::SelectJet(edm::Handle<pat::JetCollection> CA8JetswithQjets,
			   edm::Handle<pat::JetCollection> CA8JetsPruned,
			   vector<pat::JetCollection::const_iterator> & SelectedJet,
			   vector<float> & SelectedPrunedMass, int & Njet){

  SelectedJet.clear(); SelectedPrunedMass.clear();
  for(pat::JetCollection::const_iterator jet = CA8JetswithQjets->begin(); jet != CA8JetswithQjets->end(); ++jet) {
    float dRmin = 9999.; float mass = 0.;
    for(pat::JetCollection::const_iterator jetPruned = CA8JetsPruned->begin(); jetPruned != CA8JetsPruned->end(); ++jetPruned) {
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
    if(jet->userFloat("tau2")/jet->userFloat("tau1")>0.75) continue;
    Njet = Njet + 1;
    SelectedJet.push_back(jet);
    SelectedPrunedMass.push_back(mass);
   }
}


void SidebandFullyLeptonicAnalyzer::SelectElectron(edm::Handle<pat::ElectronCollection> eleH, vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho){
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<40) continue;
    if(!(abs(electron->eta())<1.4442 || (abs(electron->eta())>1.5666 && abs(electron->eta())<2.5))) continue;
    if(abs(electron->phi())>3.2) continue;
    if(electron->userInt("HEEPId")!=0) continue;
    Nele = Nele + 1;
    SelectedEle.push_back(electron);
  }
}


void SidebandFullyLeptonicAnalyzer::SelectElectronCutBased(edm::Handle<pat::ElectronCollection> eleH, vector<pat::ElectronCollection::const_iterator> & SelectedEle, int & Nele, float rho, reco::Vertex primaryVertex){
  for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
    if(electron->pt()<10) continue;
    if(fabs(electron->superCluster()->eta())<=1.479){
      if(electron->deltaEtaSuperClusterTrackAtVtx()>=0.004) continue;
      if(electron->deltaPhiSuperClusterTrackAtVtx()>=0.030) continue;
      if(electron->sigmaIetaIeta()>=0.01) continue;
      if(electron->hadronicOverEm()>=0.12) continue;
      if(fabs(electron->gsfTrack()->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(electron->gsfTrack()->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(electron->passConversionVeto()==0) continue;
      if(electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()!=0) continue;
    }
    if(fabs(electron->superCluster()->eta())>1.479 && fabs(electron->superCluster()->eta())<2.5){
      if(electron->deltaEtaSuperClusterTrackAtVtx()>=0.005) continue;
      if(electron->deltaPhiSuperClusterTrackAtVtx()>=0.020) continue;
      if(electron->sigmaIetaIeta()>=0.03) continue;
      if(electron->hadronicOverEm()>=0.10) continue;
      if(fabs(electron->gsfTrack()->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(electron->gsfTrack()->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(electron->passConversionVeto()==0) continue;
      if(electron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()!=0) continue;
    }
    Nele = Nele + 1;
    SelectedEle.push_back(electron);
  }
}


void SidebandFullyLeptonicAnalyzer::SelectTightMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(muon->pt()<10) continue;
    if(abs(muon->eta())>2.4) continue;
    if(abs(muon->phi())>3.2) continue;
    if(!(muon->isGlobalMuon())) continue;
    if(!(muon->isPFMuon())) continue;
    if(muon->globalTrack()->normalizedChi2()>=10) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(muon->dB())>=0.2 ) continue;
    if(fabs(muon->muonBestTrack()->dz(primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    Nmuo = Nmuo + 1;
    SelectedMuo.push_back(muon);
  }
}



void SidebandFullyLeptonicAnalyzer::SelectHighptMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo, int & Nmuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(abs(cktTrack->eta())>2.4) continue;
    if(abs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    SelectedMuo.push_back(muon);
    Nmuo = Nmuo + 1;
  }
}


void SidebandFullyLeptonicAnalyzer::SelectTrackerMuon(edm::Handle<pat::MuonCollection> muoH, vector<pat::MuonCollection::const_iterator> & SelectedMuo,int & Nmuo, reco::Vertex primaryVertex){
  for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
    if(muon->pt()<10) continue;
    if(abs(muon->eta())>2.4) continue;
    if(abs(muon->phi())>3.2) continue;
    if(!(muon->isTrackerMuon())) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(muon->muonBestTrack()->dz(primaryVertex.position()))>=0.5) continue;
    if(fabs(muon->dB())>=0.2 ) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    if((muon->muonBestTrack()->ptError()/muon->muonBestTrack()->pt())>0.3) continue;
    SelectedMuo.push_back(muon);
    Nmuo = Nmuo + 1;
  }
}

void SidebandFullyLeptonicAnalyzer::SelectTrackerGlobalID(vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex, bool & hasAtLeastOneHighPtMuo){
  for(unsigned int i=0; i<SelectedMuo.size(); i++){
    if(!(SelectedMuo[i]->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = (muon::tevOptimized(*(SelectedMuo[i]), 200, 17., 40., 0.25)).first;
    if(cktTrack->pt()<10) continue;
    if(abs(cktTrack->eta())>2.4) continue;
    if(abs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(SelectedMuo[i]->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(SelectedMuo[i]->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(SelectedMuo[i]->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(SelectedMuo[i]->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    hasAtLeastOneHighPtMuo = true;
  }
}


float SidebandFullyLeptonicAnalyzer::ElectronPFIso(pat::ElectronCollection::const_iterator electron, float rho){
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

float SidebandFullyLeptonicAnalyzer::MuonPFIso(pat::MuonCollection::const_iterator muon, bool highpt){
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

float SidebandFullyLeptonicAnalyzer::MuonDETIso(pat::MuonCollection::const_iterator SelectedMuo, vector<pat::MuonCollection::const_iterator> SelectedHighptMuo){
  float isovar = SelectedMuo->trackIso()/SelectedMuo->pt();
  for(unsigned int j = 0; j< SelectedHighptMuo.size();++j){
    if(SelectedMuo==SelectedHighptMuo[j]) continue;
    double dR = ROOT::Math::VectorUtil::DeltaR(SelectedMuo->p4(),SelectedHighptMuo[j]->p4());
    if(dR < 0.3) isovar = isovar - ((SelectedHighptMuo[j]->track()->pt())/SelectedMuo->pt());
  }
  return isovar;
}

bool SidebandFullyLeptonicAnalyzer::ElectronDETIso(pat::ElectronCollection::const_iterator electron, float rho){
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

void SidebandFullyLeptonicAnalyzer::BtagVeto(edm::Handle<pat::JetCollection> ak5jetCands, int & nbtagsL, int & nbtagsM, int & nbtagsT, 
				     pat::JetCollection::const_iterator SelectedJet){
  for(pat::JetCollection::const_iterator ak5 = ak5jetCands->begin(); ak5 != ak5jetCands->end(); ++ak5) {
    if(ROOT::Math::VectorUtil::DeltaR(ak5->p4(),SelectedJet->p4())<0.8) continue;
    double discCSV = ak5->bDiscriminator("combinedSecondaryVertexBJetTags");
    if(discCSV>0.244) nbtagsL++;  //loose working point
    if(discCSV>0.679) nbtagsM++;  //medium working point
    if(discCSV>0.898) nbtagsT++;  //tight working point
  }
}

float SidebandFullyLeptonicAnalyzer::SVFitMass(edm::Handle<pat::METCollection> met, edm::Handle<pat::METCollection> metRaw, TLorentzVector & SVFitTauTau,
			    LorentzVector lep1, LorentzVector lep2){
  TMatrixD covMET(2, 2); // PFMET significance matrix
  covMET[0][0] = (metRaw->front() ).getSignificanceMatrix()(0,0);
  covMET[1][0] = (metRaw->front() ).getSignificanceMatrix()(1,0);
  covMET[0][1] = (metRaw->front() ).getSignificanceMatrix()(0,1);
  covMET[1][1] = (metRaw->front() ).getSignificanceMatrix()(1,1);
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep2));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->begin()->momentum(), covMET, 0);
  algo.addLogM(false);
  algo.integrateMarkovChain();
  SVFitTauTau.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
  double MassSVFit = algo.getMass(); 
  return MassSVFit;
}


void SidebandFullyLeptonicAnalyzer::CollinearApproximation(edm::Handle<pat::METCollection> met, TLorentzVector & CATauTau, TLorentzVector PrunedJet, LorentzVector lep1, LorentzVector lep2, bool & CA){
  float a = (lep2.py()*met->begin()->px()-lep2.px()*met->begin()->py())/
    (lep1.px()*lep2.py()-lep1.py()*lep2.px());
  float b = (lep1.py()*met->begin()->px()-lep1.px()*met->begin()->py())/
    (lep2.px()*lep1.py()-lep2.py()*lep1.px());
  if(((1+a)*(1+b))>0) {
    CATauTau.SetPxPyPzE((1+a)*lep1.px()+(1+b)*lep2.px(),
			 (1+a)*lep1.py()+(1+b)*lep2.py(),
			 (1+a)*lep1.pz()+(1+b)*lep2.pz(),
			 (1+a)*lep1.energy()+(1+b)*lep2.energy());
    CA=true;
  }
  else CA = false;
}


void SidebandFullyLeptonicAnalyzer::jetPlot(edm::Handle<pat::METCollection> met, vector<pat::JetCollection::const_iterator> SelectedJet, vector<float> SelectedPrunedMass, vector<pat::MuonCollection::const_iterator> SelectedMuo, vector<pat::ElectronCollection::const_iterator> SelectedEle, int jetInd){
  jetPtSelected->Fill(SelectedJet[jetInd]->pt());
  jetEtaSelected->Fill(SelectedJet[jetInd]->eta());
  jetMassSelected->Fill(SelectedPrunedMass[jetInd]);
  jetSubjettinessSelected->Fill(SelectedJet[jetInd]->userFloat("tau2")/SelectedJet[jetInd]->userFloat("tau1"));
  if(SelectedEle.size()==1 || SelectedEle.size()==2) jetDRLepSelected->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedJet[jetInd]->p4()));
  if(SelectedEle.size()==2) jetDRLepSelected->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedEle[1]->p4(),SelectedJet[jetInd]->p4()));
  if(SelectedMuo.size()==1 || SelectedMuo.size()==2) jetDRLepSelected->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedJet[jetInd]->p4()));
  if(SelectedMuo.size()==2) jetDRLepSelected->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuo[1]->p4(),SelectedJet[jetInd]->p4()));
  jetDRMetSelected->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedJet[jetInd]->p4()));
  jetDPhiMetSelected->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedJet[jetInd]->p4()));
  TLorentzVector VisZ;
  if(SelectedEle.size()==2) {VisZ.SetPxPyPzE(SelectedEle[0]->px()+SelectedEle[1]->px(),SelectedEle[0]->py()+SelectedEle[1]->py(),
					     SelectedEle[0]->pz()+SelectedEle[1]->pz(),SelectedEle[0]->energy()+SelectedEle[1]->energy());}
  if(SelectedMuo.size()==2) {VisZ.SetPxPyPzE(SelectedMuo[0]->px()+SelectedMuo[1]->px(),SelectedMuo[0]->py()+SelectedMuo[1]->py(),
					     SelectedMuo[0]->pz()+SelectedMuo[1]->pz(),SelectedMuo[0]->energy()+SelectedMuo[1]->energy());}
  if(SelectedEle.size()==1) {VisZ.SetPxPyPzE(SelectedEle[0]->px()+SelectedMuo[0]->px(),SelectedEle[0]->py()+SelectedMuo[0]->py(),
					     SelectedEle[0]->pz()+SelectedMuo[0]->pz(),SelectedEle[0]->energy()+SelectedMuo[0]->energy());}
  jetDRVisZSelected->Fill(ROOT::Math::VectorUtil::DeltaR(VisZ,SelectedJet[jetInd]->p4()));
}

void SidebandFullyLeptonicAnalyzer::muonPlot(edm::Handle<pat::METCollection> met, vector<pat::MuonCollection::const_iterator> SelectedMuo, reco::Vertex primaryVertex){
  if(SelectedMuo.size()==2){
    muonDR->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedMuo[0]->p4(),SelectedMuo[1]->p4()));
    if(SelectedMuo[0]->pt()>SelectedMuo[1]->pt()){
      muonDetIsoLead->Fill(MuonDETIso(SelectedMuo[0], SelectedMuo));
      muonPtLead->Fill(SelectedMuo[0]->pt());
      muonEtaLead->Fill(SelectedMuo[0]->eta());
      muonPFIsoLead->Fill(MuonPFIso(SelectedMuo[0],false));
      muonDPhiMetLead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
      muonDRMetLead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[0]->p4()));
      muonDetIsoSublead->Fill(MuonDETIso(SelectedMuo[1], SelectedMuo));
      muonPtSublead->Fill(SelectedMuo[1]->pt());
      muonEtaSublead->Fill(SelectedMuo[1]->eta());
      muonPFIsoSublead->Fill(MuonPFIso(SelectedMuo[1],false));
      muonDPhiMetSublead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()));
      muonDRMetSublead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[1]->p4()));
    } else {
      muonDetIsoLead->Fill(MuonDETIso(SelectedMuo[1], SelectedMuo));
      muonPtLead->Fill(SelectedMuo[1]->pt());
      muonEtaLead->Fill(SelectedMuo[1]->eta());
      muonPFIsoLead->Fill(MuonPFIso(SelectedMuo[1],false));
      muonDPhiMetLead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[1]->p4()));
      muonDRMetLead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[1]->p4()));
      muonDetIsoSublead->Fill(MuonDETIso(SelectedMuo[0], SelectedMuo));
      muonPtSublead->Fill(SelectedMuo[0]->pt());
      muonEtaSublead->Fill(SelectedMuo[0]->eta());
      muonPFIsoSublead->Fill(MuonPFIso(SelectedMuo[0],false));
      muonDPhiMetSublead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
      muonDRMetSublead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[0]->p4()));
    } 
  } else {
      muonDetIsoEMU->Fill(MuonDETIso(SelectedMuo[0], SelectedMuo));
      muonPtEMU->Fill(SelectedMuo[0]->pt());
      muonEtaEMU->Fill(SelectedMuo[0]->eta());
      muonPFIsoEMU->Fill(MuonPFIso(SelectedMuo[0],true));
      muonDPhiMetEMU->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedMuo[0]->p4()));
      muonDRMetEMU->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedMuo[0]->p4()));
  }
}


void SidebandFullyLeptonicAnalyzer::electronPlot(edm::Handle<pat::METCollection> met, vector<pat::ElectronCollection::const_iterator> SelectedEle, float rho){
  if(SelectedEle.size()==2){
    electronDR->Fill(ROOT::Math::VectorUtil::DeltaR(SelectedEle[0]->p4(),SelectedEle[1]->p4()));
    if(SelectedEle[0]->pt()>SelectedEle[1]->pt()){
      electronPtLead->Fill(SelectedEle[0]->pt());
      electronEtaLead->Fill(SelectedEle[0]->eta());
      electronPFIsoLead->Fill(ElectronPFIso(SelectedEle[0],rho));
      electronDetIsoLead->Fill((SelectedEle[0]->userIso(0)+SelectedEle[0]->userIso(1)+SelectedEle[0]->userIso(2))/SelectedEle[0]->pt());
      electronDRMetLead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[0]->p4()));
      electronDPhiMetLead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
      electronPtSublead->Fill(SelectedEle[1]->pt());
      electronEtaSublead->Fill(SelectedEle[1]->eta());
      electronPFIsoSublead->Fill(ElectronPFIso(SelectedEle[1],rho));
      electronDetIsoSublead->Fill((SelectedEle[1]->userIso(0)+SelectedEle[1]->userIso(1)+SelectedEle[1]->userIso(2))/SelectedEle[1]->pt());
      electronDRMetSublead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[1]->p4()));
      electronDPhiMetSublead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()));
    } else {
      electronPtLead->Fill(SelectedEle[1]->pt());
      electronEtaLead->Fill(SelectedEle[1]->eta());
      electronPFIsoLead->Fill(ElectronPFIso(SelectedEle[1],rho));
      electronDetIsoLead->Fill((SelectedEle[1]->userIso(0)+SelectedEle[1]->userIso(1)+SelectedEle[1]->userIso(2))/SelectedEle[1]->pt());
      electronDRMetLead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[1]->p4()));
      electronDPhiMetLead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[1]->p4()));
      electronPtSublead->Fill(SelectedEle[0]->pt());
      electronEtaSublead->Fill(SelectedEle[0]->eta());
      electronPFIsoSublead->Fill(ElectronPFIso(SelectedEle[0],rho));
      electronDetIsoSublead->Fill((SelectedEle[0]->userIso(0)+SelectedEle[0]->userIso(1)+SelectedEle[0]->userIso(2))/SelectedEle[0]->pt());
      electronDRMetSublead->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[0]->p4()));
      electronDPhiMetSublead->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
    }
  } else {
    electronPtEMU->Fill(SelectedEle[0]->pt());
    electronEtaEMU->Fill(SelectedEle[0]->eta());
    electronPFIsoEMU->Fill(ElectronPFIso(SelectedEle[0],rho));
    electronDetIsoEMU->Fill((SelectedEle[0]->userIso(0)+SelectedEle[0]->userIso(1)+SelectedEle[0]->userIso(2))/SelectedEle[0]->pt());
    electronDRMetEMU->Fill(ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),SelectedEle[0]->p4()));
    electronDPhiMetEMU->Fill(ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),SelectedEle[0]->p4()));
  }
}

void SidebandFullyLeptonicAnalyzer::massTauPlot(edm::Handle<pat::METCollection> met, vector<pat::ElectronCollection::const_iterator> SelectedEle, vector<pat::MuonCollection::const_iterator> SelectedMuo, float MassSVFit, float MassCA){
  if(SelectedEle.size()==2){
    math::PtEtaPhiELorentzVector diele; 
    math::PtEtaPhiELorentzVector dielemet; 
    diele = SelectedEle[0]->p4()+SelectedEle[1]->p4();
    dielemet = SelectedEle[0]->p4()+SelectedEle[1]->p4()+met->begin()->p4();
    MassVisEleEle->Fill(diele.mass());
    MassEffEleEle->Fill(dielemet.mass());
    MassSvfitEleEle->Fill(MassSVFit);
    MassCAEleEle->Fill(MassCA);
  }
  if(SelectedMuo.size()==2){
    math::PtEtaPhiELorentzVector dimuo; 
    math::PtEtaPhiELorentzVector dimuomet; 
    dimuo = SelectedMuo[0]->p4()+SelectedMuo[1]->p4();
    dimuomet = SelectedMuo[0]->p4()+SelectedMuo[1]->p4()+met->begin()->p4();
    MassVisMuoMuo->Fill(dimuo.mass());
    MassEffMuoMuo->Fill(dimuomet.mass());
    MassSvfitMuoMuo->Fill(MassSVFit);
    MassCAMuoMuo->Fill(MassCA);
  }
  if(SelectedEle.size()==1 && SelectedMuo.size()==1){
    math::PtEtaPhiELorentzVector elemuo; 
    math::PtEtaPhiELorentzVector elemuomet; 
    elemuo = SelectedEle[0]->p4()+SelectedMuo[0]->p4();
    elemuomet = SelectedEle[0]->p4()+SelectedMuo[0]->p4()+met->begin()->p4();
    MassVisEleMuo->Fill(elemuo.mass());
    MassEffEleMuo->Fill(elemuomet.mass());
    MassSvfitEleMuo->Fill(MassSVFit);
    MassCAEleMuo->Fill(MassCA);
  }
}

void SidebandFullyLeptonicAnalyzer::massXPlot(edm::Handle<pat::METCollection> met, vector<pat::ElectronCollection::const_iterator> SelectedEle, vector<pat::MuonCollection::const_iterator> SelectedMuo, math::PtEtaPhiMLorentzVector PrunedJet_prov, float XmassCA, float XmassSVFit){
  if(SelectedEle.size()==2){
    math::PtEtaPhiELorentzVector dielejet; 
    math::PtEtaPhiELorentzVector dielemetjet; 
    dielejet = SelectedEle[0]->p4()+SelectedEle[1]->p4()+PrunedJet_prov;
    dielemetjet = SelectedEle[0]->p4()+SelectedEle[1]->p4()+met->begin()->p4()+PrunedJet_prov;
    XMassVisEleEle->Fill(dielejet.mass());
    XMassEffEleEle->Fill(dielemetjet.mass());
    XMassSVFitEleEle->Fill(XmassSVFit);
    XMassCAEleEle->Fill(XmassCA);
  }
  if(SelectedMuo.size()==2){
    math::PtEtaPhiELorentzVector dimuojet; 
    math::PtEtaPhiELorentzVector dimuometjet; 
    dimuojet = SelectedMuo[0]->p4()+SelectedMuo[1]->p4()+PrunedJet_prov;
    dimuometjet = SelectedMuo[0]->p4()+SelectedMuo[1]->p4()+met->begin()->p4()+PrunedJet_prov;
    XMassVisMuoMuo->Fill(dimuojet.mass());
    XMassEffMuoMuo->Fill(dimuometjet.mass());
    XMassSVFitMuoMuo->Fill(XmassSVFit);
    XMassCAMuoMuo->Fill(XmassCA);
  }
  if(SelectedEle.size()==1 && SelectedMuo.size()==1){
    math::PtEtaPhiELorentzVector elemuojet; 
    math::PtEtaPhiELorentzVector elemuojetmet; 
    elemuojet = SelectedEle[0]->p4()+SelectedMuo[0]->p4()+PrunedJet_prov;
    elemuojetmet = SelectedEle[0]->p4()+SelectedMuo[0]->p4()+met->begin()->p4()+PrunedJet_prov;
    XMassVisEleMuo->Fill(elemuojet.mass());
    XMassEffEleMuo->Fill(elemuojetmet.mass());
    XMassSVFitEleMuo->Fill(XmassSVFit);
    XMassCAEleMuo->Fill(XmassCA);
  }
}


void SidebandFullyLeptonicAnalyzer::otherPlot(edm::Handle<pat::METCollection> met, int Njet, int Nele, int Nmuo, float dRJetZ, int nbtagsL, int nbtagsM, int nbtagsT, bool selected){
  if(selected){
    metPtSelected->Fill(met->begin()->pt());
    NJetSelected->Fill(Njet);
    NEleSelected->Fill(Nele);
    NMuoSelected->Fill(Nmuo);
    NLepSelected->Fill(Nele+Nmuo);
    NbtagsLSelected->Fill(nbtagsL);
    NbtagsMSelected->Fill(nbtagsM);
    NbtagsTSelected->Fill(nbtagsT);
    JetZdRSelected->Fill(dRJetZ);
  } else {
   metPt->Fill(met->begin()->pt());
   NJet->Fill(Njet);
   NEle->Fill(Nele);
   NMuo->Fill(Nmuo);
   NLep->Fill(Nele+Nmuo);
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(SidebandFullyLeptonicAnalyzer);
