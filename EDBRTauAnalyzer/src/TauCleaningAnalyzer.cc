// -*- C++ -*-
//
// Package:    TauCleaningAnalyzer
// Class:      TauCleaningAnalyzer
// 
/**\class TauCleaningAnalyzer TauCleaningAnalyzer.cc BulkG_TauTau/TauCleaningAnalyzer/src/TauCleaningAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aniello Spiezia,21 1-007,+41227676459,
//         Created:  Mon Jul 15 17:08:06 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TMath.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include <vector>
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "TLorentzVector.h"



//
// class declaration
//

class TauCleaningAnalyzer : public edm::EDAnalyzer {
public:
  explicit TauCleaningAnalyzer(const edm::ParameterSet&);
  ~TauCleaningAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  
  double deltaR_;
  double massMin_; double massMax_;
  bool signal_;

  TH1D* Ntau;  TH1D* Nele;  TH1D* Nmuo;  TH1D* tauDecay;

  TH1D* VertexFake;
  TH1D* VertexRho;
  TH1D* VertexNdof;
  TH1D* VertexZ;
  
  TH1D* X_pt; TH1D* X_mass; TH1D* X_phi; TH1D* X_eta;
  TH1D* Ztau_mass; TH1D* Ztau_pt; TH1D* Ztau_phi; TH1D* Ztau_eta;
  
  TH1D* genTauPt1; TH1D* genTauPt2;
  TH1D* ele_pt; TH1D* muo_pt;
  TH1D* genTauPt_mutau; TH1D* genTauPtWithNeutrino_mutau;
  TH1D* genMuoPt_mutau;
  TH1D* genTauPhi1; TH1D* genTau1_eta; TH1D* genTau2_phi; TH1D* genTau2_eta;
  TH1D* genTau_DeltaR;
  TH1D* genTau_DeltaRWithNeutrino; TH1D* genTau_DeltaPhi; TH1D* genTau_DeltaEta;
  TH1D* deltaR1GenReco; TH1D* deltaR2GenReco;
  TH1D* deltaRHad;

  TH1D* h_mutau_decayModeFinding;  TH1D* h_mutau_againstMuonMedium;        TH1D* h_mutau_againstElectronMedium;     TH1D* h_mutau_ByVLooseCombinedIsolationDBSumPtCorr; 
  TH1D* h_mutau_againstMuonLoose;  TH1D* h_mutau_againstElectronLoose;     TH1D* h_mutau_againstMuonLoose2;         TH1D* h_mutau_againstMuonMedium2;
  TH1D* h_mutau_againstMuonLoose3; TH1D* h_mutau_againstMuonTight3; TH1D* h_mutau_againstElectronLooseMVA3; TH1D* h_mutau_againstElectronMediumMVA3; 
  TH1D* h_mutau_byCombinedIsolationDeltaBetaCorrRaw3Hits; TH1D* h_mutau_byLooseCombinedIsolationDeltaBetaCorr3Hits; TH1D* h_mutau_byMediumCombinedIsolationDeltaBetaCorr3Hits;

  TH1D* deltaRGenReco_mutau; 
  TH1D* recoTauPt_mutau; TH1D* recoTauMass_mutau; TH1D* recoTauDecayMode; TH1D* recoMuonPt;
  TH1D* recoTauleadPFCandPt; TH1D* recoTauPlusMuonPt; TH1D* recoLeadMinusMuonPtNorm; TH1D* recoLeadMinusMuonPt;
  TH1D* ptResol_mutau; TH1D* ptRatio_mutau; TH2D* Response_mutau;
  TH2D* ResponseVSdR_mutau; 
  TProfile* Response_mutau_profile;
  TProfile* Response_mutau_profile_pt;

  TH1D* tauSelection_mutau; 
  TH1D* match_2tauh; TH2D* match_2tauh_VS_dR;

  TH1D* DRGenTauPt; TH1D* DRRecoTauPt; TH1D* DRGenTauEta; TH1D* DRRecoTauEta; TH1D* DRMuonTau;
  TH1D* Nmatch;

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
TauCleaningAnalyzer::TauCleaningAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  deltaR_ = iConfig.getUntrackedParameter<double>("deltaR");
  massMin_ = iConfig.getUntrackedParameter<double>("massMin");
  massMax_ = iConfig.getUntrackedParameter<double>("massMax");
  signal_ = iConfig.getUntrackedParameter<bool>("signal");

}


TauCleaningAnalyzer::~TauCleaningAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauCleaningAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);
  Handle<reco::PFTauCollection> tauHandle;
  iEvent.getByLabel("hpsPFTauProducer",tauHandle);
  Handle<reco::PFTauDiscriminator> decayModeFinding;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",decayModeFinding);
  Handle<reco::PFTauDiscriminator> againstMuonMedium;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumMuonRejection",againstMuonMedium);
  Handle<reco::PFTauDiscriminator> againstElectronMedium;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumElectronRejection",againstElectronMedium);
  Handle<reco::PFTauDiscriminator> ByVLooseCombinedIsolationDBSumPtCorr;
  iEvent.getByLabel("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr",ByVLooseCombinedIsolationDBSumPtCorr);
  Handle<reco::PFTauDiscriminator> againstMuonLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection",againstMuonLoose);
  Handle<reco::PFTauDiscriminator> againstElectronLoose;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseElectronRejection",againstElectronLoose);
  Handle<reco::PFTauDiscriminator> againstElectronLooseMVA3;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA3LooseElectronRejection",againstElectronLooseMVA3);
  Handle<reco::PFTauDiscriminator> againstElectronMediumMVA3;
  iEvent.getByLabel("hpsPFTauDiscriminationByMVA3MediumElectronRejection",againstElectronMediumMVA3);
  Handle<reco::PFTauDiscriminator> againstMuonLoose2;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection2",againstMuonLoose2);
  Handle<reco::PFTauDiscriminator> againstMuonMedium2;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumMuonRejection2",againstMuonMedium2);
  Handle<reco::PFTauDiscriminator> againstMuonLoose3;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseMuonRejection3",againstMuonLoose3);
  Handle<reco::PFTauDiscriminator> againstMuonTight3;
  iEvent.getByLabel("hpsPFTauDiscriminationByTightMuonRejection3",againstMuonTight3);
  Handle<reco::PFTauDiscriminator> byCombinedIsolationDeltaBetaCorrRaw3Hits;
  iEvent.getByLabel("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits",byCombinedIsolationDeltaBetaCorrRaw3Hits);
  Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDeltaBetaCorr3Hits;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits",byLooseCombinedIsolationDeltaBetaCorr3Hits);
  Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDeltaBetaCorr3Hits;
  iEvent.getByLabel("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits",byMediumCombinedIsolationDeltaBetaCorr3Hits);
 						     
  Handle<reco::MuonCollection> muH;
  iEvent.getByLabel("muons", muH);  
  Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel("gsfElectrons", elH);

  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);
  VertexFake->Fill(primaryVertex.isFake());
  VertexRho->Fill(primaryVertex.position().rho());
  VertexNdof->Fill(primaryVertex.ndof());
  VertexZ->Fill(primaryVertex.position().Z());

  bool tauele1 = false; bool taumuo1 = false;
  bool tauele2 = false; bool taumuo2 = false;
  vector<math::PtEtaPhiELorentzVector> gentau;
  vector<math::PtEtaPhiELorentzVector> gentauHad;
  vector<math::PtEtaPhiELorentzVector> gentauWithNeutrino;
  vector<math::PtEtaPhiELorentzVector> gentauHadWithNeutrino;
  vector<math::PtEtaPhiELorentzVector> genele;
  vector<math::PtEtaPhiELorentzVector> genmuo;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if(abs(genPart.pdgId())==2212) continue;
    const reco::Candidate * mom = genPart.mother();
    if(abs(genPart.pdgId())==15 && genPart.status()!=3 && (abs(mom->pdgId())==23 || abs(mom->pdgId())==15)){
      math::PtEtaPhiELorentzVector gentau_prov;
      math::PtEtaPhiELorentzVector gentauHad_prov;
      math::PtEtaPhiELorentzVector gentauWithNeutrino_prov;
      math::PtEtaPhiELorentzVector gentauHadWithNeutrino_prov;
      math::PtEtaPhiELorentzVector gen_prov;
      bool genTauHadBool = false;
      if((!genPart.pt()>10)) continue;
      gentau_prov=genPart.p4();
      gentauWithNeutrino_prov=genPart.p4();
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if(genPart.pdgId()== 15 && abs(daughter->pdgId())==11) tauele1 = true;
	if(genPart.pdgId()== 15 && abs(daughter->pdgId())==13) taumuo1 = true;
	if(genPart.pdgId()==-15 && abs(daughter->pdgId())==11) tauele2 = true;
	if(genPart.pdgId()==-15 && abs(daughter->pdgId())==13) taumuo2 = true;
	if(abs(daughter->pdgId())!=11 && abs(daughter->pdgId())!=12 && abs(daughter->pdgId())!=13 && abs(daughter->pdgId())!=14 
	   && abs(daughter->pdgId())!=15 && abs(daughter->pdgId())!=16){
	  gentauHad_prov=genPart.p4();
	  gentauHadWithNeutrino_prov=genPart.p4();
	  genTauHadBool = true;
	}
	gen_prov=daughter->p4();
	if(abs(daughter->pdgId())==11) genele.push_back(gen_prov);
	if(abs(daughter->pdgId())==13) genmuo.push_back(gen_prov);
      }
      for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	const reco::Candidate * daughter = genPart.daughter(ndaugh);
	if((abs(daughter->pdgId())==12 || abs(daughter->pdgId())==14 ||abs(daughter->pdgId())==16) && daughter->status()==1){
	  gentau_prov    = gentau_prov    - daughter->p4();
	  if(genTauHadBool) gentauHad_prov = gentauHad_prov - daughter->p4();
	}
      }
      if(gentau_prov.pt()>10){
	gentau.push_back(gentau_prov);
	gentauWithNeutrino.push_back(gentauWithNeutrino_prov);
	if(genTauHadBool){
	  gentauHad.push_back(gentauHad_prov);
	  gentauHadWithNeutrino.push_back(gentauWithNeutrino_prov);
	}
      }
    }
  }
	
			
  if(signal_ && gentau.size()==2 && ((taumuo1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && taumuo2))){
    Ntau->Fill(gentau.size());
    Nele->Fill(genele.size());
    Nmuo->Fill(genmuo.size());
    genTau_DeltaR->Fill(ROOT::Math::VectorUtil::DeltaR(gentau[0],gentau[1]));
    genTau_DeltaRWithNeutrino->Fill(ROOT::Math::VectorUtil::DeltaR(gentauWithNeutrino[0],gentauWithNeutrino[1]));
    genTau_DeltaPhi->Fill(ROOT::Math::VectorUtil::DeltaPhi(gentau[0],gentau[1]));
    genTau_DeltaEta->Fill(gentau[0].eta()-gentau[1].eta());
    if(gentau[0].pt()>gentau[1].pt()){
      genTauPt1->Fill(gentau[0].pt());   genTauPt2->Fill(gentau[1].pt());
      genTauPhi1->Fill(gentau[0].phi()); genTau2_phi->Fill(gentau[1].phi());
      genTau1_eta->Fill(gentau[0].eta()); genTau2_eta->Fill(gentau[1].eta());
    } else {
      genTauPt1->Fill(gentau[1].pt());   genTauPt2->Fill(gentau[0].pt());
      genTauPhi1->Fill(gentau[1].phi()); genTau2_phi->Fill(gentau[0].phi());
      genTau1_eta->Fill(gentau[1].eta()); genTau2_eta->Fill(gentau[0].eta());
    }
    for(unsigned int k=0; k<genele.size(); k++){ele_pt->Fill(genele[k].pt());}
    for(unsigned int k=0; k<genmuo.size(); k++){muo_pt->Fill(genmuo[k].pt());}
    Ztau_mass->Fill((gentau[0]+gentau[1]).M());
    Ztau_pt->Fill((gentau[0]+gentau[1]).pt());
    Ztau_phi->Fill((gentau[0]+gentau[1]).phi());
    Ztau_eta->Fill((gentau[0]+gentau[1]).eta());
			    
    math::PtEtaPhiELorentzVector genZ1;
    math::PtEtaPhiELorentzVector genZ2;
    for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
      const reco::GenParticle & genPart = (*genParts)[ngenPart];
      if(genPart.pdgId()==23 && genPart.status()!=3 && genZ1.px()==0) genZ1=genPart.p4();
      if(genPart.pdgId()==23 && genPart.status()!=3 && genZ1.px()!=0) genZ2=genPart.p4();
      if((abs(genPart.pdgId())!=23) || (genPart.status()==3)) continue;
    }
    X_pt->Fill((genZ1+genZ2).pt());
    X_phi->Fill((genZ1+genZ2).phi());
    X_eta->Fill((genZ1+genZ2).eta());
    X_mass->Fill((genZ1+genZ2).M()) ;
  }
  
  if(tauele1 && tauele2)                                                          tauDecay->Fill(1);
  else if(taumuo1 && taumuo2)                                                     tauDecay->Fill(2);
  else if((tauele1 && taumuo2) || (taumuo1 && tauele2))                           tauDecay->Fill(3);
  else if((tauele1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && tauele2)) tauDecay->Fill(4);
  else if((taumuo1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && taumuo2)) tauDecay->Fill(5);
  else if(!taumuo1 && !taumuo2 && !tauele1 && !tauele2)                           tauDecay->Fill(6);
  else                                                                            tauDecay->Fill(7);
			
  //CASE MU-TAU_H
  if(signal_ && gentau.size()==2 && ((taumuo1 && !taumuo2 && !tauele2) || (!taumuo1 && !tauele1 && taumuo2))){

    float deltaR_muon=99.; int muonIndex = -1;
    for ( unsigned int i=0; i<muH->size(); ++i ){
      const reco::Muon& muon(muH->at(i));
      if(ROOT::Math::VectorUtil::DeltaR(genmuo[0],muon.p4())<deltaR_ && ROOT::Math::VectorUtil::DeltaR(genmuo[0],muon.p4())<deltaR_muon){
	deltaR_muon = ROOT::Math::VectorUtil::DeltaR(genmuo[0],muon.p4());
	muonIndex = i;
      }
    }
    if(muonIndex!=-1) {
      const reco::Muon& muon(muH->at(muonIndex));
      recoMuonPt->Fill(muon.pt());
    }

    int Nmatching=0;
    for(unsigned int i=0;i<tauHandle->size();++i) {
      reco::PFTauRef pfTauRef(tauHandle,i);
      for(unsigned int j=0; j<gentauHad.size(); j++){
	deltaRHad->Fill(ROOT::Math::VectorUtil::DeltaR(gentauHad[j],pfTauRef->p4()));
	if(ROOT::Math::VectorUtil::DeltaR(gentauHad[j],pfTauRef->p4())<deltaR_) Nmatching=Nmatching+1;
	if(ROOT::Math::VectorUtil::DeltaR(gentauHad[j],pfTauRef->p4())>0.1 && ROOT::Math::VectorUtil::DeltaR(gentauHad[j],pfTauRef->p4())<0.35){
	  DRGenTauPt->Fill(gentauHad[j].pt());
	  DRRecoTauPt->Fill(pfTauRef->pt());
	  DRGenTauEta->Fill(gentauHad[j].eta());
	  DRRecoTauEta->Fill(pfTauRef->eta());
	  if(muonIndex!=-1) {
	    const reco::Muon& muon(muH->at(muonIndex));
	    DRMuonTau->Fill(ROOT::Math::VectorUtil::DeltaR(muon.p4(),pfTauRef->p4()));
	  }
	}
      }
    } 
    Nmatch->Fill(Nmatching);

    float deltaR=99.; int tauIndex = -1;
    for(unsigned int i=0;i<tauHandle->size();++i) {
      reco::PFTauRef pfTauRef(tauHandle,i);
      if(!(pfTauRef->pt()>10)) continue;
      if(ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<deltaR_ && ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4())<deltaR){
	deltaR = ROOT::Math::VectorUtil::DeltaR(gentauHad[0],pfTauRef->p4());
	tauIndex = i;
      }
    } 
			
    if(tauIndex==-1) match_2tauh->Fill(10);
    if(tauIndex!=-1) match_2tauh->Fill(11);
    if(tauIndex==-1) match_2tauh_VS_dR->Fill(10,ROOT::Math::VectorUtil::DeltaR(gentau[0],gentau[1]));
    if(tauIndex!=-1) match_2tauh_VS_dR->Fill(11,ROOT::Math::VectorUtil::DeltaR(gentau[0],gentau[1]));
			
    deltaRGenReco_mutau->Fill(deltaR);
    if(tauIndex!=-1){
      genTauPt_mutau->Fill(gentauHad[0].pt());   
      genTauPtWithNeutrino_mutau->Fill(gentauHadWithNeutrino[0].pt());   
      genMuoPt_mutau->Fill(genmuo[0].pt());  
      reco::PFTauRef PFTau(tauHandle,tauIndex);
      float outputDiscmnt1  = (*decayModeFinding)[PFTau];
      float outputDiscmnt2  = (*againstMuonMedium)[PFTau];
      float outputDiscmnt3  = (*againstElectronMedium)[PFTau];
      float outputDiscmnt4  = (*ByVLooseCombinedIsolationDBSumPtCorr)[PFTau];
      float outputDiscmnt5  = (*againstMuonLoose)[PFTau];
      float outputDiscmnt6  = (*againstElectronLoose)[PFTau];
      float outputDiscmnt7  = (*againstElectronLooseMVA3)[PFTau];
      float outputDiscmnt8  = (*againstElectronMediumMVA3)[PFTau];
      float outputDiscmnt9  = (*againstMuonLoose2)[PFTau];
      float outputDiscmnt10 = (*againstMuonMedium2)[PFTau];
      float outputDiscmnt11 = (*againstMuonLoose3)[PFTau];
      float outputDiscmnt12 = (*againstMuonTight3)[PFTau];
      float outputDiscmnt13 = (*byCombinedIsolationDeltaBetaCorrRaw3Hits)[PFTau];
      float outputDiscmnt14 = (*byLooseCombinedIsolationDeltaBetaCorr3Hits)[PFTau];
      float outputDiscmnt15 = (*byMediumCombinedIsolationDeltaBetaCorr3Hits)[PFTau];
			
      recoTauPt_mutau->Fill(PFTau->pt());
      recoTauMass_mutau->Fill(PFTau->mass());
      recoTauDecayMode->Fill(PFTau->decayMode());
      if(PFTau->leadPFCand().isNonnull()) recoTauleadPFCandPt->Fill(PFTau->leadPFCand()->pt());
      tauSelection_mutau->Fill(1);
      if(PFTau->pt()>20) tauSelection_mutau->Fill(2);
      if(PFTau->pt()>20 && outputDiscmnt1>0.5) tauSelection_mutau->Fill(3);
      if(PFTau->pt()>20 && outputDiscmnt1>0.5 && outputDiscmnt5>0.5) tauSelection_mutau->Fill(4);
      if(PFTau->pt()>20 && outputDiscmnt1>0.5 && outputDiscmnt5>0.5 && outputDiscmnt6>0.5) tauSelection_mutau->Fill(5);
      if(PFTau->pt()>20 && outputDiscmnt1>0.5 && outputDiscmnt5>0.5 && outputDiscmnt6>0.5 && outputDiscmnt4>0.5) tauSelection_mutau->Fill(6);
      
      if(muonIndex!=-1) {
	const reco::Muon& muon(muH->at(muonIndex));
	TLorentzVector TauPlusMuon; TauPlusMuon.SetPxPyPzE(PFTau->px()+muon.px(),PFTau->py()+muon.py(),PFTau->pz()+muon.pz(),PFTau->energy()+muon.energy());
	recoTauPlusMuonPt->Fill(TauPlusMuon.Pt());
	if(PFTau->leadPFCand().isNonnull()) recoLeadMinusMuonPtNorm->Fill((muon.pt()-PFTau->leadPFCand()->pt())/(muon.pt()));
	if(PFTau->leadPFCand().isNonnull()) recoLeadMinusMuonPt->Fill(fabs(muon.pt()-PFTau->leadPFCand()->pt()));
      }
 		      
      if(PFTau->pt()>20){
	h_mutau_decayModeFinding->Fill(outputDiscmnt1);
	h_mutau_againstMuonMedium->Fill(outputDiscmnt2);
	h_mutau_againstElectronMedium->Fill(outputDiscmnt3);
	h_mutau_ByVLooseCombinedIsolationDBSumPtCorr->Fill(outputDiscmnt4);
	h_mutau_againstMuonLoose->Fill(outputDiscmnt5);
	h_mutau_againstElectronLoose->Fill(outputDiscmnt6);
	h_mutau_againstElectronLooseMVA3->Fill(outputDiscmnt7);
	h_mutau_againstElectronMediumMVA3->Fill(outputDiscmnt8);
	h_mutau_againstMuonLoose2->Fill(outputDiscmnt9);
	h_mutau_againstMuonMedium2->Fill(outputDiscmnt10);
	h_mutau_againstMuonLoose3->Fill(outputDiscmnt11);
	h_mutau_againstMuonTight3->Fill(outputDiscmnt12);
	h_mutau_byCombinedIsolationDeltaBetaCorrRaw3Hits->Fill(outputDiscmnt13);
	h_mutau_byLooseCombinedIsolationDeltaBetaCorr3Hits->Fill(outputDiscmnt14);
	h_mutau_byMediumCombinedIsolationDeltaBetaCorr3Hits->Fill(outputDiscmnt15);
	ptResol_mutau->Fill((PFTau->pt() - gentauHad[0].pt())/gentauHad[0].pt());
	ptRatio_mutau->Fill(PFTau->pt()/gentauHad[0].pt());
	Response_mutau->Fill(gentauHad[0].pt(), (PFTau->pt() - gentauHad[0].pt())/gentauHad[0].pt());
	ResponseVSdR_mutau->Fill(ROOT::Math::VectorUtil::DeltaR(gentau[0],gentau[1]), (PFTau->pt() - gentauHad[0].pt())/gentauHad[0].pt());
	Response_mutau_profile->Fill(gentauHad[0].pt(), (PFTau->pt() - gentauHad[0].pt())/gentauHad[0].pt());
	if(outputDiscmnt1>0.5 && outputDiscmnt5>0.5 && outputDiscmnt6>0.5 && outputDiscmnt4>0.5){
	  Response_mutau_profile_pt->Fill(gentauHad[0].pt(), (PFTau->pt() - gentauHad[0].pt())/gentauHad[0].pt());
	}
      }
    }
  }

  /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
  */
}


// ------------ method called once each job just before starting event loop  ------------
void 
TauCleaningAnalyzer::beginJob()
{
  
  Service<TFileService> fs;
  
  Ntau = fs->make<TH1D>( "Ntau", "Ntau", 10, -0.5, 9.5 );
  Nele = fs->make<TH1D>( "Nele", "Nele", 10, -0.5, 9.5 );
  Nmuo = fs->make<TH1D>( "Nmuo", "Nmuo", 10, -0.5, 9.5 );
  tauDecay = fs->make<TH1D>( "tauDecay", "tauDecay", 10, -0.5, 9.5 );
  X_pt      = fs->make<TH1D>( "X_pt", "X_pt", 1000, 0, 1000 );
  X_phi     = fs->make<TH1D>( "X_phi", "X_phi", 800, -4, 4 );
  X_eta     = fs->make<TH1D>( "X_eta", "X_eta", 1200, -6, 6 );
  X_mass    = fs->make<TH1D>( "X_mass", "X_mass", 1000, massMin_, massMax_);
  genTauPt1   = fs->make<TH1D>( "genTauPt1", "genTauPt1", 1000, 0, 1000 );
  genTauPt2   = fs->make<TH1D>( "genTauPt2", "genTauPt2", 1000, 0, 1000 );
  ele_pt   = fs->make<TH1D>( "ele_pt", "ele_pt", 1000, 0, 1000 );
  muo_pt   = fs->make<TH1D>( "muo_pt", "muo_pt", 1000, 0, 1000 );
  genTauPt_mutau  = fs->make<TH1D>( "genTauPt_mutau", "genTauPt_mutau", 1000, 0, 1000 );
  genTauPtWithNeutrino_mutau  = fs->make<TH1D>( "genTauPtWithNeutrino_mutau","genTauPtWithNeutrino_mutau", 1000, 0, 1000 );
  genMuoPt_mutau  = fs->make<TH1D>( "genMuoPt_mutau", "genMuoPt_mutau", 1000, 0, 1000 );
  genTauPhi1  = fs->make<TH1D>( "genTauPhi1", "genTauPhi1", 800, -4, 4 );
  genTau1_eta  = fs->make<TH1D>( "genTau1_eta", "genTau1_eta", 1200, -6, 6 );
  genTau2_phi  = fs->make<TH1D>( "genTau2_phi", "genTau2_phi", 800, -4, 4 );
  genTau2_eta  = fs->make<TH1D>( "genTau2_eta", "genTau2_eta", 1200, -6, 6 );
  genTau_DeltaR    = fs->make<TH1D>( "genTau_DeltaR", "genTau_DeltaR", 1000, 0, 10 );
  genTau_DeltaRWithNeutrino    = fs->make<TH1D>( "genTau_DeltaRWithNeutrino", "genTau_DeltaRWithNeutrino", 1000, 0, 10 );
  genTau_DeltaPhi  = fs->make<TH1D>( "genTau_DeltaPhi", "genTau_DeltaPhi", 800, -4, 4 );
  genTau_DeltaEta  = fs->make<TH1D>( "genTau_DeltaEta", "genTau_DeltaEta", 1200, -6, 6 );
  Ztau_pt   = fs->make<TH1D>( "Ztau_pt", "Ztau_pt", 1000, 0, 1000 );
  Ztau_phi  = fs->make<TH1D>( "Ztau_phi", "Ztau_phi", 800, -4, 4 );
  Ztau_eta  = fs->make<TH1D>( "Ztau_eta", "Ztau_eta", 1200, -6, 6 );
  Ztau_mass = fs->make<TH1D>( "Ztau_mass", "Ztau_mass", 2000, 0, 200 );

  VertexFake = fs->make<TH1D>("VertexFake", "VertexFake", 4, -0.5, 3.5);
  VertexRho = fs->make<TH1D>("VertexRho", "VertexRho", 1000, -5, 5);
  VertexNdof = fs->make<TH1D>("VertexNdof", "VertexNdof", 100, 0, 100);
  VertexZ = fs->make<TH1D>("VertexZ", "VertexZ", 200, 0,  20);

  h_mutau_decayModeFinding = fs->make<TH1D>( "h_mutau_decayModeFinding", "h_mutau_decayModeFinding", 4, -1.5, 2.5 );
  h_mutau_againstMuonMedium = fs->make<TH1D>( "h_mutau_againstMuonMedium", "h_mutau_againstMuonMedium", 4, -1.5, 2.5 );
  h_mutau_againstElectronMedium = fs->make<TH1D>( "h_mutau_againstElectronMedium", "h_mutau_againstElectronMedium", 4, -1.5, 2.5 );
  h_mutau_ByVLooseCombinedIsolationDBSumPtCorr = fs->make<TH1D>( "h_mutau_ByVLooseCombinedIsolationDBSumPtCorr", "h_mutau_ByVLooseCombinedIsolationDBSumPtCorr", 4,-1.5, 2.5 );
  h_mutau_againstMuonLoose=fs->make<TH1D>("h_mutau_againstMuonLoose","h_mutau_againstMuonLoose", 4, -1.5, 2.5 );
  h_mutau_againstElectronLoose=fs->make<TH1D>("h_mutau_againstElectronLoose","h_mutau_againstElectronLoose", 4, -1.5, 2.5 );
  h_mutau_againstElectronLooseMVA3=fs->make<TH1D>("h_mutau_againstElectronLooseMVA3","h_mutau_againstElectronLooseMVA3", 4, -1.5, 2.5 );
  h_mutau_againstElectronMediumMVA3=fs->make<TH1D>("h_mutau_againstElectronMediumMVA3","h_mutau_againstElectronMediumMVA3", 4, -1.5, 2.5 );
  h_mutau_againstMuonLoose2=fs->make<TH1D>("h_mutau_againstMuonLoose2","h_mutau_againstMuonLoose2", 4, -1.5, 2.5 );
  h_mutau_againstMuonMedium2=fs->make<TH1D>("h_mutau_againstMuonMedium2","h_mutau_againstMuonMedium2", 4, -1.5, 2.5 );
  h_mutau_againstMuonLoose3=fs->make<TH1D>("h_mutau_againstMuonLoose3","h_mutau_againstMuonLoose3", 4, -1.5, 2.5 );
  h_mutau_againstMuonTight3=fs->make<TH1D>("h_mutau_againstMuonTight3","h_mutau_againstMuonTight3", 4, -1.5, 2.5 );
  h_mutau_byCombinedIsolationDeltaBetaCorrRaw3Hits=fs->make<TH1D>("h_mutau_byCombinedIsolationDeltaBetaCorrRaw3Hits","h_mutau_byCombinedIsolationDeltaBetaCorrRaw3Hits", 4, -1.5, 2.5 );
  h_mutau_byLooseCombinedIsolationDeltaBetaCorr3Hits=fs->make<TH1D>("h_mutau_byLooseCombinedIsolationDeltaBetaCorr3Hits","h_mutau_byLooseCombinedIsolationDeltaBetaCorr3Hits", 4, -1.5, 2.5 );
  h_mutau_byMediumCombinedIsolationDeltaBetaCorr3Hits=fs->make<TH1D>("h_mutau_byMediumCombinedIsolationDeltaBetaCorr3Hits","h_mutau_byMediumCombinedIsolationDeltaBetaCorr3Hits", 4, -1.5, 2.5 );

  deltaRHad = fs->make<TH1D>("deltaRHad", "deltaRHad", 1000, 0, 10 );
  deltaRGenReco_mutau = fs->make<TH1D>("deltaRGenReco_mutau", "deltaRGenReco_mutau", 1000, 0, 10 );

  recoTauPt_mutau = fs->make<TH1D>( "recoTauPt_mutau", "recoTauPt_mutau", 1000, 0, 1000 );
  recoTauMass_mutau = fs->make<TH1D>( "recoTauMass_mutau", "recoTauMass_mutau", 5000, 0, 50);
  recoTauDecayMode = fs->make<TH1D>("recoTauDecayMode", "recoTauDecayMode", 21, -0.5, 20.5);
  recoTauleadPFCandPt = fs->make<TH1D>("recoTauleadPFCandPt", "recoTauleadPFCandPt", 1000, 0, 1000 );
  recoTauPlusMuonPt = fs->make<TH1D>("recoTauPlusMuonPt", "recoTauPlusMuonPt", 1000, 0, 1000 );
  recoLeadMinusMuonPtNorm = fs->make<TH1D>("recoLeadMinusMuonPtNorm", "recoLeadMinusMuonPtNorm", 2000, -10, 10 );
  recoLeadMinusMuonPt = fs->make<TH1D>("recoLeadMinusMuonPt", "recoLeadMinusMuonPt", 1000, 0, 1000 );
  recoMuonPt = fs->make<TH1D>( "recoMuonPt", "recoMuonPt", 1000, 0, 1000 );

  ptResol_mutau = fs->make<TH1D>( "ptResol_mutau", "ptResol_mutau", 4000, -20, 20 );
  ptRatio_mutau = fs->make<TH1D>( "ptRatio_mutau", "ptRatio_mutau", 4000, -20, 20 );
  Response_mutau = fs->make<TH2D>( "Response_mutau", "Response_mutau", 1000,0,1000, 4000, -20, 20 );
  Response_mutau_profile  = fs->make<TProfile>("Response_mutau_profile","Response_mutau_profile",100,0,1000,20,20);
  Response_mutau_profile_pt  = fs->make<TProfile>("Response_mutau_profile_pt","Response_mutau_profile_pt",100,0,1000,20,20);
  ResponseVSdR_mutau = fs->make<TH2D>( "ResponseVSdR_mutau", "ResponseVSdR_mutau", 1000,0,10, 4000,-20,20);
  tauSelection_mutau = fs->make<TH1D>( "tauSelection_mutau", "tauSelection_mutau", 12, -0.5, 11.5 );

  match_2tauh = fs->make<TH1D>("match_2tauh", "match_2tauh", 16, -0.5, 15.5 );
  match_2tauh_VS_dR = fs->make<TH2D>("match_2tauh_VS_dR", "match_2tauh_VS_dR", 16, -0.5, 15.5, 1000, 0, 10);


  DRGenTauPt   = fs->make<TH1D>("DRGenTauPt",  "DRGenTauPt",   1000, 0, 1000);
  DRRecoTauPt  = fs->make<TH1D>("DRRecoTauPt", "DRRecoTauPt",  1000, 0, 1000);
  DRGenTauEta  = fs->make<TH1D>("DRGenTauEta", "DRGenTauEta",  120, -6, 6);
  DRRecoTauEta = fs->make<TH1D>("DRRecoTauEta","DRRecoTauEta", 120, -6, 6);
  DRMuonTau = fs->make<TH1D>("DRMuonTau", "DRMuonTau", 600, 0, 6 );
  Nmatch = fs->make<TH1D>("Nmatch", "Nmatch", 11, -0.5, 10.5);
 
  return;
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauCleaningAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TauCleaningAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauCleaningAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauCleaningAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauCleaningAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauCleaningAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauCleaningAnalyzer);
