// -*- C++ -*-
//
// Package:    PtStudy
// Class:      PtStudy
// 
/**\class PtStudy PtStudy.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/PtStudy.cc

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

class PtStudy : public edm::EDAnalyzer {
public:
  explicit PtStudy(const edm::ParameterSet&);
  ~PtStudy();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
 
  //HISTOGRAMS
  TH1D* elePt1;
  TH1D* elePt2;
  TH1D* cutElePt1;
  TH1D* cutElePt2;
  TH1D* cutEleMatching1;
  TH1D* cutEleMatching2;
  TH1D* cutEleMatching;

  TH1D* muoPt1;
  TH1D* muoPt2;
  TH1D* cutMuoPt1;
  TH1D* cutMuoPt2;
  TH1D* cutMuoMatching1;
  TH1D* cutMuoMatching2;
  TH1D* cutMuoMatching;

  TH1D* muoEMPt;
  TH1D* eleEMPt;
  TH1D* cutMuoEMPt;
  TH1D* cutEleEMPt;
  TH1D* cutMuoEMMatching;
  TH1D* cutEleEMMatching;
  TH1D* cutEleMuoEMMatching;

  TH1D* metPtDiEle;
  TH1D* cutMetPtDiEle;
  TH1D* metPtDiMuo;
  TH1D* cutMetPtDiMuo;
  TH1D* metPtEleMuo;
  TH1D* cutMetPtEleMuo;
  TH1D* metPtEleTau;
  TH1D* cutMetPtEleTau;
  TH1D* metPtMuoTau;
  TH1D* cutMetPtMuoTau;

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
PtStudy::PtStudy(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
}


PtStudy::~PtStudy()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PtStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel("primaryVertexFilter", vertices);
   reco::Vertex primaryVertex;
   primaryVertex = vertices->at(0);

   edm::Handle<pat::METCollection> met;
   iEvent.getByLabel("patMETs", met);

   edm::Handle<pat::ElectronCollection> eleH;
   iEvent.getByLabel("patElectronsWithTrigger", eleH);

   edm::Handle<pat::MuonCollection> muoH;
   iEvent.getByLabel("patMuonsWithTrigger", muoH);
   
   Handle<vector<reco::GenParticle> > genParts;
   iEvent.getByLabel("genParticles", genParts);
  
   int ele = 0; int muo = 0;
   vector<math::PtEtaPhiELorentzVector> genele;
   vector<math::PtEtaPhiELorentzVector> genmuo;
   for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
     const reco::GenParticle & genPart = (*genParts)[ngenPart];
     if(abs(genPart.pdgId())==15 && genPart.status()!=3){
       for(unsigned int ndaugh=0; ndaugh<genPart.numberOfDaughters(); ndaugh++){
	 const reco::Candidate * daughter = genPart.daughter(ndaugh);
	 if(abs(daughter->pdgId())==11 && daughter->status()==1) ele = ele + 1;
	 if(abs(daughter->pdgId())==13 && daughter->status()==1) muo = muo + 1;
	 math::PtEtaPhiELorentzVector gen_prov; gen_prov=daughter->p4();
	 if(abs(daughter->pdgId())==11 && daughter->status()==1) genele.push_back(gen_prov);
	 if(abs(daughter->pdgId())==13 && daughter->status()==1) genmuo.push_back(gen_prov);
       }
     }
   }

   if(ele==2 && muo==0){
     bool matching1=false;bool matching2=false;
     pat::ElectronCollection::const_iterator SelectedEle1;
     pat::ElectronCollection::const_iterator SelectedEle2;
     float DRele1=99.; float DRele2=99.;
     for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
       if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0])<DRele1 && ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0])<0.3){
	 DRele1 = ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0]);
	 SelectedEle1 = electron;
	 matching1=true;
       }
       if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[1])<DRele2 && ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[1])<0.3){
	 DRele2 = ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[1]);
	 SelectedEle2 = electron;
	 matching2=true;
       }
     }
     if(matching1 && matching2){
       if(SelectedEle1->pt()>SelectedEle2->pt()){
	 elePt1->Fill(SelectedEle1->pt());
	 elePt2->Fill(SelectedEle2->pt());
	 cutElePt1->Fill(SelectedEle1->pt()>20);
	 cutElePt2->Fill(SelectedEle2->pt()>10);
       } else {
	 elePt1->Fill(SelectedEle2->pt());
	 elePt2->Fill(SelectedEle1->pt());
	 cutElePt1->Fill(SelectedEle2->pt()>20);
	 cutElePt2->Fill(SelectedEle1->pt()>10);
       }
     }
     cutEleMatching1->Fill(matching1);
     cutEleMatching2->Fill(matching2);
     cutEleMatching->Fill(matching1&&matching2);
     metPtDiEle->Fill(met->begin()->pt());
     cutMetPtDiEle->Fill(met->begin()->pt()>100);
   }

   if(ele==0 && muo==2){
     bool matching1=false;bool matching2=false;
     pat::MuonCollection::const_iterator SelectedMuo1;
     pat::MuonCollection::const_iterator SelectedMuo2;
     float DRmuo1=99.; float DRmuo2=99.;
     for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
       if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0])<DRmuo1 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0])<0.3){
	 DRmuo1 = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0]);
	 SelectedMuo1 = muon;
	 matching1=true;
       }
       if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[1])<DRmuo2 && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[1])<0.3){
	 DRmuo2 = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[1]);
	 SelectedMuo2 = muon;
	 matching2=true;
       }
     }
     if(matching1 && matching2){
       if(SelectedMuo1->pt()>SelectedMuo2->pt()){
	 muoPt1->Fill(SelectedMuo1->pt());
	 muoPt2->Fill(SelectedMuo2->pt());
	 cutMuoPt1->Fill(SelectedMuo1->pt()>20);
	 cutMuoPt2->Fill(SelectedMuo2->pt()>10);
       } else {
	 muoPt1->Fill(SelectedMuo2->pt());
	 muoPt2->Fill(SelectedMuo1->pt());
	 cutMuoPt1->Fill(SelectedMuo2->pt()>20);
	 cutMuoPt2->Fill(SelectedMuo1->pt()>10);
       }
     }
     cutMuoMatching1->Fill(matching1);
     cutMuoMatching2->Fill(matching2);
     cutMuoMatching->Fill(matching1&&matching2);
     metPtDiMuo->Fill(met->begin()->pt());
     cutMetPtDiMuo->Fill(met->begin()->pt()>100);
   }

   if(ele==1 && muo==1){
     bool matching1=false;bool matching2=false;
     pat::ElectronCollection::const_iterator SelectedEle;
     pat::MuonCollection::const_iterator SelectedMuo;
     float DRele=99.; float DRmuo=99.;
     for(pat::ElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {
       if(ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0])<DRele && ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0])<0.3){
	 DRele = ROOT::Math::VectorUtil::DeltaR(electron->p4(),genele[0]);
	 SelectedEle = electron;
	 matching1=true;
       }
     }
     for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
       if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0])<DRmuo && ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0])<0.3){
	 DRmuo = ROOT::Math::VectorUtil::DeltaR(muon->p4(),genmuo[0]);
	 SelectedMuo = muon;
	 matching2=true;
       }
     }
     if(matching1 && matching2){
       eleEMPt->Fill(SelectedEle->pt());
       muoEMPt->Fill(SelectedMuo->pt());
       cutEleEMPt->Fill(SelectedEle->pt()>10);
       cutMuoEMPt->Fill(SelectedMuo->pt()>10);
     }
     cutEleEMMatching->Fill(matching1);
     cutMuoEMMatching->Fill(matching2);
     cutEleMuoEMMatching->Fill(matching1&&matching2);
     metPtEleMuo->Fill(met->begin()->pt());
     cutMetPtEleMuo->Fill(met->begin()->pt()>100);
   }

   if(ele==0 && muo==1){
     metPtMuoTau->Fill(met->begin()->pt());
     cutMetPtMuoTau->Fill(met->begin()->pt()>50);
   }

   if(ele==1 && muo==0){
     metPtEleTau->Fill(met->begin()->pt());
     cutMetPtEleTau->Fill(met->begin()->pt()>50);
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
PtStudy::beginJob()
{
  Service<TFileService> fs;
  elePt1          = fs->make<TH1D>("elePt1",          "elePt1",         2000, 0, 2000);
  elePt2          = fs->make<TH1D>("elePt2",          "elePt2",         2000, 0, 2000);
  cutElePt1       = fs->make<TH1D>("cutElePt1",       "cutElePt1",       2, -0.5, 1.5);
  cutElePt2       = fs->make<TH1D>("cutElePt2",       "cutElePt2",       2, -0.5, 1.5);
  cutEleMatching1 = fs->make<TH1D>("cutEleMatching1", "cutEleMatching1", 2, -0.5, 1.5);
  cutEleMatching2 = fs->make<TH1D>("cutEleMatching2", "cutEleMatching2", 2, -0.5, 1.5);
  cutEleMatching  = fs->make<TH1D>("cutEleMatching",  "cutEleMatching",  2, -0.5, 1.5);

  muoPt1          = fs->make<TH1D>("muoPt1",          "muoPt1",         2000, 0, 2000);
  muoPt2          = fs->make<TH1D>("muoPt2",          "muoPt2",         2000, 0, 2000);
  cutMuoPt1       = fs->make<TH1D>("cutMuoPt1",       "cutMuoPt1",       2, -0.5, 1.5);
  cutMuoPt2       = fs->make<TH1D>("cutMuoPt2",       "cutMuoPt2",       2, -0.5, 1.5);
  cutMuoMatching1 = fs->make<TH1D>("cutMuoMatching1", "cutMuoMatching1", 2, -0.5, 1.5);
  cutMuoMatching2 = fs->make<TH1D>("cutMuoMatching2", "cutMuoMatching2", 2, -0.5, 1.5);
  cutMuoMatching  = fs->make<TH1D>("cutMuoMatching",  "cutMuoMatching",  2, -0.5, 1.5);

  eleEMPt              = fs->make<TH1D>("eleEMPt",             "eleEMPt",            2000, 0, 2000);
  muoEMPt              = fs->make<TH1D>("muoEMPt",             "muoEMPt",            2000, 0, 2000);
  cutEleEMPt           = fs->make<TH1D>("cutEleEMPt",          "cutEleEMPt",          2, -0.5, 1.5);
  cutMuoEMPt           = fs->make<TH1D>("cutMuoEMPt",          "cutMuoEMPt",          2, -0.5, 1.5);
  cutEleEMMatching     = fs->make<TH1D>("cutEleEMMatching",    "cutEleEMMatching",    2, -0.5, 1.5);
  cutMuoEMMatching     = fs->make<TH1D>("cutMuoMatching",      "cutMuoEMMatching",    2, -0.5, 1.5);
  cutEleMuoEMMatching  = fs->make<TH1D>("cutEleMuoEMMatching", "cutEleMuoEMMatching", 2, -0.5, 1.5);

  metPtDiEle      = fs->make<TH1D>("metPtDiEle",      "metPtDiEle",     2000, 0, 2000);
  cutMetPtDiEle   = fs->make<TH1D>("cutMetPtDiEle",   "cutMetPtDiEle",   2, -0.5, 1.5);
  metPtDiMuo      = fs->make<TH1D>("metPtDiMuo",      "metPtDiMuo",     2000, 0, 2000);
  cutMetPtDiMuo   = fs->make<TH1D>("cutMetPtDiMuo",   "cutMetPtDiMuo",   2, -0.5, 1.5);
  metPtEleMuo     = fs->make<TH1D>("metPtEleMuo",     "metPtEleMuo",    2000, 0, 2000);
  cutMetPtEleMuo  = fs->make<TH1D>("cutMetPtEleMuo",  "cutMetPtEleMuo",  2, -0.5, 1.5);
  metPtEleTau     = fs->make<TH1D>("metPtEleTau",     "metPtEleTau",    2000, 0, 2000);
  cutMetPtEleTau  = fs->make<TH1D>("cutMetPtEleTau",  "cutMetPtEleTau",  2, -0.5, 1.5);
  metPtMuoTau     = fs->make<TH1D>("metPtMuoTau",     "metPtMuoTau",    2000, 0, 2000);
  cutMetPtMuoTau  = fs->make<TH1D>("cutMetPtMuoTau",  "cutMetPtMuoTau",  2, -0.5, 1.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PtStudy::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
PtStudy::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PtStudy::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PtStudy::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PtStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PtStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(PtStudy);
