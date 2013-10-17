import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ca8PFJetsCHS = ak5PFJets.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )
jetSource = 'ca8PFJetsCHS'

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
ca8PFJetsCHSpruned = ak5PFJetsPruned.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )
jetSource = 'ca8PFJetsCHSpruned'


pfNoPileUp = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

ca8PFJetsCHSwithNsub = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("ca8PFJetsCHS"),
                              cone=cms.double(0.8)
                              )


CA8JetsSequence = cms.Sequence(pfNoPileUp+ca8PFJetsCHS+ca8PFJetsCHSpruned+ca8PFJetsCHSwithNsub)
