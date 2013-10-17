import FWCore.ParameterSet.Config as cms

ak5PFJetsNoMu = cms.EDProducer(
    "CleanJetsProducer"
    )

CleanJetsSequence = cms.Sequence(ak5PFJetsNoMu)
