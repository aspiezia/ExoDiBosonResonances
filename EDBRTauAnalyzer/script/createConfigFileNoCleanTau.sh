#! /bin/bash

index_min=$1
index_max=$2
index_OutputFile=$3
dataset=$4
index=0

echo "import FWCore.ParameterSet.Config as cms"
echo ""
echo 'process = cms.Process("Demo")'
echo ""
echo 'process.load("FWCore.MessageService.MessageLogger_cfi")'
echo "process.MessageLogger.cerr.FwkReport.reportEvery = 10000"
echo ""
echo "process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )"
echo ""
echo 'process.source = cms.Source("PoolSource",'
echo "    # replace 'myfile.root' with the source file you want to use"
echo "    fileNames = cms.untracked.vstring("


while read LINE
  do
  let index=$index+1
  if [  "$index" -lt "$index_max" ] && [ "$index" -ge "$index_min" ]
  then
  lin=" $LINE "
  echo $LINE
  fi
done < $dataset\_fileList.txt

echo "    )"
echo ")"
echo ""
echo "process.demo = cms.EDAnalyzer('TauCleaningAnalyzer',"
echo "                              deltaR = cms.untracked.double(0.4),"
echo "                              massMin = cms.untracked.double(1900),"
echo "                              massMax = cms.untracked.double(2100),"
echo "                              signal = cms.untracked.bool(True)"
echo ")"
echo "#process.demo = cms.EDAnalyzer('Test'"
echo "#)"
echo ""
echo 'process.TFileService = cms.Service("TFileService",'
echo "        fileName = cms.string('output$index_OutputFile.root')"
echo ")"
echo ""
echo '#process.load("ExoDiBosonResonances.EDBRTauAnalyzer.CleanJets_cff")'
echo "#process.p = cms.Path(process.CleanJetsSequence + process.demo)"
echo 'process.load("Configuration.Geometry.GeometryIdeal_cff")'
echo "process.load('Configuration.StandardSequences.Services_cff')"
echo "process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')"
echo "process.GlobalTag.globaltag = 'START53_V7F::All'"
echo 'process.load("Configuration.StandardSequences.MagneticField_cff")'
echo 'process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")'
echo "process.hpsPFTauDiscriminationByLooseMuonRejection3.dRmuonMatch = cms.double(0.05)"
echo "process.hpsPFTauDiscriminationByTightMuonRejection3.dRmuonMatch = cms.double(0.05)"
echo ""
echo 'process.primaryVertexFilter = cms.EDFilter("VertexSelector",'
echo '    src = cms.InputTag("offlinePrimaryVertices"),'
echo '    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"),'
echo "    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection."
echo ")"
echo ""
echo 'process.load("ExoDiBosonResonances.EDBRTauAnalyzer.CA8Jet_cff")'
echo "process.p = cms.Path(process.primaryVertexFilter + process.CA8JetsSequence + process.PFTau * process.demo)"
echo "process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )"
