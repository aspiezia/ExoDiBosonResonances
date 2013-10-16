import FWCore.ParameterSet.Config as cms
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ('PoolSource',fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_10_1_SsM.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_11_1_Ng2.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_1_1_fJS.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_2_1_0iG.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_3_1_oQO.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_4_1_T0d.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_5_1_Elm.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_6_1_IgP.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_7_1_Agu.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_8_1_30G.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_20131015_RSGZZPythia_20131015_090830/bonato/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola/EDBR_PATtuple_edbr_20131015_RSGZZPythia/1b325ddfb984c14533be7920e22baeef/RSGravitonToZZToLLJJ_kMpl005_M-1000_TuneZ2star_8TeV-pythia6-tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_9_1_eeY.root' 
]);
