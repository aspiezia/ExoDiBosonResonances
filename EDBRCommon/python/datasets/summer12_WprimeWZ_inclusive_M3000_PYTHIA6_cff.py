import FWCore.ParameterSet.Config as cms
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ('PoolSource',fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_zz_Summer12MC_WprimeWZ_PYTHIA_20131013_190656/bonato/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6/EDBR_PATtuple_edbr_zz_20131013_Summer12MC_WprimeWZ_PYTHIA/1b325ddfb984c14533be7920e22baeef/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_1_1_iSM.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_zz_Summer12MC_WprimeWZ_PYTHIA_20131013_190656/bonato/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6/EDBR_PATtuple_edbr_zz_20131013_Summer12MC_WprimeWZ_PYTHIA/1b325ddfb984c14533be7920e22baeef/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_2_1_UTg.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_zz_Summer12MC_WprimeWZ_PYTHIA_20131013_190656/bonato/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6/EDBR_PATtuple_edbr_zz_20131013_Summer12MC_WprimeWZ_PYTHIA/1b325ddfb984c14533be7920e22baeef/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_3_1_AOJ.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_zz_Summer12MC_WprimeWZ_PYTHIA_20131013_190656/bonato/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6/EDBR_PATtuple_edbr_zz_20131013_Summer12MC_WprimeWZ_PYTHIA/1b325ddfb984c14533be7920e22baeef/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_4_1_TbV.root',
       '/store/cmst3/group/exovv/bonato/EDBR_PATtuple_edbr_zz_Summer12MC_WprimeWZ_PYTHIA_20131013_190656/bonato/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6/EDBR_PATtuple_edbr_zz_20131013_Summer12MC_WprimeWZ_PYTHIA/1b325ddfb984c14533be7920e22baeef/WprimeToWZ_M-3000_TuneZ2star_8TeV-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM_5_1_0zE.root'
]);
