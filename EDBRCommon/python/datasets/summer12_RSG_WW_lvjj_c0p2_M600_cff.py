import FWCore.ParameterSet.Config as cms

cmgFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",
                    noEventSort = cms.untracked.bool(True),
                    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                    fileNames = cmgFiles
                   )

cmgFiles.extend([
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_101_1_0Rm.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_10_1_t2F.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_11_1_caC.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_12_1_FXr.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_13_1_2L0.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_14_1_y6O.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_15_1_uxM.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_16_1_YTQ.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_17_1_fgc.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_18_1_MOP.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_19_1_FCW.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_1_1_5VJ.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_20_1_bTp.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_21_1_Dmk.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_22_1_Iev.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_23_1_0Nk.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_24_1_3V4.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_25_1_AIs.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_26_1_J1D.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_27_1_2yp.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_29_1_0Wn.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_2_1_0yH.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_30_1_j2c.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_31_1_Oum.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_32_1_dQ2.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_33_1_n28.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_34_1_JgL.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_35_1_h22.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_36_1_sG1.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_37_1_gCP.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_38_1_mBY.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_39_1_9es.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_3_1_Lwu.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_40_1_2nR.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_41_1_YgR.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_42_1_8wC.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_43_1_BW3.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_45_1_n3h.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_46_1_GRh.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_47_1_qJC.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_48_1_B0C.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_49_1_r7C.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_4_1_SEu.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_50_1_I1i.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_51_1_f2z.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_52_1_iua.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_53_1_IkN.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_54_1_QLw.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_55_1_VH0.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_56_1_RaL.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_57_1_iiq.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_58_1_3fE.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_59_1_A8A.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_60_1_uKC.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_61_1_7kO.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_65_1_Bkb.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_66_1_s4m.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_68_1_fxJ.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_69_1_Gmw.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_70_1_oA9.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_71_1_Zqh.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_72_1_xAN.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_73_1_Fpz.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_74_1_VEd.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_75_1_P0i.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_76_1_Uyc.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_77_1_oCv.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_78_1_QI1.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_79_1_BF4.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_7_1_GR2.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_80_1_9bJ.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_81_1_s58.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_82_1_USM.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_85_1_8cG.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_86_1_1LF.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_87_1_ULf.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_89_1_htm.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_8_1_Y9l.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_90_1_xrH.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_91_1_5T5.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_92_1_S0A.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_95_1_rUg.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_96_1_qDw.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_97_1_537.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_98_1_MPT.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_99_1_KvM.root',
    '/store/cmst3/group/exovv/shuai/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_RSGravitonsWW_20131111_042217/shuai/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/RSGraviton_WWlvjj_kMpl02_M-600_TuneZ2star_8TeV-MG__qili-RSWW_600_02_AODSIM-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_9_1_70C.root',
    ])
