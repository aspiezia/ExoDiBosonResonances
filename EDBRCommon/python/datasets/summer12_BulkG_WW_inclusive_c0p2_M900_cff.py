import FWCore.ParameterSet.Config as cms

cmgFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",
                    noEventSort = cms.untracked.bool(True),
                    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                    fileNames = cmgFiles
                   )

cmgFiles.extend([
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_10_1_Iza.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_11_1_VOQ.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_12_1_J8E.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_13_1_jrK.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_14_1_KUF.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_15_1_YyT.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_16_1_iYp.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_17_1_ggK.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_18_1_CBU.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_19_1_Fwq.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_1_1_2MN.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_20_1_mRs.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_21_1_9Ci.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_2_1_uoS.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_3_1_VhN.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_4_1_r1c.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_5_1_QAn.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_6_1_Bww.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_7_1_UgE.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_8_1_5CB.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M900_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M900_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M900_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_9_1_9Ju.root',
    ])
