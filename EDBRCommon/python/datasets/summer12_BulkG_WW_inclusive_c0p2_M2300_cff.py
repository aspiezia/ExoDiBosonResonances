import FWCore.ParameterSet.Config as cms

cmgFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",
                    noEventSort = cms.untracked.bool(True),
                    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                    fileNames = cmgFiles
                   )

cmgFiles.extend([
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_10_1_KGO.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_11_1_SZh.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_12_1_aSx.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_13_1_kmQ.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_14_1_wER.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_15_1_bkY.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_16_1_ClT.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_17_1_Bku.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_18_1_y3z.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_19_1_j8u.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_1_1_cA7.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_20_1_RYF.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_21_1_jUC.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_2_1_dod.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_3_1_G7u.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_4_1_k47.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_5_1_la1.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_6_1_Sff.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_7_1_AAq.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_8_1_cmI.root',
    '/store/cmst3/group/exovv/santanas/EDBR_PATtuple_edbr_zz_20130605_Summer12MC_BulkGravitonsWW_20130916_175629/santanas/BulkG_WW_inclusive_c0p2_M2300_GENSIM/EDBR_PATtuple_edbr_zz_20130605/1b325ddfb984c14533be7920e22baeef/BulkG_WW_inclusive_c0p2_M2300_GENSIM__shuai-BulkG_WW_inclusive_c0p2_M2300_AODSIM-2c74483358b1f8805e5601fc325d256c__USER_9_1_A80.root',
    ])
