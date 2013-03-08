import FWCore.ParameterSet.Config as cms

cmgFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",
                                                noEventSort = cms.untracked.bool(True),
                                                duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                                                fileNames = cmgFiles
                                                )

cmgFiles.extend([

       '/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/CMGtuple/santanas/production07032013_edbr_vv_20130313_plusBtagFix/Summer12/CA8//RSG_WW_lvjj_c0p2_M600_xww/cmgTuple_0.root',
       '/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/CMGtuple/santanas/production07032013_edbr_vv_20130313_plusBtagFix/Summer12/CA8//RSG_WW_lvjj_c0p2_M600_xww/cmgTuple_1.root',
       '/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/CMGtuple/santanas/production07032013_edbr_vv_20130313_plusBtagFix/Summer12/CA8//RSG_WW_lvjj_c0p2_M600_xww/cmgTuple_2.root',
       '/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/CMGtuple/santanas/production07032013_edbr_vv_20130313_plusBtagFix/Summer12/CA8//RSG_WW_lvjj_c0p2_M600_xww/cmgTuple_3.root',
       '/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/CMGtuple/santanas/production07032013_edbr_vv_20130313_plusBtagFix/Summer12/CA8//RSG_WW_lvjj_c0p2_M600_xww/cmgTuple_4.root',
    ])
