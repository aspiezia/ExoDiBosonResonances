import FWCore.ParameterSet.Config as cms
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ('PoolSource',fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_1_1_msq.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_2_1_HPv.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_3_1_4B2.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_4_1_9Vq.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_5_1_9Nt.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_6_1_L5W.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_7_1_jnB.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_8_1_Weh.root',
'/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/santanas/EDBR_PATtuple_edbr_zz_20130113_Summer12MC_BulkGravitonsZZ_20130206_194506/santanas/JHU_Bulk600_ZZ_c1_v2/EDBR_PATtuple_edbr_zz_20130113/a2df48025a27dd8dc244a890c390dec0/JHU_Bulk600_ZZ_c1_v2__tomei-JHU_Bulk600_ZZ_c1_STEP3-c8f8ed334db8a7d6f56c62266b1dfa5b__USER_9_1_bXy.root'
]);