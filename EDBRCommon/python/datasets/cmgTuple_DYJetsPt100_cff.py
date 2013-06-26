import FWCore.ParameterSet.Config as cms

cmgFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",
                    noEventSort = cms.untracked.bool(True),
                    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                    fileNames = cmgFiles
                   )

cmgFiles.extend([
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_0.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_1.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_10.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_11.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_12.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_13.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_14.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_15.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_16.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_17.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_18.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_19.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_2.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_20.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_21.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_22.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_23.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_24.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_25.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_26.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_27.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_28.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_29.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_3.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_30.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_31.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_32.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_33.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_34.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_35.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_36.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_37.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_38.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_39.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_4.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_40.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_41.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_42.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_43.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_44.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_45.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_46.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_47.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_48.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_49.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_5.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_50.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_51.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_52.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_53.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_54.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_55.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_56.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_57.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_58.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_59.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_6.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_60.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_61.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_62.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_63.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_64.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_7.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_8.root',
    '/store/cmst3/group/exovv/CMGtuple/productionv2b/Summer12/DYJetsPt100/cmgTuple_9.root',
    ])
