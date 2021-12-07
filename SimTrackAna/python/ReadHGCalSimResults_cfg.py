import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D83Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/user/i/idas/MuDeltaPt100GeV/step1_0.root'
    )
)

process.prodEE = cms.EDAnalyzer('ReadHGCalSimResults',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             Detector   = cms.string("HGCalEESensitive"),
                         )

process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
)

process.prodHEB = process.prodHEF.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
    Detector   = cms.string("HGCalHEScintillatorSensitive"),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('HGCalSimResults.root')
 )

process.p = cms.Path(process.prodEE*process.prodHEF*process.prodHEB)
