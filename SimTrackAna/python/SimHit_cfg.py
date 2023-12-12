
import FWCore.ParameterSet.Config as cms

# from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
# process = cms.Process('PROD',Phase2C11)

#process.load("SimGeneral.HepPDTESSource.pdt_cfi")
#process.load("Configuration.Geometry.GeometryExtended2026D76Reco_cff")

#from Configuration.Eras.Modifier_phase2_hgcalV12_cff import phase2_hgcalV12
#process = cms.Process('PROD',phase2_hgcalV12)

# from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
# process = cms.Process('PROD',Phase2C9)

# from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
# process = cms.Process('PROD',Phase2C11I13M9)

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D83Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.source = cms.Source("PoolSource",
                            dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
                            #fileNames = cms.untracked.vstring('file:/home/indra/Data/store/relval/CMSSW_13_0_0_pre2/23293.0_CloseByParticleGun+2026D88_NoNoise/step2.root'),
                            #fileNames = cms.untracked.vstring('file:/home/indra/Data/store/relval/CMSSW_13_0_0_pre2/Single-Muon-2026D88_NoNoise_noise_MIP_sdPixels/step2.root'),
                            fileNames = cms.untracked.vstring('file:/home/indra/Data/store/relval/CMSSW_13_0_0_pre2/Single-Muon-2026D88_NoNoise_noise_MIP/step2.root'),
                            #fileNames = cms.untracked.vstring('file:/home/indra/Data/store/relval/CMSSW_13_0_0_pre2/Single-Muon-2026D88_NoNoise_sdPixels/step2.root'),
                            inputCommands = cms.untracked.vstring(
                                # 'keep *',
                                # 'drop PSimHits_*_*_*',
                                'drop *',
                                'keep *_g4SimHits__SIM',
                                'keep *_g4SimHits_HGCHitsEE_SIM',
                                'keep *_g4SimHits_HGCHitsHEfront_SIM',
                                'keep *_g4SimHits_HGCHitsHEback_SIM',
                                'keep *_simHGCalUnsuppressedDigis_EE_*',
                                'keep *_simHGCalUnsuppressedDigis_HEfront_*',
                                'keep *_simHGCalUnsuppressedDigis_HEback_*',
                            ),
                            secondaryFileNames = cms.untracked.vstring()
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.prodEE = cms.EDAnalyzer('SimHit',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             Detector   = cms.string("HGCalEESensitive"),
                             digihits = cms.untracked.InputTag("simHGCalUnsuppressedDigis","EE"),
                             SampleIndx = cms.untracked.int32(2),
                         )

process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
    digihits = cms.untracked.InputTag("simHGCalUnsuppressedDigis","HEfront"),
)

process.prodHEB = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
    Detector   = cms.string("HGCalHEScintillatorSensitive"),
    digihits = cms.untracked.InputTag("simHGCalUnsuppressedDigis","HEback"),
)

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('geantoutput.root')
 )

#process.p = cms.Path(process.prodEE*process.prodHEF*process.prodHEB)
#process.p = cms.Path(process.prodEE*process.prodHEF)
process.p = cms.Path(process.prodEE*process.prodHEB)
