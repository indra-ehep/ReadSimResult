###############################################################################
# Way to use this:
#   cmsRun runHGCalRecHitStudy_cfg.py geometry=D82
#
#   Options for geometry D88, D92, D93
#
###############################################################################
import FWCore.ParameterSet.Config as cms
import os, sys, imp, re
import FWCore.ParameterSet.VarParsing as VarParsing

####################################################################
### Run as
## $cmsRun RecoResults_cfg.py geometry=D88 iter=$i

### SETUP OPTIONS
options = VarParsing.VarParsing('standard')
options.register('geometry',
                 "D92",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "geometry of operations: D88, D92, D93")

### get and parse the command line arguments
# options.parseArguments()

# print(options)

# options = VarParsing.VarParsing('standard')
options.register('iter',
                 "1",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "iterations 0 1 2 ....")

### get and parse the command line arguments
options.parseArguments()

print(options)

####################################################################
# Use the options
from Configuration.Eras.Era_Phase2C11M9_cff import Phase2C11M9
process = cms.Process('ReadHGCalReco',Phase2C11M9)

geomFile = "Configuration.Geometry.GeometryExtended2026"+ options.geometry +"Reco_cff"
fileInput = "root://se01.indiacms.res.in//cms/store/user/idas/SimOut/geomval/etaphi_debug_reeval/CMSSW_12_6_X_2022-09-27-2300/Extended2026" + options.geometry + "/step3_" + options.iter + ".root"
fileName = "hgcRecHit"+ options.geometry +"_"+ options.iter +".root"

print("Geometry file: ", geomFile)
print("Input file:    ", fileInput)
print("Output file:   ", fileName)

process.load(geomFile)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fileInput) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.load('Validation.HGCalValidation.hgcalRecHitStudy_cff')

process.hgcalRecHitStudyEE = cms.EDAnalyzer('ReadHGCalRecoResults',
                                            detectorName = cms.string('HGCalEESensitive'),
                                            source = cms.InputTag('HGCalRecHit', 'HGCEERecHits')
)

process.hgcalRecHitStudyFH = process.hgcalRecHitStudyEE.clone(detectorName  = cms.string("HGCalHESiliconSensitive"),
                                                              source        = cms.InputTag("HGCalRecHit", "HGCHEFRecHits")
                                                          )

process.hgcalRecHitStudyBH = process.hgcalRecHitStudyEE.clone( detectorName  = cms.string("HGCalHEScintillatorSensitive"),
                                                       source        = cms.InputTag("HGCalRecHit", "HGCHEBRecHits")
                                                   )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(fileName),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

#SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.p = cms.Path(process.hgcalRecHitStudyEE+process.hgcalRecHitStudyFH+process.hgcalRecHitStudyBH)

