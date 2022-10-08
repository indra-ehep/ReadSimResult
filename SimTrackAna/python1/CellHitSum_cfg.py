####################################################################
import FWCore.ParameterSet.Config as cms
import os, sys, imp, re
import FWCore.ParameterSet.VarParsing as VarParsing


####################################################################
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

geomFile = "Configuration.Geometry.GeometryExtended2026"+ options.geometry +"Reco_cff"
#fileInput = "root://se01.indiacms.res.in//cms/store/user/idas/SimOut/geomval/etaphi_debug_reeval/CMSSW_12_6_X_2022-09-27-2300/Extended2026" + options.geometry + "/step1_" + options.iter + ".root"
#fileInput = "file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/etaphi_debug_reeval/CMSSW_12_5_0_pre5/Extended2026" + options.geometry + "/step1_" + options.iter + ".root"
fileInput = "file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/etaphi_debug_reeval/CMSSW_12_4_0_pre4/Extended2026" + options.geometry + "/step1_" + options.iter + ".root"
fileName = "geantoutput"+ options.geometry +"_"+ options.iter +".root"

print("Geometry file: ", geomFile)
print("Input file:    ", fileInput)
print("Output file:   ", fileName)

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load(geomFile)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(fileInput)

)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.prodEE = cms.EDAnalyzer('CellHitSum',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             Detector   = cms.string("HGCalEESensitive"),
                         )


process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
)

process.prodHEB = process.prodHEF.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback"),
    Detector   = cms.string("HGCalHEScintillatorSensitive"),
)

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(fileName)
 )

process.p = cms.Path(process.prodEE*process.prodHEF*process.prodHEB)
