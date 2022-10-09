// -*- C++ -*-
//
// Package:    Demo/ReadHGCalRecoResults
//
/**\class ReadHGCalRecoResults
// Derived from : HGCalRecHitStudy.cc

 Description: To read the recohits of HGCAL

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Indranil Das
//         Created:  Thu, 06 Oct 2022 06:18:11 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <fstream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include "CoralBase/Exception.h"

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToModule.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TGraph.h>

//
// class declaration
//

using reco::TrackCollection;

class ReadHGCalRecoResults : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit ReadHGCalRecoResults(const edm::ParameterSet&);
  ~ReadHGCalRecoResults();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  
  // ----------member data ---------------------------
  // edm::EDGetTokenT<edm::SimTrackContainer> tSimTrackContainer; 
  // edm::EDGetTokenT<edm::PCaloHitContainer> tSimCaloHitContainer; 
  // std::string name;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_; 
  hgcal::RecHitTools rhtools_;
  const std::string nameDetector_;
  const edm::EDGetTokenT<HGCRecHitCollection> recHitSource_;
  const edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> tok_hgcaldd_;
  const edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> tok_hgcGeom_;
  
  //TH1D *hPt;
  //TH2D *hXYhits;
  // For rechittool z positions. The 0 and 1 are for -ve and +ve, respectively.
  TH1D *hEF;
  TH1D *hECN;
  TH1D *hECK;
  TH1D *hESc;

  TH1D **hELossLayer0;
  TH1D **hELossLayer1;
  TH1D **hNCellsLayer0;
  TH1D **hNCellsLayer1;

  TH2D **hXYhitsF0;
  TH2D **hXYhitsCN0;
  TH2D **hXYhitsCK0;
  TH2D **hXYhitsB0;

  TH2D **hXYhitsF1;
  TH2D **hXYhitsCN1;
  TH2D **hXYhitsCK1;
  TH2D **hXYhitsB1;

  TH2D **hEPhitsF0;
  TH2D **hEPhitsCN0;
  TH2D **hEPhitsCK0;
  TH2D **hEPhitsB0;

  TH2D **hEPhitsF1;
  TH2D **hEPhitsCN1;
  TH2D **hEPhitsCK1;
  TH2D **hEPhitsB1;

  TH2D **hXYFailhitsF0;
  TH2D **hXYFailhitsCN0;
  TH2D **hXYFailhitsCK0;
  TH2D **hXYFailhitsB0;

  TH2D **hXYFailhitsF1;
  TH2D **hXYFailhitsCN1;
  TH2D **hXYFailhitsCK1;
  TH2D **hXYFailhitsB1;

  TH2D **hEPFailhitsF0;
  TH2D **hEPFailhitsCN0;
  TH2D **hEPFailhitsCK0;
  TH2D **hEPFailhitsB0;

  TH2D **hEPFailhitsF1;
  TH2D **hEPFailhitsCN1;
  TH2D **hEPFailhitsCK1;
  TH2D **hEPFailhitsB1;

  TH1D **hELossLayerF0;
  TH1D **hELossLayerCN0;
  TH1D **hELossLayerCK0;
  TH1D **hELossLayerB0;

  TH1D **hELossLayerF1;
  TH1D **hELossLayerCN1;
  TH1D **hELossLayerCK1;
  TH1D **hELossLayerB1;

  TGraph **grXYhitsF0;
  TGraph **grXYhitsCN0;
  TGraph **grXYhitsCK0;
  TGraph **grXYhitsAR0;
  TGraph **grXYhitsB0;
  int ixyF0[50], ixyCN0[50], ixyCK0[50], ixyAR0[50], ixyB0[50];

  TGraph **grXYhitsF1;
  TGraph **grXYhitsCN1;
  TGraph **grXYhitsCK1;
  TGraph **grXYhitsAR1;
  TGraph **grXYhitsB1;
  int ixyF1[50], ixyCN1[50], ixyCK1[50], ixyAR1[50], ixyB1[50];
  /////////////////////////////////

  // For rechittool z positions. The 0 and 1 are for -ve and +ve, respectively.
  TGraph **grEtaPhihitsF0;
  TGraph **grEtaPhihitsCN0;
  TGraph **grEtaPhihitsCK0;
  TGraph **grEtaPhihitsB0;
  int iepF0[50], iepCN0[50], iepCK0[50], iepB0[50];

  TGraph **grEtaPhihitsF1;
  TGraph **grEtaPhihitsCN1;
  TGraph **grEtaPhihitsCK1;
  TGraph **grEtaPhihitsB1;
  int iepF1[50], iepCN1[50], iepCK1[50], iepB1[50];

};

//
// constructors and destructor
//
ReadHGCalRecoResults::ReadHGCalRecoResults(const edm::ParameterSet& iConfig)
  :
  // tSimTrackContainer(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simtrack"))),
  // tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simhits")))
  nameDetector_(iConfig.getParameter<std::string>("detectorName")),
  recHitSource_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("source"))),
  tok_hgcaldd_(esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(
											     edm::ESInputTag{"", nameDetector_})),
  tok_hgcGeom_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", nameDetector_}))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  //hPt = fs->make<TH1D>("hPt" , "hPt" , 1000 , 0. , 1000. );
  //hXYhits = fs->make<TH2D>("hXYhits","Hits in XY", 600, -300., 300., 600, -300., 300.);

  hEF = fs->make<TH1D>("hEF" , "hEF" , 1000 , 0. , 50. );
  hECN = fs->make<TH1D>("hECN" , "hECN" , 1000 , 0. , 50. );
  hECK = fs->make<TH1D>("hECK" , "hECK" , 1000 , 0. , 50. );
  hESc = fs->make<TH1D>("hESc" , "hESc" , 1000 , 0. , 50. );

  hELossLayer0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayer1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hNCellsLayer0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hNCellsLayer1 = new TH1D*[50]; // for 50 layers in earch +z-direction

  hXYhitsF0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYhitsCN0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYhitsCK0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYhitsB0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYhitsF1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYhitsCN1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYhitsCK1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYhitsB1 =  new TH2D*[50]; // for 50 layers in earch + z-direction

  hEPhitsF0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPhitsCN0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPhitsCK0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPhitsB0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPhitsF1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPhitsCN1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPhitsCK1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPhitsB1 =  new TH2D*[50]; // for 50 layers in earch + z-direction

  hXYFailhitsF0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYFailhitsCN0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYFailhitsCK0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYFailhitsB0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hXYFailhitsF1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYFailhitsCN1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYFailhitsCK1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hXYFailhitsB1 =  new TH2D*[50]; // for 50 layers in earch + z-direction

  hEPFailhitsF0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPFailhitsCN0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPFailhitsCK0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPFailhitsB0 =  new TH2D*[50]; // for 50 layers in earch - z-direction
  hEPFailhitsF1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPFailhitsCN1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPFailhitsCK1 =  new TH2D*[50]; // for 50 layers in earch + z-direction
  hEPFailhitsB1 =  new TH2D*[50]; // for 50 layers in earch + z-direction

  hELossLayerF0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayerCN0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayerCK0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayerB0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayerF1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hELossLayerCN1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hELossLayerCK1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hELossLayerB1 = new TH1D*[50]; // for 50 layers in earch +z-direction

  grXYhitsF0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsCN0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsCK0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsAR0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsB0 =  new TGraph*[50];
  grXYhitsF1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsCN1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsCK1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsAR1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grXYhitsB1 =  new TGraph*[50];

  grEtaPhihitsF0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsCN0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsCK0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsB0 =  new TGraph*[50];
  grEtaPhihitsF1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsCN1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsCK1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  grEtaPhihitsB1 =  new TGraph*[50];

  for(int i=1;i<=50;i++)
    hELossLayer0[i] = fs->make<TH1D>(Form("hELossLayer0_%02d",i),Form("Rec:hELossLayer0_%02d",i), 500, 0., 5000.);
  for(int i=1;i<=50;i++)
    hELossLayer1[i] = fs->make<TH1D>(Form("hELossLayer1_%02d",i),Form("Rec:hELossLayer1_%02d",i), 500, 0., 5000.);
  for(int i=1;i<=50;i++)
    hNCellsLayer0[i] = fs->make<TH1D>(Form("hNCellsLayer0_%02d",i),Form("Rec:hNCellsLayer0_%02d",i), 200, -0.5, 199.5);
  for(int i=1;i<=50;i++)
    hNCellsLayer1[i] = fs->make<TH1D>(Form("hNCellsLayer1_%02d",i),Form("Rec:hNCellsLayer1_%02d",i), 200, -0.5, 199.5);
  for(int i=1;i<=50;i++)
    hXYhitsF0[i] = fs->make<TH2D>(Form("hXYhitsF0_layer_%02d",i),Form("Rec:HitsF0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCN0[i] = fs->make<TH2D>(Form("hXYhitsCN0_layer_%02d",i),Form("Rec:HitsCN0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCK0[i] = fs->make<TH2D>(Form("hXYhitsCK0_layer_%02d",i),Form("Rec:HitsCK0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsB0[i] = fs->make<TH2D>(Form("hXYhitsB0_layer_%02d",i),Form("Rec:HitsB0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsF1[i] = fs->make<TH2D>(Form("hXYhitsF1_layer_%02d",i),Form("Rec:HitsF1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCN1[i] = fs->make<TH2D>(Form("hXYhitsCN1_layer_%02d",i),Form("Rec:HitsCN1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCK1[i] = fs->make<TH2D>(Form("hXYhitsCK1_layer_%02d",i),Form("Rec:HitsCK1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsB1[i] = fs->make<TH2D>(Form("hXYhitsB1_layer_%02d",i),Form("Rec:HitsB1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);

  for(int i=1;i<=50;i++)
    hEPhitsF0[i] = fs->make<TH2D>(Form("hEPhitsF0_layer_%02d",i),Form("Rec:HitsF0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCN0[i] = fs->make<TH2D>(Form("hEPhitsCN0_layer_%02d",i),Form("Rec:HitsCN0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCK0[i] = fs->make<TH2D>(Form("hEPhitsCK0_layer_%02d",i),Form("Rec:HitsCK0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsB0[i] = fs->make<TH2D>(Form("hEPhitsB0_layer_%02d",i),Form("Rec:HitsB0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsF1[i] = fs->make<TH2D>(Form("hEPhitsF1_layer_%02d",i),Form("Rec:HitsF1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCN1[i] = fs->make<TH2D>(Form("hEPhitsCN1_layer_%02d",i),Form("Rec:HitsCN1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCK1[i] = fs->make<TH2D>(Form("hEPhitsCK1_layer_%02d",i),Form("Rec:HitsCK1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsB1[i] = fs->make<TH2D>(Form("hEPhitsB1_layer_%02d",i),Form("Rec:HitsB1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);

  for(int i=1;i<=50;i++)
    hXYFailhitsF0[i] = fs->make<TH2D>(Form("hXYFailhitsF0_layer_%02d",i),Form("Rec:FailhitsF0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCN0[i] = fs->make<TH2D>(Form("hXYFailhitsCN0_layer_%02d",i),Form("Rec:FailhitsCN0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCK0[i] = fs->make<TH2D>(Form("hXYFailhitsCK0_layer_%02d",i),Form("Rec:FailhitsCK0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsB0[i] = fs->make<TH2D>(Form("hXYFailhitsB0_layer_%02d",i),Form("Rec:FailhitsB0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsF1[i] = fs->make<TH2D>(Form("hXYFailhitsF1_layer_%02d",i),Form("Rec:FailhitsF1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCN1[i] = fs->make<TH2D>(Form("hXYFailhitsCN1_layer_%02d",i),Form("Rec:FailhitsCN1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCK1[i] = fs->make<TH2D>(Form("hXYFailhitsCK1_layer_%02d",i),Form("Rec:FailhitsCK1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsB1[i] = fs->make<TH2D>(Form("hXYFailhitsB1_layer_%02d",i),Form("Rec:FailhitsB1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  
  for(int i=1;i<=50;i++)
    hEPFailhitsF0[i] = fs->make<TH2D>(Form("hEPFailhitsF0_layer_%02d",i),Form("Rec:FailhitsF0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCN0[i] = fs->make<TH2D>(Form("hEPFailhitsCN0_layer_%02d",i),Form("Rec:FailhitsCN0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCK0[i] = fs->make<TH2D>(Form("hEPFailhitsCK0_layer_%02d",i),Form("Rec:FailhitsCK0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsB0[i] = fs->make<TH2D>(Form("hEPFailhitsB0_layer_%02d",i),Form("Rec:FailhitsB0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsF1[i] = fs->make<TH2D>(Form("hEPFailhitsF1_layer_%02d",i),Form("Rec:FailhitsF1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCN1[i] = fs->make<TH2D>(Form("hEPFailhitsCN1_layer_%02d",i),Form("Rec:FailhitsCN1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCK1[i] = fs->make<TH2D>(Form("hEPFailhitsCK1_layer_%02d",i),Form("Rec:FailhitsCK1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsB1[i] = fs->make<TH2D>(Form("hEPFailhitsB1_layer_%02d",i),Form("Rec:FailhitsB1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);

  for(int i=1;i<=50;i++)
    hELossLayerF0[i] = fs->make<TH1D>(Form("hELossLayerF0_%02d",i),Form("Rec:hELossLayerF0_%02d",i), 50, 0., 50.);
  for(int i=1;i<=50;i++)
    hELossLayerCN0[i] = fs->make<TH1D>(Form("hELossLayerCN0_%02d",i),Form("Rec:hELossLayerCN0_%02d",i), 50, 0., 50.);
  for(int i=1;i<=50;i++)
    hELossLayerCK0[i] = fs->make<TH1D>(Form("hELossLayerCK0_%02d",i),Form("Rec:hELossLayerCK0_%02d",i), 100, 0., 100.);
  for(int i=1;i<=50;i++)
    hELossLayerB0[i] = fs->make<TH1D>(Form("hELossLayerB0_%02d",i),Form("Rec:hELossLayerB0_%02d",i), 100, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerF1[i] = fs->make<TH1D>(Form("hELossLayerF1_%02d",i),Form("Rec:hELossLayerF1_%02d",i), 50, 0., 50.);
  for(int i=1;i<=50;i++)
    hELossLayerCN1[i] = fs->make<TH1D>(Form("hELossLayerCN1_%02d",i),Form("Rec:hELossLayerCN1_%02d",i), 50, 0., 50.);
  for(int i=1;i<=50;i++)
    hELossLayerCK1[i] = fs->make<TH1D>(Form("hELossLayerCK1_%02d",i),Form("Rec:hELossLayerCK1_%02d",i), 100, 0., 100.);
  for(int i=1;i<=50;i++)
    hELossLayerB1[i] = fs->make<TH1D>(Form("hELossLayerB1_%02d",i),Form("Rec:hELossLayerB1_%02d",i), 100, 0., 500.);

  for(int i=1;i<=50;i++){
    grXYhitsF0[i] = fs->make<TGraph>(0);
    grXYhitsF0[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("Rec:HitsF0 in XY for layer %d",i));
    grXYhitsCN0[i] = fs->make<TGraph>(0);
    grXYhitsCN0[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("Rec:HitsCN0 in XY for layer %d",i));
    grXYhitsCK0[i] = fs->make<TGraph>(0);
    grXYhitsCK0[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("Rec:HitsCK0 in XY for layer %d",i));
    grXYhitsAR0[i] = fs->make<TGraph>(0);
    grXYhitsAR0[i]->SetNameTitle(Form("grXYhitsAR0_layer_%02d",i),Form("Rec:HitsAr0 in XY for layer %d",i));
    grXYhitsB0[i] = fs->make<TGraph>(0);
    grXYhitsB0[i]->SetNameTitle(Form("grXYhitsB0_layer_%02d",i),Form("Rec:HitsB0 in XY for layer %d",i));
    ixyF0[i-1] = 0; ixyCN0[i-1] = 0; ixyCK0[i-1] = 0; ixyAR0[i-1] = 0; ixyB0[i-1] = 0; 

    grXYhitsF1[i] = fs->make<TGraph>(0);
    grXYhitsF1[i]->SetNameTitle(Form("grXYhitsF1_layer_%02d",i),Form("Rec:HitsF1 in XY for layer %d",i));
    grXYhitsCN1[i] = fs->make<TGraph>(0);
    grXYhitsCN1[i]->SetNameTitle(Form("grXYhitsCN1_layer_%02d",i),Form("Rec:HitsCN1 in XY for layer %d",i));
    grXYhitsCK1[i] = fs->make<TGraph>(0);
    grXYhitsCK1[i]->SetNameTitle(Form("grXYhitsCK1_layer_%02d",i),Form("Rec:HitsCK1 in XY for layer %d",i));
    grXYhitsAR1[i] = fs->make<TGraph>(0);
    grXYhitsAR1[i]->SetNameTitle(Form("grXYhitsAR1_layer_%02d",i),Form("Rec:HitsAr1 in XY for layer %d",i));
    grXYhitsB1[i] = fs->make<TGraph>(0);
    grXYhitsB1[i]->SetNameTitle(Form("grXYhitsB1_layer_%02d",i),Form("Rec:HitsB1 in XY for layer %d",i));
    ixyF1[i-1] = 0; ixyCN1[i-1] = 0; ixyCK1[i-1] = 0; ixyAR1[i-1] = 0;  ixyB1[i-1] = 0; 

    grEtaPhihitsF0[i] = fs->make<TGraph>(0);
    grEtaPhihitsF0[i]->SetNameTitle(Form("grEtaPhihitsF0_layer_%02d",i),Form("Rec:HitsF0 in XY for layer %d",i));
    grEtaPhihitsCN0[i] = fs->make<TGraph>(0);
    grEtaPhihitsCN0[i]->SetNameTitle(Form("grEtaPhihitsCN0_layer_%02d",i),Form("Rec:HitsCN0 in XY for layer %d",i));
    grEtaPhihitsCK0[i] = fs->make<TGraph>(0);
    grEtaPhihitsCK0[i]->SetNameTitle(Form("grEtaPhihitsCK0_layer_%02d",i),Form("Rec:HitsCK0 in XY for layer %d",i));
    grEtaPhihitsB0[i] = fs->make<TGraph>(0);
    grEtaPhihitsB0[i]->SetNameTitle(Form("grEtaPhihitsB0_layer_%02d",i),Form("Rec:HitsB0 in XY for layer %d",i));
    iepF0[i-1] = 0; iepCN0[i-1] = 0; iepCK0[i-1] = 0; iepB0[i-1] = 0; 

    grEtaPhihitsF1[i] = fs->make<TGraph>(0);
    grEtaPhihitsF1[i]->SetNameTitle(Form("grEtaPhihitsF1_layer_%02d",i),Form("Rec:HitsF1 in XY for layer %d",i));
    grEtaPhihitsCN1[i] = fs->make<TGraph>(0);
    grEtaPhihitsCN1[i]->SetNameTitle(Form("grEtaPhihitsCN1_layer_%02d",i),Form("Rec:HitsCN1 in XY for layer %d",i));
    grEtaPhihitsCK1[i] = fs->make<TGraph>(0);
    grEtaPhihitsCK1[i]->SetNameTitle(Form("grEtaPhihitsCK1_layer_%02d",i),Form("Rec:HitsCK1 in XY for layer %d",i));
    grEtaPhihitsB1[i] = fs->make<TGraph>(0);
    grEtaPhihitsB1[i]->SetNameTitle(Form("grEtaPhihitsB1_layer_%02d",i),Form("Rec:HitsB1 in XY for layer %d",i));
    iepF1[i-1] = 0; iepCN1[i-1] = 0; iepCK1[i-1] = 0;  iepB1[i-1] = 0; 


  }
  
  // name = iConfig.getParameter<std::string>("Detector");
  caloGeomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();
  
}

ReadHGCalRecoResults::~ReadHGCalRecoResults()
{

}

//
// member functions
//

// ------------ method called for each event  ------------
void
ReadHGCalRecoResults::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  const CaloGeometry &geomCalo = iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(geomCalo);

  // Here you access the generated particle information.  
  //bool ok(true);
  // bool ifNose_(false);
  // unsigned int ntot(0), nused(0);
  int verbosity_ = 0;
  double ElossLayer0[51],ElossLayer1[51];
  int nCellsLayer0[51],nCellsLayer1[51];
  for(int ilayer=0;ilayer<=50;ilayer++){
    ElossLayer0[ilayer] = ElossLayer1[ilayer] = 0.0;
    nCellsLayer0[ilayer] = nCellsLayer1[ilayer] = 0;
  }
  
  const edm::ESHandle<HGCalGeometry>& geom = iSetup.getHandle(tok_hgcGeom_);
  if (!geom.isValid())
    edm::LogWarning("HGCalValidation") << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
  const HGCalGeometry* geom0 = geom.product();

  const edm::Handle<HGCRecHitCollection>& theRecHitContainers = iEvent.getHandle(recHitSource_);
  if (theRecHitContainers.isValid()) {
    if (verbosity_ > 0)
      edm::LogVerbatim("HGCalValidation") << nameDetector_ << " with " << theRecHitContainers->size() << " element(s)";

    for (const auto& it : *(theRecHitContainers.product())) {
      // ntot++;
      // nused++;
      DetId detId = it.id();
      // int layer = (ifNose_ ? HFNoseDetId(detId).layer()
      //                      : ((detId.det() == DetId::HGCalHSc) ? HGCScintillatorDetId(detId).layer()
      // 			      : HGCSiliconDetId(detId).layer()));
      //recHitValidation(detId, layer, geom0, &it);
      //GlobalPoint global = geom0->getPosition(detId);
      GlobalPoint global1 = rhtools_.getPosition(detId);
      double energy = it.energy()*1.e3;
      
      // float globalx = global.x();
      // float globaly = global.y();
      //float globalz = global.z();
      //h_RZ_->Fill(std::abs(globalz), global.perp());
      //hXYhits->Fill(global.x(),global.y());
      
      if (geom0->topology().valid(detId)) {

	if(rhtools_.isSilicon(detId)){
	  HGCSiliconDetId id(it.id());
	  int il = rhtools_.getLayerWithOffset(detId);
	  if(id.type()==HGCSiliconDetId::HGCalFine){
	    hEF->Fill(energy);

	    if(global1.z()<0.0){
	      grXYhitsF0[il]->SetPoint(ixyF0[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsF0[il]->SetPoint(iepF0[il-1]++,global1.eta(), global1.phi());
	      hXYhitsF0[il]->Fill(global1.x(),global1.y());
	      hEPhitsF0[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerF0[il]->Fill(energy);
	      ElossLayer0[il] += energy;
	      nCellsLayer0[il]++;
	    }else{
	      grXYhitsF1[il]->SetPoint(ixyF1[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsF1[il]->SetPoint(iepF1[il-1]++,global1.eta(), global1.phi());
	      hXYhitsF1[il]->Fill(global1.x(),global1.y());
	      hEPhitsF1[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerF1[il]->Fill(energy);
	      ElossLayer1[il] += energy;
	      nCellsLayer1[il]++;
	    }

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	    hECN->Fill(energy);

	    if(global1.z()<0.0){
	      grXYhitsCN0[il]->SetPoint(ixyCN0[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsCN0[il]->SetPoint(iepCN0[il-1]++,global1.eta(), global1.phi());
	      hXYhitsCN0[il]->Fill(global1.x(),global1.y());
	      hEPhitsCN0[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerCN0[il]->Fill(energy);
	      ElossLayer0[il] += energy;
	      nCellsLayer0[il]++;
	    }else{
	      grXYhitsCN1[il]->SetPoint(ixyCN1[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsCN1[il]->SetPoint(iepCN1[il-1]++,global1.eta(), global1.phi());
	      hXYhitsCN1[il]->Fill(global1.x(),global1.y());
	      hEPhitsCN1[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerCN1[il]->Fill(energy);
	      ElossLayer1[il] += energy;
	      nCellsLayer1[il]++;
	    }

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	    hECK->Fill(energy);
	    
	    if(global1.z()<0.0){
	      grXYhitsCK0[il]->SetPoint(ixyCK0[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsCK0[il]->SetPoint(iepCK0[il-1]++,global1.eta(), global1.phi());
	      hXYhitsCK0[il]->Fill(global1.x(),global1.y());
	      hEPhitsCK0[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerCK0[il]->Fill(energy);
	      ElossLayer0[il] += energy;
	      nCellsLayer0[il]++;
	    }else{
	      grXYhitsCK1[il]->SetPoint(ixyCK1[il-1]++,global1.x(),global1.y());
	      grEtaPhihitsCK1[il]->SetPoint(iepCK1[il-1]++,global1.eta(), global1.phi());
	      hXYhitsCK1[il]->Fill(global1.x(),global1.y());
	      hEPhitsCK1[il]->Fill(global1.eta(),global1.phi());
	      hELossLayerCK1[il]->Fill(energy);
	      ElossLayer1[il] += energy;
	      nCellsLayer1[il]++;
	    }
	    
	  }
	  //The following line by Pruthvi to number the cells with id0 and detId
	  if(rhtools_.getCell(detId).first + rhtools_.getCell(detId).second <= 2){
	    if(global1.z()<0.0)
	      grXYhitsAR0[il]->SetPoint(ixyAR0[il-1]++,global1.x(),global1.y());
	    else
	      grXYhitsAR1[il]->SetPoint(ixyAR1[il-1]++,global1.x(),global1.y());
	  }
	}else if(rhtools_.isScintillator(detId)){
	
	  //HGCScintillatorDetId id(itHit->id());
	  int il = rhtools_.getLayerWithOffset(detId);
	  hESc->Fill(energy);

	  if (global1.z() < 0.0){
	    grXYhitsB0[il]->SetPoint(ixyB0[il]++, global1.x(), global1.y());
	    grEtaPhihitsB0[il]->SetPoint(iepB0[il]++, global1.eta(), global1.phi());
	    hXYhitsB0[il]->Fill(global1.x(),global1.y());
	    hEPhitsB0[il]->Fill(global1.eta(),global1.phi());
	    hELossLayerB0[il]->Fill(energy);
	    ElossLayer0[il] += energy;
	    nCellsLayer0[il]++;
	  }else{
	    grXYhitsB1[il]->SetPoint(ixyB1[il]++, global1.x(), global1.y());
	    grEtaPhihitsB1[il]->SetPoint(iepB1[il]++, global1.eta(), global1.phi());
	    hXYhitsB1[il]->Fill(global1.x(),global1.y());
	    hEPhitsB1[il]->Fill(global1.eta(),global1.phi());
	    hELossLayerB1[il]->Fill(energy);
	    ElossLayer1[il] += energy;
	    nCellsLayer1[il]++;
	  }
	  
	}//Silicon or scintillator
      
	///.................
      }else{//valid topology else invalid

	if(rhtools_.isSilicon(detId)){
	  HGCSiliconDetId id(it.id());
	  int il = rhtools_.getLayerWithOffset(detId);
	  if(id.type()==HGCSiliconDetId::HGCalFine){

	    if(global1.z()<0.0){
	      hXYFailhitsF0[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsF0[il]->Fill(global1.eta(),global1.phi());
	    }else{
	      hXYFailhitsF1[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsF1[il]->Fill(global1.eta(),global1.phi());
	    }

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThin){

	    if(global1.z()<0.0){
	      hXYFailhitsCN0[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsCN0[il]->Fill(global1.eta(),global1.phi());
	    }else{
	      hXYFailhitsCN1[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsCN1[il]->Fill(global1.eta(),global1.phi());
	    }

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	    
	    if(global1.z()<0.0){
	      hXYFailhitsCK0[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsCK0[il]->Fill(global1.eta(),global1.phi());
	    }else{
	      hXYFailhitsCK1[il]->Fill(global1.x(),global1.y());
	      hEPFailhitsCK1[il]->Fill(global1.eta(),global1.phi());
	    }
	    
	  }
	}else if(rhtools_.isScintillator(detId)){
	
	  int il = rhtools_.getLayerWithOffset(detId);
	  if (global1.z() < 0.0){
	    hXYFailhitsB0[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsB0[il]->Fill(global1.eta(),global1.phi());
	  }else{
	    hXYFailhitsB1[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsB1[il]->Fill(global1.eta(),global1.phi());
	  }
	  
	}//Silicon or scintillator

      }//invalid topology
    }//loop over iterator
  }//is Valid container

  for(int i=1;i<=50;i++){
    if(ElossLayer0[i]>0.0)
      hELossLayer0[i]->Fill(ElossLayer0[i]);
    if(nCellsLayer0[i]>0)
      hNCellsLayer0[i]->Fill(nCellsLayer0[i]);
 
    if(ElossLayer1[i]>0.0)
      hELossLayer1[i]->Fill(ElossLayer1[i]);
    if(nCellsLayer1[i]>0)
      hNCellsLayer1[i]->Fill(nCellsLayer1[i]);
  }

}



// ------------ method called once each job just before starting event loop  ------------
void
ReadHGCalRecoResults::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ReadHGCalRecoResults::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ReadHGCalRecoResults::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ReadHGCalRecoResults);
