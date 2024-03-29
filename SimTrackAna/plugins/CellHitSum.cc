// -*- C++ -*-
//
// Package:    Demo/CellHitSum
// Class:      CellHitSum
//
/**\class CellHitSum TrackAnalyzer.cc Track/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Indranil Das
//         Created:  Wed, 25 Aug 2021 06:18:11 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <fstream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TMath.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class CellHitSum : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  //Implemented following Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
  struct energysum {
    energysum() {
      etotal = 0;
      for (int i = 0; i < 6; ++i)
        eTime[i] = 0.;
    }
    double eTime[6], etotal;
  };
  
  struct waferinfo {
    waferinfo() {      
      layer = u = v = type = -999;
    }
    int layer, u, v, type;
  };
  
  struct hitsinfo {
    hitsinfo() {
      x = y = z = phi = eta = trkpt = trketa = trkphi = 0.0;
      cell = cell2 = sector = sector2 = type = layer = pdg = charge = 0;
      hitid = nhits = 0;
      isMu = false;
    }
    double x, y, z, phi, eta, trkpt, trketa, trkphi;
    int cell, cell2, sector, sector2, type, layer, pdg, charge;
    unsigned int hitid, nhits;
    bool isMu ;
  };
  
  
  explicit CellHitSum(const edm::ParameterSet&);
  ~CellHitSum();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<edm::SimTrackContainer> tSimTrackContainer; 
  edm::EDGetTokenT<edm::PCaloHitContainer> tSimCaloHitContainer; 
  std::string name;
  edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomToken_;
  
  TH1D *hCharge;
  TH1D *hChargeLowELoss;

  TH1D *hPt; 
  TH1D *hPtNoGen;
  TH1D *hPtLowELoss;

  TH1D *hEta; 
  TH1D *hEtaCell;
  TH1D *hEtaLowELoss;

  TH1D *hPhi; 
  TH1D *hPhiCell;
  TH1D *hPhiLowELoss;

  TH1D *hPDG;
  TH1D *hPDGLowELoss;
  
  TH1D *hELossEE;
  TH1D *hELossEEF ;
  TH1D *hELossEECN ;
  TH1D *hELossEECK ;
  TH1D *hELossHEF ;
  TH1D *hELossHEFF ;
  TH1D *hELossHEFCN ;
  TH1D *hELossHEFCK ;
  TH1D *hELossHEB ;

  TH1D *hELossCSinBunchEE;
  TH1D *hELossCSinBunchEEF ;
  TH1D *hELossCSinBunchEECN ;
  TH1D *hELossCSinBunchEECK ;
  TH1D *hELossCSinBunchHEF ;
  TH1D *hELossCSinBunchHEFF ;
  TH1D *hELossCSinBunchHEFCN ;
  TH1D *hELossCSinBunchHEFCK ;
  TH1D *hELossCSinBunchHEFCNFiltered ;
  TH1D *hELossCSinBunchHEFCNNoise ;
  
  TH1D *hELossCSmissedEE;
  TH1D *hELossCSmissedEEF ;
  TH1D *hELossCSmissedEECN ;
  TH1D *hELossCSmissedEECK ;
  TH1D *hELossCSmissedHEF ;
  TH1D *hELossCSmissedHEFF ;
  TH1D *hELossCSmissedHEFCN ;
  TH1D *hELossCSmissedHEFCK ;

  TH1D *hELossCSMaxEE;
  TH1D *hELossCSMaxEEF ;
  TH1D *hELossCSMaxEECN ;
  TH1D *hELossCSMaxEECK ;
  TH1D *hELossCSMaxHEF ;
  TH1D *hELossCSMaxHEFF ;
  TH1D *hELossCSMaxHEFCN ;
  TH1D *hELossCSMaxHEFCK ;
  
  TH1D *hHxELossCSMaxF ;
  TH1D *hHxELossCSMaxCN ;
  TH1D *hHxELossCSMaxCK ;
  TH1D *hNHxELossCSMaxF ;
  TH1D *hNHxELossCSMaxCN ;
  TH1D *hNHxELossCSMaxCK ;
  
  TH1D **hELossDQMEqV ;
  TH1D **hELossLayer0;
  TH1D **hELossLayer1;

  TH1D **hNCellsLayer0;
  TH1D **hNCellsLayer1;
  
  // TH2D *hYZhits;
  TH2D **hXYhits;
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

  TH2D **hXYhitsLELCN;
  TH2D **hXYhitsHELCN;
  TH2D **hXYhitsLELCK;
  TH2D **hXYhitsHELCK;
  TH2D **hNHxXYhitsF;
  TH2D **hNHxXYhitsCN;
  TH2D **hNHxXYhitsCK;

  // For rechittool z positions. The 0 and 1 are for -ve and +ve, respectively.
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

  //////////////////////////////////////////

  TH1D **hELCSMaxF ;
  TH1D **hELCSMaxCN ;
  TH1D **hELCSMaxCK ;

  TH1D **hHxELCSMaxF ;
  TH1D **hHxELCSMaxCN ;
  TH1D **hHxELCSMaxCK ;
  TH1D **hNHxELCSMaxF ;
  TH1D **hNHxELCSMaxCN ;
  TH1D **hNHxELCSMaxCK ;

  TH2D *hXYLowELosshitsF;
  TH2D *hXYLowELosshitsCN;
  TH2D *hXYLowELosshitsCK;
  TH2D *hXYmissedhits;
  TH2D *hYZLowELosshitsF;
  TH2D *hYZLowELosshitsCN;
  TH2D *hYZLowELosshitsCK;
  TH2D *hYZLLowELosshitsHEFCN;
  TH2D *hYZmissedhits;

  TH1D *hXLowELosshitsHEFCN;
  TH1D *hYLowELosshitsHEFCN;
  TH1D *hZLowELosshitsHEFCN;
  
  TH2D *hYZhitsEE;
  TH2D *hYZhitsHEF;
  TH2D *hYZhitsHEB;

  TH2D *hYZhitsEEF;
  TH2D *hYZhitsEECN;
  TH2D *hYZhitsEECK;

  TH2D *hYZhitsHEFF;
  TH2D *hYZhitsHEFCN;
  TH2D *hYZhitsHEFCK;
  
  TH2D *hRHTXYhits;
  TH2D *hRHTYZhitsEE;
  TH2D *hRHTYZhitsHEF;
  TH2D *hRHTYZhitsHEB;
  TH2D *hRHTYZhitsEEF;
  TH2D *hRHTYZhitsEECN;
  TH2D *hRHTYZhitsEECK;
  TH2D *hRHTYZhitsHEFF;
  TH2D *hRHTYZhitsHEFCN;
  TH2D *hRHTYZhitsHEFCK;

  TH2D *hRHTRZhitsEE;
  TH2D *hRHTRZhitsHEF;
  TH2D *hRHTRZhitsHEB;
  TH2D *hRHTRZhitsEEF;
  TH2D *hRHTRZhitsEECN;
  TH2D *hRHTRZhitsEECK;
  TH2D *hRHTRZhitsHEFF;
  TH2D *hRHTRZhitsHEFCN;
  TH2D *hRHTRZhitsHEFCK;

  TH2D *hRHTGlbRZhitsF ;
  TH2D *hRHTGlbRZhitsCN ;
  TH2D *hRHTGlbRZhitsCK ;
  TH2D *hRHTGlbRZhitsSci ;

  TH1D *hDiffX ;
  TH1D *hDiffY ;
  TH1D *hDiffZ ;

  TH1D *hCellThickness ;  
  
  std::vector<waferinfo> winfo;
  int evt = 0;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_; 
  hgcal::RecHitTools rhtools_;
  //edm::ConsumesCollector iC;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CellHitSum::CellHitSum(const edm::ParameterSet& iConfig)
  :
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("simhits"))),
  tSimTrackContainer(consumes<edm::SimTrackContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simtrack"))),
  tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iConfig.getUntrackedParameter<edm::InputTag>("simhits")))
  //name(iConfig.getParameter<std::string>("Detector")),
  //geomToken_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", name}))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  hCharge = fs->make<TH1D>("charge" , "Charges" , 200 , -20 , 20 );
  hChargeLowELoss = fs->make<TH1D>("charge LowELoss" , "Charges LowELoss" , 200 , -20 , 20 );

  hPDG = fs->make<TH1D>("hPDG" , "hPDG" , 10000 , -5000 , 5000 );
  hPDGLowELoss = fs->make<TH1D>("hPDGLowELoss" , "hPDGLowELoss" ,10000,-5000,5000);  

  hPt = fs->make<TH1D>("hPt" , "hPt" , 1000 , 0. , 1000. );
  hPtNoGen = fs->make<TH1D>("hPtNoGen" , "hPtNoGen" , 1000 , 0. , 1000. );
  hPtLowELoss = fs->make<TH1D>("hPtLowELoss" , "hPtLowELoss" , 1000 , 0. , 1000. );

  hEta = fs->make<TH1D>("hEta" , "hEta" , 100 , -5. , 5. );
  hEtaCell = fs->make<TH1D>("hEtaCell" , "hEtaCell" , 100 , -5. , 5.  );
  hEtaLowELoss = fs->make<TH1D>("hEtaLowELoss" , "hEtaLowELoss" , 100 , -5. , 5. );

  hPhi = fs->make<TH1D>("hPhi" , "hPhi" , 100 , -5. , 5.  );
  hPhiCell = fs->make<TH1D>("hPhiCell" , "hPhiCell" , 100 , -5. , 5.  );
  hPhiLowELoss = fs->make<TH1D>("hPhiLowELoss" , "hPhiLowELoss" , 100 , -5. , 5.  );

  hELossEE = fs->make<TH1D>("hELossEE","hELossEE", 1000, 0., 1000.);
  hELossEEF = fs->make<TH1D>("hELossEEF","hELossEEF", 1000, 0., 1000.);
  hELossEECN = fs->make<TH1D>("hELossEECN","hELossEECN", 1000, 0., 1000.);
  hELossEECK = fs->make<TH1D>("hELossEECK","hELossEECK", 1000, 0., 1000.);

  hELossHEF = fs->make<TH1D>("hELossHEF","hELossHEF", 1000, 0., 1000.);
  hELossHEFF = fs->make<TH1D>("hELossHEFF","hELossHEFF", 1000, 0., 1000.);
  hELossHEFCN = fs->make<TH1D>("hELossHEFCN","hELossHEFCN", 1000, 0., 1000.);
  hELossHEFCK = fs->make<TH1D>("hELossHEFCK","hELossHEFCK", 1000, 0., 1000.);
  
  hELossHEB = fs->make<TH1D>("hELossHEB","hELossHEB", 1000, 0., 1000.);

  hELossCSinBunchEE = fs->make<TH1D>("hELossCSinBunchEE","hELossCSinBunchEE", 1000, 0., 1000.);
  hELossCSinBunchEEF = fs->make<TH1D>("hELossCSinBunchEEF","hELossCSinBunchEEF", 1000, 0., 1000.);
  hELossCSinBunchEECN = fs->make<TH1D>("hELossCSinBunchEECN","hELossCSinBunchEECN", 1000, 0., 1000.);
  hELossCSinBunchEECK = fs->make<TH1D>("hELossCSinBunchEECK","hELossCSinBunchEECK", 1000, 0., 1000.);
  hELossCSinBunchHEF = fs->make<TH1D>("hELossCSinBunchHEF","hELossCSinBunchHEF", 1000, 0., 1000.);
  hELossCSinBunchHEFF = fs->make<TH1D>("hELossCSinBunchHEFF","hELossCSinBunchHEFF", 1000, 0., 1000.);
  hELossCSinBunchHEFCN = fs->make<TH1D>("hELossCSinBunchHEFCN","hELossCSinBunchHEFCN", 1000, 0., 1000.);
  hELossCSinBunchHEFCK = fs->make<TH1D>("hELossCSinBunchHEFCK","hELossCSinBunchHEFCK", 1000, 0., 1000.);
  hELossCSinBunchHEFCNFiltered = fs->make<TH1D>("hELossCSinBunchHEFCNFiltered","hELossCSinBunchHEFCNFiltered", 1000, 0., 1000.);
  hELossCSinBunchHEFCNNoise = fs->make<TH1D>("hELossCSinBunchHEFCNNoise","hELossCSinBunchHEFCNNoise", 1000, 0., 1000.);

  hELossCSmissedEE = fs->make<TH1D>("hELossCSmissedEE","hELossCSmissedEE", 1000, 0., 1000.);
  hELossCSmissedEEF = fs->make<TH1D>("hELossCSmissedEEF","hELossCSmissedEEF", 1000, 0., 1000.);
  hELossCSmissedEECN = fs->make<TH1D>("hELossCSmissedEECN","hELossCSmissedEECN", 1000, 0., 1000.);
  hELossCSmissedEECK = fs->make<TH1D>("hELossCSmissedEECK","hELossCSmissedEECK", 1000, 0., 1000.);
  hELossCSmissedHEF = fs->make<TH1D>("hELossCSmissedHEF","hELossCSmissedHEF", 1000, 0., 1000.);
  hELossCSmissedHEFF = fs->make<TH1D>("hELossCSmissedHEFF","hELossCSmissedHEFF", 1000, 0., 1000.);
  hELossCSmissedHEFCN = fs->make<TH1D>("hELossCSmissedHEFCN","hELossCSmissedHEFCN", 1000, 0., 1000.);
  hELossCSmissedHEFCK = fs->make<TH1D>("hELossCSmissedHEFCK","hELossCSmissedHEFCK", 1000, 0., 1000.);

  hELossCSMaxEE = fs->make<TH1D>("hELossCSMaxEE","hELossCSMaxEE", 1000, 0., 1000.);
  hELossCSMaxEEF = fs->make<TH1D>("hELossCSMaxEEF","hELossCSMaxEEF", 1000, 0., 1000.);
  hELossCSMaxEECN = fs->make<TH1D>("hELossCSMaxEECN","hELossCSMaxEECN", 1000, 0., 1000.);
  hELossCSMaxEECK = fs->make<TH1D>("hELossCSMaxEECK","hELossCSMaxEECK", 1000, 0., 1000.);
  hELossCSMaxHEF = fs->make<TH1D>("hELossCSMaxHEF","hELossCSMaxHEF", 1000, 0., 1000.);
  hELossCSMaxHEFF = fs->make<TH1D>("hELossCSMaxHEFF","hELossCSMaxHEFF", 1000, 0., 1000.);
  hELossCSMaxHEFCN = fs->make<TH1D>("hELossCSMaxHEFCN","hELossCSMaxHEFCN", 1000, 0., 1000.);
  hELossCSMaxHEFCK = fs->make<TH1D>("hELossCSMaxHEFCK","hELossCSMaxHEFCK", 1000, 0., 1000.);
  
  hHxELossCSMaxF = fs->make<TH1D>("hHxELossCSMaxF","hHxELossCSMaxF", 1000, 0., 1000.);
  hHxELossCSMaxCN = fs->make<TH1D>("hHxELossCSMaxCN","hHxELossCSMaxCN", 1000, 0., 1000.);
  hHxELossCSMaxCK = fs->make<TH1D>("hHxELossCSMaxCK","hHxELossCSMaxCK", 1000, 0., 1000.);
  
  hNHxELossCSMaxF = fs->make<TH1D>("hNHxELossCSMaxF","hNHxELossCSMaxF", 1000, 0., 1000.);
  hNHxELossCSMaxCN = fs->make<TH1D>("hNHxELossCSMaxCN","hNHxELossCSMaxCN", 1000, 0., 1000.);
  hNHxELossCSMaxCK = fs->make<TH1D>("hNHxELossCSMaxCK","hNHxELossCSMaxCK", 1000, 0., 1000.);


  hELCSMaxF =  new TH1D*[50]; 
  hELCSMaxCN =  new TH1D*[50]; 
  hELCSMaxCK =  new TH1D*[50]; 
  for(int i=1;i<=50;i++)
    hELCSMaxF[i] = fs->make<TH1D>(Form("hELCSMaxF_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSMaxCN[i] = fs->make<TH1D>(Form("hELCSMaxCN_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSMaxCK[i] = fs->make<TH1D>(Form("hELCSMaxCK_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  
  hHxELCSMaxF =  new TH1D*[50]; 
  hHxELCSMaxCN =  new TH1D*[50]; 
  hHxELCSMaxCK =  new TH1D*[50]; 
  for(int i=1;i<=50;i++)
    hHxELCSMaxF[i] = fs->make<TH1D>(Form("hHxELCSMaxF_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hHxELCSMaxCN[i] = fs->make<TH1D>(Form("hHxELCSMaxCN_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hHxELCSMaxCK[i] = fs->make<TH1D>(Form("hHxELCSMaxCK_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  
  hNHxELCSMaxF =  new TH1D*[50]; 
  hNHxELCSMaxCN =  new TH1D*[50]; 
  hNHxELCSMaxCK =  new TH1D*[50]; 
  for(int i=1;i<=50;i++)
    hNHxELCSMaxF[i] = fs->make<TH1D>(Form("hNHxELCSMaxF_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hNHxELCSMaxCN[i] = fs->make<TH1D>(Form("hNHxELCSMaxCN_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hNHxELCSMaxCK[i] = fs->make<TH1D>(Form("hNHxELCSMaxCK_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  
  hELossDQMEqV = new TH1D*[50]; // for 50 layers in earch +/- z-direction
  hELossLayer0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hELossLayer1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hNCellsLayer0 = new TH1D*[50]; // for 50 layers in earch -z-direction
  hNCellsLayer1 = new TH1D*[50]; // for 50 layers in earch +z-direction
  hXYhits =  new TH2D*[50]; // for 50 layers in earch +/- z-direction

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

  hXYhitsLELCN =  new TH2D*[50]; 
  hXYhitsHELCN =  new TH2D*[50]; 
  hXYhitsLELCK =  new TH2D*[50]; 
  hXYhitsHELCK =  new TH2D*[50]; 
  
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

  //gXYhitsF0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  //gXYhitsCN0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  //gXYhitsCK0 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  //gXYhitsF1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  //gXYhitsCN1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction
  //gXYhitsCK1 =  new TGraph*[50]; // for 50 layers in earch +/- z-direction


  hNHxXYhitsF =  new TH2D*[50]; 
  hNHxXYhitsCN =  new TH2D*[50]; 
  hNHxXYhitsCK =  new TH2D*[50]; 


  for(int i=1;i<=50;i++)
    hELossDQMEqV[i] = fs->make<TH1D>(Form("hELossDQMEqV_layer_%02d",i),Form("hELossDQMEqV_layer_%02d",i), 100, 0, 0.1);
  for(int i=1;i<=50;i++)
    hELossLayer0[i] = fs->make<TH1D>(Form("hELossLayer0_%02d",i),Form("Gen:hELossLayer0_%02d",i), 1000, 0., 1000.);
  for(int i=1;i<=50;i++)
    hELossLayer1[i] = fs->make<TH1D>(Form("hELossLayer1_%02d",i),Form("Gen:hELossLayer1_%02d",i), 1000, 0., 1000.);
  for(int i=1;i<=50;i++)
    hNCellsLayer0[i] = fs->make<TH1D>(Form("hNCellsLayer0_%02d",i),Form("Gen:hNCellsLayer0_%02d",i), 20, -0.5, 19.5);
  for(int i=1;i<=50;i++)
    hNCellsLayer1[i] = fs->make<TH1D>(Form("hNCellsLayer1_%02d",i),Form("Gen:hNCellsLayer1_%02d",i), 20, -0.5, 19.5);
  for(int i=1;i<=50;i++)
    hXYhits[i] = fs->make<TH2D>(Form("hXYhits_layer_%02d",i),Form("Gen:Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsF0[i] = fs->make<TH2D>(Form("hXYhitsF0_layer_%02d",i),Form("Gen:HitsF0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCN0[i] = fs->make<TH2D>(Form("hXYhitsCN0_layer_%02d",i),Form("Gen:HitsCN0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCK0[i] = fs->make<TH2D>(Form("hXYhitsCK0_layer_%02d",i),Form("Gen:HitsCK0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsB0[i] = fs->make<TH2D>(Form("hXYhitsB0_layer_%02d",i),Form("Gen:HitsB0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsF1[i] = fs->make<TH2D>(Form("hXYhitsF1_layer_%02d",i),Form("Gen:HitsF1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCN1[i] = fs->make<TH2D>(Form("hXYhitsCN1_layer_%02d",i),Form("Gen:HitsCN1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsCK1[i] = fs->make<TH2D>(Form("hXYhitsCK1_layer_%02d",i),Form("Gen:HitsCK1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsB1[i] = fs->make<TH2D>(Form("hXYhitsB1_layer_%02d",i),Form("Gen:HitsB1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);

  for(int i=1;i<=50;i++)
    hEPhitsF0[i] = fs->make<TH2D>(Form("hEPhitsF0_layer_%02d",i),Form("Gen:HitsF0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCN0[i] = fs->make<TH2D>(Form("hEPhitsCN0_layer_%02d",i),Form("Gen:HitsCN0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCK0[i] = fs->make<TH2D>(Form("hEPhitsCK0_layer_%02d",i),Form("Gen:HitsCK0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsB0[i] = fs->make<TH2D>(Form("hEPhitsB0_layer_%02d",i),Form("Gen:HitsB0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsF1[i] = fs->make<TH2D>(Form("hEPhitsF1_layer_%02d",i),Form("Gen:HitsF1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCN1[i] = fs->make<TH2D>(Form("hEPhitsCN1_layer_%02d",i),Form("Gen:HitsCN1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsCK1[i] = fs->make<TH2D>(Form("hEPhitsCK1_layer_%02d",i),Form("Gen:HitsCK1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPhitsB1[i] = fs->make<TH2D>(Form("hEPhitsB1_layer_%02d",i),Form("Gen:HitsB1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);

  for(int i=1;i<=50;i++)
    hXYFailhitsF0[i] = fs->make<TH2D>(Form("hXYFailhitsF0_layer_%02d",i),Form("Gen:FailhitsF0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCN0[i] = fs->make<TH2D>(Form("hXYFailhitsCN0_layer_%02d",i),Form("Gen:FailhitsCN0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCK0[i] = fs->make<TH2D>(Form("hXYFailhitsCK0_layer_%02d",i),Form("Gen:FailhitsCK0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsB0[i] = fs->make<TH2D>(Form("hXYFailhitsB0_layer_%02d",i),Form("Gen:FailhitsB0 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsF1[i] = fs->make<TH2D>(Form("hXYFailhitsF1_layer_%02d",i),Form("Gen:FailhitsF1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCN1[i] = fs->make<TH2D>(Form("hXYFailhitsCN1_layer_%02d",i),Form("Gen:FailhitsCN1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsCK1[i] = fs->make<TH2D>(Form("hXYFailhitsCK1_layer_%02d",i),Form("Gen:FailhitsCK1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYFailhitsB1[i] = fs->make<TH2D>(Form("hXYFailhitsB1_layer_%02d",i),Form("Gen:FailhitsB1 in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  
  for(int i=1;i<=50;i++)
    hEPFailhitsF0[i] = fs->make<TH2D>(Form("hEPFailhitsF0_layer_%02d",i),Form("Gen:FailhitsF0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCN0[i] = fs->make<TH2D>(Form("hEPFailhitsCN0_layer_%02d",i),Form("Gen:FailhitsCN0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCK0[i] = fs->make<TH2D>(Form("hEPFailhitsCK0_layer_%02d",i),Form("Gen:FailhitsCK0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsB0[i] = fs->make<TH2D>(Form("hEPFailhitsB0_layer_%02d",i),Form("Gen:FailhitsB0 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsF1[i] = fs->make<TH2D>(Form("hEPFailhitsF1_layer_%02d",i),Form("Gen:FailhitsF1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCN1[i] = fs->make<TH2D>(Form("hEPFailhitsCN1_layer_%02d",i),Form("Gen:FailhitsCN1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsCK1[i] = fs->make<TH2D>(Form("hEPFailhitsCK1_layer_%02d",i),Form("Gen:FailhitsCK1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  for(int i=1;i<=50;i++)
    hEPFailhitsB1[i] = fs->make<TH2D>(Form("hEPFailhitsB1_layer_%02d",i),Form("Gen:FailhitsB1 in EP for layer %d",i), 640, -3.2, 3.2, 640, -3.2, 3.2);
  
  for(int i=1;i<=50;i++)
    hELossLayerF0[i] = fs->make<TH1D>(Form("hELossLayerF0_%02d",i),Form("Gen:hELossLayerF0_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerCN0[i] = fs->make<TH1D>(Form("hELossLayerCN0_%02d",i),Form("Gen:hELossLayerCN0_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerCK0[i] = fs->make<TH1D>(Form("hELossLayerCK0_%02d",i),Form("Gen:hELossLayerCK0_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerB0[i] = fs->make<TH1D>(Form("hELossLayerB0_%02d",i),Form("Gen:hELossLayerB0_%02d",i), 1000, 0., 1000.);
  for(int i=1;i<=50;i++)
    hELossLayerF1[i] = fs->make<TH1D>(Form("hELossLayerF1_%02d",i),Form("Gen:hELossLayerF1_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerCN1[i] = fs->make<TH1D>(Form("hELossLayerCN1_%02d",i),Form("Gen:hELossLayerCN1_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerCK1[i] = fs->make<TH1D>(Form("hELossLayerCK1_%02d",i),Form("Gen:hELossLayerCK1_%02d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELossLayerB1[i] = fs->make<TH1D>(Form("hELossLayerB1_%02d",i),Form("Gen:hELossLayerB1_%02d",i), 1000, 0., 1000.);
  
  for(int i=1;i<=50;i++)
    hXYhitsLELCN[i] = fs->make<TH2D>(Form("hXYhitsLELCN_layer_%02d",i),Form("LELCN in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsHELCN[i] = fs->make<TH2D>(Form("hXYhitsHELCN_layer_%02d",i),Form("HELCN in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsLELCK[i] = fs->make<TH2D>(Form("hXYhitsLELCK_layer_%02d",i),Form("LELCK in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsHELCK[i] = fs->make<TH2D>(Form("hXYhitsHELCK_layer_%02d",i),Form("HELCK in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);

  for(int i=1;i<=50;i++){
    grXYhitsF0[i] = fs->make<TGraph>(0);
    grXYhitsF0[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grXYhitsCN0[i] = fs->make<TGraph>(0);
    grXYhitsCN0[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grXYhitsCK0[i] = fs->make<TGraph>(0);
    grXYhitsCK0[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    grXYhitsAR0[i] = fs->make<TGraph>(0);
    grXYhitsAR0[i]->SetNameTitle(Form("grXYhitsAR0_layer_%02d",i),Form("HitsAr0 in XY for layer %d",i));
    grXYhitsB0[i] = fs->make<TGraph>(0);
    grXYhitsB0[i]->SetNameTitle(Form("grXYhitsB0_layer_%02d",i),Form("HitsB0 in XY for layer %d",i));
    ixyF0[i-1] = 0; ixyCN0[i-1] = 0; ixyCK0[i-1] = 0; ixyAR0[i-1] = 0; ixyB0[i-1] = 0; 

    grXYhitsF1[i] = fs->make<TGraph>(0);
    grXYhitsF1[i]->SetNameTitle(Form("grXYhitsF1_layer_%02d",i),Form("HitsF1 in XY for layer %d",i));
    grXYhitsCN1[i] = fs->make<TGraph>(0);
    grXYhitsCN1[i]->SetNameTitle(Form("grXYhitsCN1_layer_%02d",i),Form("HitsCN1 in XY for layer %d",i));
    grXYhitsCK1[i] = fs->make<TGraph>(0);
    grXYhitsCK1[i]->SetNameTitle(Form("grXYhitsCK1_layer_%02d",i),Form("HitsCK1 in XY for layer %d",i));
    grXYhitsAR1[i] = fs->make<TGraph>(0);
    grXYhitsAR1[i]->SetNameTitle(Form("grXYhitsAR1_layer_%02d",i),Form("HitsAr1 in XY for layer %d",i));
    grXYhitsB1[i] = fs->make<TGraph>(0);
    grXYhitsB1[i]->SetNameTitle(Form("grXYhitsB1_layer_%02d",i),Form("HitsB1 in XY for layer %d",i));
    ixyF1[i-1] = 0; ixyCN1[i-1] = 0; ixyCK1[i-1] = 0; ixyAR1[i-1] = 0;  ixyB1[i-1] = 0; 

    grEtaPhihitsF0[i] = fs->make<TGraph>(0);
    grEtaPhihitsF0[i]->SetNameTitle(Form("grEtaPhihitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grEtaPhihitsCN0[i] = fs->make<TGraph>(0);
    grEtaPhihitsCN0[i]->SetNameTitle(Form("grEtaPhihitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grEtaPhihitsCK0[i] = fs->make<TGraph>(0);
    grEtaPhihitsCK0[i]->SetNameTitle(Form("grEtaPhihitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    grEtaPhihitsB0[i] = fs->make<TGraph>(0);
    grEtaPhihitsB0[i]->SetNameTitle(Form("grEtaPhihitsB0_layer_%02d",i),Form("HitsB0 in XY for layer %d",i));
    iepF0[i-1] = 0; iepCN0[i-1] = 0; iepCK0[i-1] = 0; iepB0[i-1] = 0; 

    grEtaPhihitsF1[i] = fs->make<TGraph>(0);
    grEtaPhihitsF1[i]->SetNameTitle(Form("grEtaPhihitsF1_layer_%02d",i),Form("HitsF1 in XY for layer %d",i));
    grEtaPhihitsCN1[i] = fs->make<TGraph>(0);
    grEtaPhihitsCN1[i]->SetNameTitle(Form("grEtaPhihitsCN1_layer_%02d",i),Form("HitsCN1 in XY for layer %d",i));
    grEtaPhihitsCK1[i] = fs->make<TGraph>(0);
    grEtaPhihitsCK1[i]->SetNameTitle(Form("grEtaPhihitsCK1_layer_%02d",i),Form("HitsCK1 in XY for layer %d",i));
    grEtaPhihitsB1[i] = fs->make<TGraph>(0);
    grEtaPhihitsB1[i]->SetNameTitle(Form("grEtaPhihitsB1_layer_%02d",i),Form("HitsB1 in XY for layer %d",i));
    iepF1[i-1] = 0; iepCN1[i-1] = 0; iepCK1[i-1] = 0;  iepB1[i-1] = 0; 

    //gXYhitsF0[i] = fs->make<TGraph>(0);
    //gXYhitsF0[i]->SetNameTitle(Form("gXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    //gXYhitsCN0[i] = fs->make<TGraph>(0);
    //gXYhitsCN0[i]->SetNameTitle(Form("gXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    //gXYhitsCK0[i] = fs->make<TGraph>(0);
    //gXYhitsCK0[i]->SetNameTitle(Form("gXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    //ixydF0[i-1] = 0; ixydCN0[i-1] = 0; ixydCK0[i-1] = 0;     

    //gXYhitsF1[i] = fs->make<TGraph>(0);
    //gXYhitsF1[i]->SetNameTitle(Form("gXYhitsF1_layer_%02d",i),Form("HitsF1 in XY for layer %d",i));
    //gXYhitsCN1[i] = fs->make<TGraph>(0);
    //gXYhitsCN1[i]->SetNameTitle(Form("gXYhitsCN1_layer_%02d",i),Form("HitsCN1 in XY for layer %d",i));
    //gXYhitsCK1[i] = fs->make<TGraph>(0);
    //gXYhitsCK1[i]->SetNameTitle(Form("gXYhitsCK1_layer_%02d",i),Form("HitsCK1 in XY for layer %d",i));
    //ixydF1[i-1] = 0; ixydCN1[i-1] = 0; ixydCK1[i-1] = 0;     

  }

  for(int i=1;i<=50;i++)
    hNHxXYhitsF[i] = fs->make<TH2D>(Form("hNHxXYhitsF_layer_%02d",i),Form("NHx HitsF in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hNHxXYhitsCN[i] = fs->make<TH2D>(Form("hNHxXYhitsCN_layer_%02d",i),Form("NHx HitsCN in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hNHxXYhitsCK[i] = fs->make<TH2D>(Form("hNHxXYhitsCK_layer_%02d",i),Form("NHx HitsCK in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);

  hXYmissedhits = fs->make<TH2D>("hXYmissedhits","hXYmissedhits", 600, -300., 300., 600, -300., 300.);
  hXYLowELosshitsF = fs->make<TH2D>("hXYLowELosshitsF","hXYLowELosshitsF", 600, -300., 300., 600, -300., 300.);
  hXYLowELosshitsCN = fs->make<TH2D>("hXYLowELosshitsCN","hXYLowELosshitsCN", 600, -300., 300., 600, -300., 300.);
  hXYLowELosshitsCK = fs->make<TH2D>("hXYLowELosshitsCK","hXYLowELosshitsCK", 600, -300., 300., 600, -300., 300.);

  hYZmissedhits = fs->make<TH2D>("hYZmissedhits","hYZmissedhits", 250, 300., 550., 300, 0., 300.);
  hYZLowELosshitsF = fs->make<TH2D>("hYZLowELosshitsF","hYZLowELosshitsF", 250, 300., 550., 300, 0., 300.);
  hYZLowELosshitsCN = fs->make<TH2D>("hYZLowELosshitsCN","hYZLowELosshitsCN", 250, 300., 550., 300, 0., 300.);
  hYZLowELosshitsCK = fs->make<TH2D>("hYZLowELosshitsCK","hYZLowELosshitsCK", 250, 300., 550., 300, 0., 300.);
  hYZLLowELosshitsHEFCN = fs->make<TH2D>("hYZLLowELosshitsHEFCN","hYZLLowELosshitsHEFCN", 600, -50., 550., 350, -50., 300.);

  
  hXLowELosshitsHEFCN = fs->make<TH1D>("hXLowELosshitsHEFCN","hXLowELosshitsHEFCN", 600, -300., 300.);
  hYLowELosshitsHEFCN = fs->make<TH1D>("hYLowELosshitsHEFCN","hYLowELosshitsHEFCN", 600, -300., 300.);
  hZLowELosshitsHEFCN = fs->make<TH1D>("hZLowELosshitsHEFCN","hZLowELosshitsHEFCN", 2400, -1200., 1200.);
  
  //hYZhits = fs->make<TH2D>("hYZhits","hYZhits", 1200, -600., 600., 1200, -600., 600.);
  //hXYhits = fs->make<TH2D>("hXYhits","Hits in XY", 600, -300., 300., 600, -300., 300.);

  hYZhitsEE = fs->make<TH2D>("hYZhitsEE","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsHEF = fs->make<TH2D>("hYZhitsHEF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsHEB = fs->make<TH2D>("hYZhitsHEB","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  
  hYZhitsEEF = fs->make<TH2D>("hYZhitsEEF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsEECN = fs->make<TH2D>("hYZhitsEECN","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsEECK = fs->make<TH2D>("hYZhitsEECK","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);

  hYZhitsHEFF = fs->make<TH2D>("hYZhitsHEFF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsHEFCN = fs->make<TH2D>("hYZhitsHEFCN","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hYZhitsHEFCK = fs->make<TH2D>("hYZhitsHEFCK","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);

  hRHTXYhits = fs->make<TH2D>("hRHTXYhits","Hits in XY", 600, -300., 300., 600, -300., 300.);
  hRHTYZhitsEE = fs->make<TH2D>("hRHTYZhitsEE","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsHEF = fs->make<TH2D>("hRHTYZhitsHEF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsHEB = fs->make<TH2D>("hRHTYZhitsHEB","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsEEF = fs->make<TH2D>("hRHTYZhitsEEF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsEECN = fs->make<TH2D>("hRHTYZhitsEECN","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsEECK = fs->make<TH2D>("hRHTYZhitsEECK","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsHEFF = fs->make<TH2D>("hRHTYZhitsHEFF","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsHEFCN = fs->make<TH2D>("hRHTYZhitsHEFCN","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTYZhitsHEFCK = fs->make<TH2D>("hRHTYZhitsHEFCK","Hits in YZ plane for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  
  hRHTRZhitsEE = fs->make<TH2D>("hRHTRZhitsEE","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsHEF = fs->make<TH2D>("hRHTRZhitsHEF","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsHEB = fs->make<TH2D>("hRHTRZhitsHEB","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsEEF = fs->make<TH2D>("hRHTRZhitsEEF","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsEECN = fs->make<TH2D>("hRHTRZhitsEECN","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsEECK = fs->make<TH2D>("hRHTRZhitsEECK","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsHEFF = fs->make<TH2D>("hRHTRZhitsHEFF","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsHEFCN = fs->make<TH2D>("hRHTRZhitsHEFCN","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  hRHTRZhitsHEFCK = fs->make<TH2D>("hRHTRZhitsHEFCK","Hits for R_{xy} vs z-axis for |X| < 20 cm", 250, 300., 550., 300, 0., 300.);
  
  hRHTGlbRZhitsF = fs->make<TH2D>("hRHTGlbRZhitsF","Hits for R_{xy} vs z-axis", 250, 300., 550., 300, 0., 300.);
  hRHTGlbRZhitsCN = fs->make<TH2D>("hRHTGlbRZhitsCN","Hits for R_{xy} vs z-axis", 250, 300., 550., 300, 0., 300.);
  hRHTGlbRZhitsCK = fs->make<TH2D>("hRHTGlbRZhitsCK","Hits for R_{xy} vs z-axis", 250, 300., 550., 300, 0., 300.);
  hRHTGlbRZhitsSci = fs->make<TH2D>("hRHTGlbRZhitsSci","Hits for R_{xy} vs z-axis", 250, 300., 550., 300, 0., 300.);

  hDiffX = fs->make<TH1D>("hDiffX" , "Difference of x-position (testHGCalGeometry - RecHitTools)" , 200 , -20 , 20 );
  hDiffX->GetXaxis()->SetTitle("x-axis (cm)");
  hDiffY = fs->make<TH1D>("hDiffY" , "Difference of y-position (testHGCalGeometry - RecHitTools)" , 200 , -20 , 20 );
  hDiffY->GetXaxis()->SetTitle("y-axis (cm)");
  hDiffZ = fs->make<TH1D>("hDiffZ" , "Difference of z-position (testHGCalGeometry - RecHitTools)" , 200 , -20 , 20 );
  hDiffZ->GetXaxis()->SetTitle("z-axis (cm)");

  hCellThickness = fs->make<TH1D>("hCellThickness" , "Cell Thickness" , 500 , 0 , 500 );
  hDiffZ->GetXaxis()->SetTitle("thickness (#mum)");

  evt = 0;
  winfo.clear();

  name = iConfig.getParameter<std::string>("Detector");
  //name = iConfig.getUntrackedParameter<string>("Detector", "");
  geomToken_ = esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", name});
  
  
  caloGeomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
}


CellHitSum::~CellHitSum()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CellHitSum::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  if(evt == 0){
    std::ifstream fin("/afs/cern.ch/work/i/idas/CMSSW/CMSSW_12_1_X_2021-10-24-2300/src/ReadSimResult/wafertype/wafer.csv");
    std::string s;
    waferinfo wafer;
    //std::cout << "evt : " << evt << std::endl;
    while(std::getline(fin,s)){
      //std::cout << "line " << s.data() << std::endl;
      sscanf(s.c_str(),"%d,%d,%d,%d",&wafer.layer,&wafer.u,&wafer.v,&wafer.type);
      //printf("%d | %d | %d | %d\n",wafer.layer,wafer.u,wafer.v,wafer.type);
      //wafer.layer = 
      winfo.push_back(wafer);
    }
    fin.close();
  }
  evt++;
  
  Handle<SimTrackContainer> simtrack;
  iEvent.getByToken(tSimTrackContainer, simtrack);
  int itrk = 0;
  //double muonpt = 0.0;
  SimTrackContainer::const_iterator itTrack;
  for(itTrack = simtrack->begin(); itTrack != simtrack->end(); ++itTrack) {
    // int charge = itTrack->charge();
    int charge = itTrack->charge();  
    hCharge->Fill( charge );
    if(!itTrack->noGenpart()){
      hPt->Fill(itTrack->momentum().pt());
      hEta->Fill(itTrack->momentum().eta());
      hPhi->Fill(itTrack->momentum().phi());
    }
    hPDG->Fill(itTrack->type());
    
    if(itTrack->noGenpart())
      hPtNoGen->Fill(itTrack->momentum().pt());
    
    // if(itTrack->charge()==0){
    //   printf("index : %d, (id,type,charge) : (%06u,%5d,%04.1f), (vert,genpat) : (%02d, %02d), \n\t\t(pt,eta,phi,e) : (%7.5lf, %7.5lf, %7.5lf,%7.5lf), Surface (x,y,z) : (%7.5lf,%7.5lf,%7.5lf)\n",
    // 	     itrk,itTrack->trackId(), itTrack->type(), itTrack->charge(),
    // 	     itTrack->vertIndex(), itTrack->genpartIndex(),
    // 	     itTrack->momentum().pt(),itTrack->momentum().eta(),itTrack->momentum().phi(),itTrack->momentum().e(),
    // 	     itTrack->trackerSurfacePosition().x(),itTrack->trackerSurfacePosition().y(),itTrack->trackerSurfacePosition().z()
    // 	     );
    // }
    itrk++;
  }
  

  const CaloGeometry &geomCalo = iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(geomCalo);
  
  const auto& geomR = iSetup.getData(geomToken_);
  const HGCalGeometry* geom = &geomR;
  //DetId::Detector det;
  if (geom->topology().waferHexagon6()) {
    ForwardSubdetector subdet;
    if (name == "HGCalHESiliconSensitive")
      subdet = HGCHEF;
    else if (name == "HGCalHEScintillatorSensitive")
      subdet = HGCHEB;
    else
      subdet = HGCEE;
    std::cout << "a) Perform test for " << name << " Detector:Subdetector " << DetId::Forward << ":" << subdet << " Mode "
              << geom->topology().dddConstants().geomMode() << std::endl;
  }else{
    //DetId::Detector det;
    // if (name == "HGCalHESiliconSensitive")
    //   det = DetId::HGCalHSi;
    // else if (name == "HGCalHEScintillatorSensitive")
    //   det = DetId::HGCalHSc;
    // else
    //   det = DetId::HGCalEE;
    // std::cout << "b) Perform test for " << name << " Detector " << det << " Mode "
    //           << geom->topology().dddConstants().geomMode() << std::endl;
  }
  
  std::map<uint32_t, std::pair<hitsinfo, energysum> > map_hits;
  map_hits.clear();
  unsigned int nofSiHits = 0;
  Handle<PCaloHitContainer> simhit;
  iEvent.getByToken(tSimCaloHitContainer, simhit);
  for(PCaloHitContainer::const_iterator itHit= simhit->begin(); itHit!= simhit->end(); ++itHit) {

    if(name == "HGCalEESensitive" or name == "HGCalHESiliconSensitive"){
      
      HGCSiliconDetId id(itHit->id());
      
      if(name == "HGCalEESensitive"){
	hELossEE->Fill(itHit->energy()*1.e6);
	if(id.type()==HGCSiliconDetId::HGCalFine)
	  hELossEEF->Fill(itHit->energy()*1.e6); //in keV
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	  hELossEECN->Fill(itHit->energy()*1.e6); //in keV
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	  hELossEECK->Fill(itHit->energy()*1.e6); //in keV
      }

      if(name == "HGCalHESiliconSensitive"){
	hELossHEF->Fill(itHit->energy()*1.e6);
	if(id.type()==HGCSiliconDetId::HGCalFine)
	  hELossHEFF->Fill(itHit->energy()*1.e6); //in keV
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	  hELossHEFCN->Fill(itHit->energy()*1.e6); //in keV
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	  hELossHEFCK->Fill(itHit->energy()*1.e6); //in keV
      }
    }
    
    if (name == "HGCalHEScintillatorSensitive")
      hELossHEB->Fill(itHit->energy()*1.e6);
    
    DetId id1 = static_cast<DetId>(itHit->id());
    GlobalPoint global2 = rhtools_.getPosition(id1);
    double RXY = TMath::Sqrt(global2.x()*global2.x() + global2.y()*global2.y());
    
    // std::cout << "DetId (" << det << ": position ("<< global2.x() << ", " << global2.y() << ", " << global2.z() 
    // 	      << "), Si thickness "<< rhtools_.getSiThickness(id1) 
    // 	      << ", IsSi "<< rhtools_.isSilicon(id1)
    // 	      << ", IsSci "<< rhtools_.isScintillator(id1)
    // 	      << ", Layer1 "<< rhtools_.getLayer(id1)
    // 	      << ", Layer2 "<< rhtools_.getLayerWithOffset(id1)
    // 	      << ", lastLayerEE  "<< rhtools_.lastLayerEE()
    // 	      << ", lastLayerFH  "<< rhtools_.lastLayerFH()
    // 	      << ", firstLayerBH  "<< rhtools_.firstLayerBH()
    // 	      << ", lastLayerBH  "<< rhtools_.lastLayerBH()
    // 	      << ", lastLayer  "<< rhtools_.lastLayer()
    // 	      << std::endl;
    
    //if((rhtools_.isSilicon(id1) or rhtools_.isScintillator(id1)) and TMath::Abs(global2.x())<20.0){
    if(rhtools_.isSilicon(id1) or rhtools_.isScintillator(id1)){
      
      if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),120.,1.e-7))
	hRHTGlbRZhitsF->Fill(TMath::Abs(global2.z()), RXY);
      else if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),200.,1.e-7))
	hRHTGlbRZhitsCN->Fill(TMath::Abs(global2.z()), RXY);
      else if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),300.,1.e-7))
	hRHTGlbRZhitsCK->Fill(TMath::Abs(global2.z()), RXY);	    
      else
	hRHTGlbRZhitsSci->Fill(TMath::Abs(global2.z()), RXY);
      
      //if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),0.,1.e-7))
	

    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if(rhtools_.isSilicon(id1) or rhtools_.isScintillator(id1)){
      uint32_t id_ = itHit->id();
      
      energysum esum;
      hitsinfo hinfo;
    
      if (map_hits.count(id_) != 0) {
    	hinfo = map_hits[id_].first;
    	esum = map_hits[id_].second;
      } else {
    	hinfo.hitid =  nofSiHits;
    	hinfo.x = global2.x();
    	hinfo.y = global2.y();
    	hinfo.z = global2.z();
    	// hinfo.sector = sector;
    	// hinfo.sector2 = subsector;
    	// hinfo.cell = rhtools_.getCell(id1).first ;
    	// hinfo.cell2 = rhtools_.getCell(id1).second ;
    	// hinfo.type = type;
    	hinfo.layer = rhtools_.getLayerWithOffset(id1);
    	hinfo.phi = rhtools_.getPhi(id1);
    	hinfo.eta = rhtools_.getEta(id1);
    	for(itTrack = simtrack->begin(); itTrack != simtrack->end(); ++itTrack) {
    	  if(itTrack->trackId() == UInt_t(itHit->geantTrackId())){ //!itTrack->noGenpart() to confirm that it is coming from muon
    	    hinfo.trkpt = itTrack->momentum().pt();
    	    hinfo.trketa = itTrack->momentum().eta();
    	    hinfo.trkphi = itTrack->momentum().phi();
    	    hinfo.charge = itTrack->charge();
    	    hinfo.pdg = itTrack->type();
    	  }
    	}
      }
      esum.etotal += itHit->energy();
      hinfo.nhits++;

      HepGeom::Point3D<float> gcoord = HepGeom::Point3D<float>(global2.x(), global2.y(), global2.z());
      double tof = (gcoord.mag() * CLHEP::cm) / CLHEP::c_light;
      double time = itHit->time() ;
      time -= tof ; 
      
      for (unsigned int k = 0; k < 2; ++k) {
        if (time > 0 && time < 25.)
          esum.eTime[k] += itHit->energy();
    	else{
    	  esum.eTime[k+2] += itHit->energy();
    	  // printf("Failed Det : %s, trackid : %d, ihit : %d, layer : %d, id : %u, time-tof : %lf, Eloss : %5.2lf (keV), (x,y,z) : (%lf,%lf,%lf)\n", 
    	  // 	 name.c_str(), itHit->geantTrackId(), nofSiHits, hinfo.layer, id_, time, itHit->energy()*1.e6, hinfo.x, hinfo.y, hinfo.z);
    	}
      }
      
      // if (verbosity_ > 1)
      //   edm::LogVerbatim("HGCalValidation") << " -----------------------   gx = " << hinfo.x << " gy = " << hinfo.y
      //                                       << " gz = " << hinfo.z << " phi = " << hinfo.phi << " eta = " << hinfo.eta;
      
      // printf("Det : %s, ihit:%d, layer:%d, id:%u, time : %5.3lf, tof : %5.3lf, ti-to : %5.3lf, Eloss : %5.2lf (keV), (x,y,z) : (%5.2lf,%5.2lf,%5.2lf)\n", 
      // 	     name.c_str(), nofSiHits, hinfo.layer, id_, itHit->time(), tof, time, itHit->energy()*1.e6, hinfo.x, hinfo.y, hinfo.z);
      
      map_hits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);
      nofSiHits++;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    hEtaCell->Fill(rhtools_.getEta(id1));
    hPhiCell->Fill(rhtools_.getPhi(id1));
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    GlobalPoint global1 = geom->getPosition(id1);
    
    if (geom->topology().valid(id1)) {

      
      //std::cout << "DetId (" << det << ": position ("<< global1.x() << ", " << global1.y() << ", " << global1.z() << ") " << std::endl;
      
      //hYZhits->Fill(global1.z(),global1.y());
      if(TMath::Abs(global1.x())<20.0){

	if(name == "HGCalEESensitive" or name == "HGCalHESiliconSensitive"){
	  HGCSiliconDetId id(itHit->id());
	  
	  if(name == "HGCalEESensitive"){
	    hYZhitsEE->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalFine)
	      hYZhitsEEF->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	      hYZhitsEECN->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	      hYZhitsEECK->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));	    
	  }

	  if(name == "HGCalHESiliconSensitive"){
	    hYZhitsHEF->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalFine)
	      hYZhitsHEFF->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	      hYZhitsHEFCN->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
	    if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	      hYZhitsHEFCK->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));	    	    
	  }
	  
	}
	
	if(name == "HGCalHEScintillatorSensitive")
	  hYZhitsHEB->Fill(TMath::Abs(global1.z()),TMath::Abs(global1.y()));
      }


      
      /// Using rechit tools
      //===============================================================================
      if(TMath::Abs(global2.x())<20.0){
	
	if(rhtools_.isSilicon(id1)){

	  if(rhtools_.getLayerWithOffset(id1) <= rhtools_.lastLayerEE()){
	    
	    hRHTYZhitsEE->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	    hRHTRZhitsEE->Fill(TMath::Abs(global2.z()), RXY);

	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),120.,1.e-7)){
	      hRHTYZhitsEEF->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	      hRHTRZhitsEEF->Fill(TMath::Abs(global2.z()), RXY);
	    }
	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),200.,1.e-7)){
	      hRHTYZhitsEECN->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	      hRHTRZhitsEECN->Fill(TMath::Abs(global2.z()), RXY);
	    }
	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),300.,1.e-7)){
	      hRHTYZhitsEECK->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));	    
	      hRHTRZhitsEECK->Fill(TMath::Abs(global2.z()), RXY);	    
	    }

	  }else{
	    
	    hRHTYZhitsHEF->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	    hRHTRZhitsHEF->Fill(TMath::Abs(global2.z()), RXY);

	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),120.,1.e-7)){
	      hRHTYZhitsHEFF->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	      hRHTRZhitsHEFF->Fill(TMath::Abs(global2.z()), RXY);
	    }
	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),200.,1.e-7)){
	      hRHTYZhitsHEFCN->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	      hRHTRZhitsHEFCN->Fill(TMath::Abs(global2.z()), RXY);
	    }
	    if(TMath::AreEqualAbs(rhtools_.getSiThickness(id1),300.,1.e-7)){
	      hRHTYZhitsHEFCK->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));	    	    
	      hRHTRZhitsHEFCK->Fill(TMath::Abs(global2.z()), RXY);	    	    
	    }
	  }
	  
	}//is Si
	
	if(rhtools_.isScintillator(id1)){
	  hRHTYZhitsHEB->Fill(TMath::Abs(global2.z()),TMath::Abs(global2.y()));
	  hRHTRZhitsHEB->Fill(TMath::Abs(global2.z()), RXY);
	}      

      }
      //===============================================================================
      hRHTXYhits->Fill(global2.x(),global2.y());	    	    

      hXYhits[rhtools_.getLayerWithOffset(id1)]->Fill(global2.x(),global2.y());
      
      if(rhtools_.isSilicon(id1)){
	HGCSiliconDetId id(itHit->id());
	HGCalDetId hid(itHit->id());
	int il = rhtools_.getLayerWithOffset(id1);
	if(id.type()==HGCSiliconDetId::HGCalFine){
	  
	  if(global1.z()<0.0){
	    grXYhitsF0[il]->SetPoint(ixyF0[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsF0[il]->SetPoint(iepF0[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsF0[il]->Fill(global2.x(),global2.y());
	    // hELossLayerF0[il]->Fill(itHit->energy()*1.0e6);
	  }else{
	    grXYhitsF1[il]->SetPoint(ixyF1[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsF1[il]->SetPoint(iepF1[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsF1[il]->Fill(global2.x(),global2.y());
	    // hELossLayerF1[il]->Fill(itHit->energy()*1.0e6);
	  }

	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	  
	  if(global1.z()<0.0){
	    grXYhitsCN0[il]->SetPoint(ixyCN0[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsCN0[il]->SetPoint(iepCN0[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsCN0[il]->Fill(global2.x(),global2.y());
	    // hELossLayerCN0[il]->Fill(itHit->energy()*1.0e6);
	  }else{
	    grXYhitsCN1[il]->SetPoint(ixyCN1[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsCN1[il]->SetPoint(iepCN1[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsCN1[il]->Fill(global2.x(),global2.y());
	    // hELossLayerCN1[il]->Fill(itHit->energy()*1.0e6);
	  }

	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	  
	  if(global1.z()<0.0){
	    grXYhitsCK0[il]->SetPoint(ixyCK0[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsCK0[il]->SetPoint(iepCK0[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsCK0[il]->Fill(global2.x(),global2.y());
	    // hELossLayerCK0[il]->Fill(itHit->energy()*1.0e6);
	  }else{
	    grXYhitsCK1[il]->SetPoint(ixyCK1[il-1]++,global1.x(),global1.y());
	    grEtaPhihitsCK1[il]->SetPoint(iepCK1[il-1]++,global1.eta(), global1.phi());
	    // hXYhitsCK1[il]->Fill(global2.x(),global2.y());
	    // hELossLayerCK1[il]->Fill(itHit->energy()*1.0e6);
	  }

	}
	//The following line by Pruthvi to number the cells with id0 and id1
        if(rhtools_.getCell(id1).first + rhtools_.getCell(id1).second <= 2){
          if(global1.z()<0.0)
            grXYhitsAR0[il]->SetPoint(ixyAR0[il-1]++,global1.x(),global1.y());
          else
            grXYhitsAR1[il]->SetPoint(ixyAR1[il-1]++,global1.x(),global1.y());
        }
      }else if(rhtools_.isScintillator(id1)){
	
	//HGCScintillatorDetId id(itHit->id());
        int il = rhtools_.getLayerWithOffset(id1);
	
	if (global1.z() < 0.0){
          grXYhitsB0[il]->SetPoint(ixyB0[il]++, global1.x(), global1.y());
          grEtaPhihitsB0[il]->SetPoint(iepB0[il]++, global1.eta(), global1.phi());
	  // hXYhitsB0[il]->Fill(global1.x(),global1.y());
	  // hELossLayerB0[il]->Fill(itHit->energy()*1.0e6);
	}else{
          grXYhitsB1[il]->SetPoint(ixyB1[il]++, global1.x(), global1.y());
          grEtaPhihitsB1[il]->SetPoint(iepB1[il]++, global1.eta(), global1.phi());
	  // hXYhitsB1[il]->Fill(global1.x(),global1.y());
	  // hELossLayerB1[il]->Fill(itHit->energy()*1.0e6);
	}
	
      }//is Sci or Si
    }//is valid detID
    
    
    hDiffX->Fill(global1.x()-global2.x());
    hDiffY->Fill(global1.y()-global2.y());
    hDiffZ->Fill(global1.z()-global2.z());
    
  }//Loop over HGCAL hits

  //std::cout << "simhit size : " << simhit->size() << ", nof hits in Si : " << nofSiHits << ", map size : " << map_hits.size() << std::endl;
  
  bool isPWafer = false;
  bool isFWafer = false;
  
  std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr;
  for (itr = map_hits.begin(); itr != map_hits.end(); ++itr) {
    //uint32_t id_ = (*itr).first;
    hitsinfo hinfo = (*itr).second.first;
    energysum esum = (*itr).second.second;
    hELossDQMEqV[hinfo.layer]->Fill(esum.eTime[0]);
    
    // printf("\tCellSummed : Det : %s, first hit : %d, nhits : %u, id : %u, Edep : %5.2lf (keV), (x,y,z) : (%5.2lf,%5.2lf,%5.2lf)\n", 
    //  	   name.c_str(), hinfo.hitid, hinfo.nhits, (*itr).first, esum.eTime[0]*1.e6, hinfo.x, hinfo.y, hinfo.z);
    
    HGCSiliconDetId id((*itr).first);

    // DetId id1 = static_cast<DetId>((*itr).first);
    
    // isPWafer = false;
    // isFWafer = false;
    // if(rhtools_.isSilicon(id1)){
    //   for(unsigned int iw = 0 ; iw < winfo.size() ; iw++){
    // 	if(hinfo.layer == winfo[iw].layer and rhtools_.getWafer(id1).first == winfo[iw].u and rhtools_.getWafer(id1).second == winfo[iw].v){
    // 	  if(winfo[iw].type == 0)
    // 	    isPWafer = true;
    // 	  if(winfo[iw].type == 1)
    // 	    isFWafer = true;
    // 	}
    //   }
    // }
    
    DetId id1 = static_cast<DetId>((*itr).first);
    GlobalPoint global1 = geom->getPosition(id1);
    
    if (geom->topology().valid(id1)) {

      if(rhtools_.isSilicon(id1)){
	HGCSiliconDetId id((*itr).first);
	int il = rhtools_.getLayerWithOffset(id1);
	if(id.type()==HGCSiliconDetId::HGCalFine){
	  
	  if(global1.z()<0.0){
	    hXYhitsF0[il]->Fill(global1.x(),global1.y());
	    hEPhitsF0[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerF0[il]->Fill(esum.etotal*1.0e6);
	  }else{
	    hXYhitsF1[il]->Fill(global1.x(),global1.y());
	    hEPhitsF1[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerF1[il]->Fill(esum.etotal*1.0e6);
	  }
	  
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	  
	  if(global1.z()<0.0){
	    hXYhitsCN0[il]->Fill(global1.x(),global1.y());
	    hEPhitsCN0[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerCN0[il]->Fill(esum.etotal*1.0e6);
	  }else{
	    hXYhitsCN1[il]->Fill(global1.x(),global1.y());
	    hEPhitsCN1[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerCN1[il]->Fill(esum.etotal*1.0e6);
	  }
	  
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	  
	  if(global1.z()<0.0){
	    hXYhitsCK0[il]->Fill(global1.x(),global1.y());
	    hEPhitsCK0[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerCK0[il]->Fill(esum.etotal*1.0e6);
	  }else{
	    hXYhitsCK1[il]->Fill(global1.x(),global1.y());
	    hEPhitsCK1[il]->Fill(global1.eta(), global1.phi());
	    hELossLayerCK1[il]->Fill(esum.etotal*1.0e6);
	  }

	}
      }else if(rhtools_.isScintillator(id1)){
	
	//HGCScintillatorDetId id(itHit->id());
        int il = rhtools_.getLayerWithOffset(id1);	
	if (global1.z() < 0.0){
	  hXYhitsB0[il]->Fill(global1.x(),global1.y());
	  hEPhitsB0[il]->Fill(global1.eta(), global1.phi());
	  hELossLayerB0[il]->Fill(esum.etotal*1.0e6);
	}else{
	  hXYhitsB1[il]->Fill(global1.x(),global1.y());
	  hEPhitsB1[il]->Fill(global1.eta(), global1.phi());
	  hELossLayerB1[il]->Fill(esum.etotal*1.0e6);
	}
	
      }//is Sci or Si
    }//valid detid
    //Now for invalid detids
    else{
      
      if(rhtools_.isSilicon(id1)){
	HGCSiliconDetId id((*itr).first);
	int il = rhtools_.getLayerWithOffset(id1);
	if(id.type()==HGCSiliconDetId::HGCalFine){
	  
	  if(global1.z()<0.0){
	    hXYFailhitsF0[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsF0[il]->Fill(global1.eta(), global1.phi());
	  }else{
	    hXYFailhitsF1[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsF1[il]->Fill(global1.eta(), global1.phi());
	  }
	  
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	  
	  if(global1.z()<0.0){
	    hXYFailhitsCN0[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsCN0[il]->Fill(global1.eta(), global1.phi());
	  }else{
	    hXYFailhitsCN1[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsCN1[il]->Fill(global1.eta(), global1.phi());
	  }
	  
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	  
	  if(global1.z()<0.0){
	    hXYFailhitsCK0[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsCK0[il]->Fill(global1.eta(), global1.phi());
	  }else{
	    hXYFailhitsCK1[il]->Fill(global1.x(),global1.y());
	    hEPFailhitsCK1[il]->Fill(global1.eta(), global1.phi());
	  }

	}
      }else if(rhtools_.isScintillator(id1)){
	
	//HGCScintillatorDetId id(itHit->id());
        int il = rhtools_.getLayerWithOffset(id1);	
	if (global1.z() < 0.0){
	  hXYFailhitsB0[il]->Fill(global1.x(),global1.y());
	  hEPFailhitsB0[il]->Fill(global1.eta(), global1.phi());
	}else{
	  hXYFailhitsB1[il]->Fill(global1.x(),global1.y());
	  hEPFailhitsB1[il]->Fill(global1.eta(), global1.phi());
	}
	
      }//is Sci or Si

    }//Invalid detid
    
    if(!TMath::AreEqualAbs(esum.eTime[0]*1.e6,0.0,1.e-5)){
      if(name == "HGCalEESensitive"){
  	hELossCSinBunchEE->Fill(esum.eTime[0]*1.e6);
  	if(id.type()==HGCSiliconDetId::HGCalFine){
  	  hELossCSinBunchEEF->Fill(esum.eTime[0]*1.e6); //in keV
	}
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
  	  hELossCSinBunchEECN->Fill(esum.eTime[0]*1.e6); //in keV
	}
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
  	  hELossCSinBunchEECK->Fill(esum.eTime[0]*1.e6); //in keV
	}
      }
    
      if(name == "HGCalHESiliconSensitive"){
  	hELossCSinBunchHEF->Fill(esum.eTime[0]*1.e6);
  	if(id.type()==HGCSiliconDetId::HGCalFine){
  	  hELossCSinBunchHEFF->Fill(esum.eTime[0]*1.e6); //in keV
  	  if(esum.eTime[0]*1.e6 < 35.){
  	    hXYLowELosshitsF->Fill(hinfo.x,hinfo.y);
  	    hYZLowELosshitsF->Fill(TMath::Abs(hinfo.z),TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y));
  	  }
  	}
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
  	  hELossCSinBunchHEFCN->Fill(esum.eTime[0]*1.e6); //in keV
  	  if(TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y) > 45.0 and TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y) < 60.0 and hinfo.layer >= 38)
  	    hELossCSinBunchHEFCNNoise->Fill(esum.eTime[0]*1.e6); //in keV
  	  else
  	    hELossCSinBunchHEFCNFiltered->Fill(esum.eTime[0]*1.e6); //in keV
  	  if(esum.eTime[0]*1.e6 < 35.){
  	    hPtLowELoss->Fill(hinfo.trkpt);
  	    hEtaLowELoss->Fill(hinfo.trketa);
  	    hPhiLowELoss->Fill(hinfo.trkphi);
  	    hChargeLowELoss->Fill(hinfo.charge);
  	    hPDGLowELoss->Fill(hinfo.pdg);
  	    hXYLowELosshitsCN->Fill(hinfo.x,hinfo.y);
  	    hYZLowELosshitsCN->Fill(TMath::Abs(hinfo.z),TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y));
  	    hYZLLowELosshitsHEFCN->Fill(TMath::Abs(hinfo.z),TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y));
  	    hXLowELosshitsHEFCN->Fill(hinfo.x);
  	    hYLowELosshitsHEFCN->Fill(hinfo.y);
  	    if(TMath::Abs(hinfo.x)<20.0 && TMath::Abs(hinfo.y)<20.0)
  	      hZLowELosshitsHEFCN->Fill(hinfo.z);
  	  }
  	}
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
  	  hELossCSinBunchHEFCK->Fill(esum.eTime[0]*1.e6); //in keV
  	  if(esum.eTime[0]*1.e6 < 10.){
  	    hXYLowELosshitsCK->Fill(hinfo.x,hinfo.y);
  	    hYZLowELosshitsCK->Fill(TMath::Abs(hinfo.z),TMath::Sqrt(hinfo.x*hinfo.x + hinfo.y*hinfo.y));
  	  }
  	}
      }
    }
    
    if(!TMath::AreEqualAbs(esum.eTime[2]*1.e6,0.0,1.e-5)){
      if(name == "HGCalEESensitive"){
  	hELossCSmissedEE->Fill(esum.eTime[2]*1.e6);
  	if(id.type()==HGCSiliconDetId::HGCalFine)
  	  hELossCSmissedEEF->Fill(esum.eTime[2]*1.e6); //in keV
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
  	  hELossCSmissedEECN->Fill(esum.eTime[2]*1.e6); //in keV
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
  	  hELossCSmissedEECK->Fill(esum.eTime[2]*1.e6); //in keV
  	hXYmissedhits->Fill(hinfo.x,hinfo.y);
  	hYZmissedhits->Fill(TMath::Abs(hinfo.z),TMath::Abs(hinfo.y));
      }
      
      if(name == "HGCalHESiliconSensitive"){
  	hELossCSmissedHEF->Fill(esum.eTime[2]*1.e6);
  	if(id.type()==HGCSiliconDetId::HGCalFine)
  	  hELossCSmissedHEFF->Fill(esum.eTime[2]*1.e6); //in keV
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
  	  hELossCSmissedHEFCN->Fill(esum.eTime[2]*1.e6); //in keV
  	if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
  	  hELossCSmissedHEFCK->Fill(esum.eTime[2]*1.e6); //in keV
  	hXYmissedhits->Fill(hinfo.x,hinfo.y);
  	hYZmissedhits->Fill(TMath::Abs(hinfo.z),TMath::Abs(hinfo.y));
      }
    }
  }
  
  std::vector<uint32_t> cellMaxEdep;
  cellMaxEdep.clear();
  for(int il=1 ; il<=50 ; il++){
    double energy =  0.;
    double energy1 =  0.;
    int ncells = 0;
    int ncells1 = 0;
    uint32_t maxid = 0;
    double maxEsum = 0.0 ;
    for (itr = map_hits.begin(); itr != map_hits.end(); ++itr) {
      //uint32_t id_ = (*itr).first;
      hitsinfo hinfo = (*itr).second.first;
      energysum esum = (*itr).second.second;
      // printf("\tDet : %s, first hit : %d, nhits : %u, id : %u, Edep : %5.2lf (keV), (x,y,z) : (%lf,%lf,%lf)\n", 
      // 	   name.c_str(), hinfo.hitid, hinfo.nhits, (*itr).first, esum.etotal*1.e6, hinfo.x, hinfo.y, hinfo.z);
      
      if(hinfo.layer==il and hinfo.z < 0.){// and TMath::Abs(hinfo.pdg)==13 and hinfo.trkpt>95.0
  	//energy += esum.eTime[0];
	energy += esum.etotal;
	ncells++;
	
  	if(esum.eTime[0] > maxEsum){
  	  maxEsum = esum.eTime[0] ; 
  	  maxid = (*itr).first ; 
  	}

      }//match layer and z-direction
      
      if(hinfo.layer==il and hinfo.z > 0.){//and TMath::Abs(hinfo.pdg)==13 and hinfo.trkpt>95.0
	//energy1 += esum.eTime[0];
	energy1 += esum.etotal;
	ncells1++;
      }

    }//map loop
    
    if(maxEsum*1.e6 > 0.)
      cellMaxEdep.push_back(maxid);
    if(energy*1.e6 > 0.){
      hELossLayer0[il]->Fill(energy*1.e6); //in keV
      hNCellsLayer0[il]->Fill(ncells); 
    }
    if(energy1*1.e6 > 0.){
      hELossLayer1[il]->Fill(energy1*1.e6); //in keV
      hNCellsLayer1[il]->Fill(ncells1); 
    }

  }
  
  // waferinfo wafer;
  // winfo.push_back(wafer);
  //std::cout<<"size :: "<<winfo.size()<<std::endl;
  
  for ( unsigned int ic = 0 ; ic < cellMaxEdep.size() ; ic++ ){
    uint32_t id_ = cellMaxEdep[ic];
    energysum esum = map_hits[id_].second;
    hitsinfo hinfo = map_hits[id_].first;
    DetId id1 = static_cast<DetId>(id_);
    
    if(!rhtools_.isSilicon(id1)) continue;
    
    HGCSiliconDetId id(id_);
    HGCalDetId hid(id);
    
    isPWafer = false;
    isFWafer = false;
    for(unsigned int iw = 0 ; iw < winfo.size() ; iw++){
      if(hinfo.layer == winfo[iw].layer and rhtools_.getWafer(id1).first == winfo[iw].u and rhtools_.getWafer(id1).second == winfo[iw].v){
    	if(winfo[iw].type == 0)
    	  isPWafer = true;
    	if(winfo[iw].type == 1)
    	  isFWafer = true;
      }
    }

    // printf("\tDet : %s, wafertype : %d, layer : %d, (u,v) : (%d,%d), ishalf : %d, first hit : %d, nhits : %u, Edep : %5.2lf (keV), (x,y,z) : (%lf,%lf,%lf)\n", 
    // 	   name.c_str(), hid.waferType(), hinfo.layer, rhtools_.getWafer(id1).first, rhtools_.getWafer(id1).second, rhtools_.isHalfCell(id1), hinfo.hitid, hinfo.nhits, esum.etotal*1.e6, hinfo.x, hinfo.y, hinfo.z);

    // printf("\tDet : %s, wafertype : %d, layer : %d, (u,v) : (%d,%d), isPWafer : %d, isFWafer : %d, (x,y,z) : (%lf,%lf,%lf)\n", 
    // 	   name.c_str(), hid.waferType(), hinfo.layer, rhtools_.getWafer(id1).first, rhtools_.getWafer(id1).second, isPWafer, isFWafer, hinfo.x, hinfo.y, hinfo.z);
    
    //for
    
    if(name == "HGCalEESensitive"){
      hCellThickness->Fill(rhtools_.getSiThickness(id1));
      hELossCSMaxEE->Fill(esum.eTime[0]*1.e6);
      if(id.type()==HGCSiliconDetId::HGCalFine){	
  	hELossCSMaxEEF->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(isPWafer){
	  hNHxELossCSMaxF->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsF[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxF->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
      if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
  	hELossCSMaxEECN->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(isPWafer){
	  hNHxELossCSMaxCN->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsCN[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxCN->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
      if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
  	hELossCSMaxEECK->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(isPWafer){
	  hNHxELossCSMaxCK->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsCK[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxCK->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
    }
    
    if(name == "HGCalHESiliconSensitive"){
      hCellThickness->Fill(rhtools_.getSiThickness(id1));
      hELossCSMaxHEF->Fill(esum.eTime[0]*1.e6);
      if(id.type()==HGCSiliconDetId::HGCalFine){
  	hELossCSMaxHEFF->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(isPWafer){
	  hNHxELossCSMaxF->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsF[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxF->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxF[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
      if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
  	hELossCSMaxHEFCN->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(esum.eTime[0]*1.e6 < 30. and esum.eTime[0]*1.e6 > 10.)
	  hXYhitsLELCN[hinfo.layer]->Fill(hinfo.x,hinfo.y); 
	else
	  hXYhitsHELCN[hinfo.layer]->Fill(hinfo.x,hinfo.y); 
	if(isPWafer){
	  hNHxELossCSMaxCN->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsCN[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxCN->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxCN[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
      if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
  	hELossCSMaxHEFCK->Fill(esum.eTime[0]*1.e6); //in keV
	hELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6); //in keV
	if(esum.eTime[0]*1.e6 < 10. )
	  hXYhitsLELCK[hinfo.layer]->Fill(hinfo.x,hinfo.y); 
	else if(esum.eTime[0]*1.e6 > 50. )
	  hXYhitsHELCK[hinfo.layer]->Fill(hinfo.x,hinfo.y); 
	if(isPWafer ){
	  hNHxELossCSMaxCK->Fill(esum.eTime[0]*1.e6); 
	  hNHxELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	  hNHxXYhitsCK[hinfo.layer]->Fill(hinfo.x,hinfo.y);
	}
	if(isFWafer){
	  hHxELossCSMaxCK->Fill(esum.eTime[0]*1.e6); 
	  hHxELCSMaxCK[hinfo.layer]->Fill(esum.eTime[0]*1.e6);
	}
      }
    }
    
  }

  // cellMaxEdep.clear();
  map_hits.clear();

  // #ifdef THIS_IS_AN_EVENT_EXAMPLE
  //    Handle<ExampleData> pIn;
  //    iEvent.getByLabel("example",pIn);
  // #endif

  // #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //    ESHandle<SetupData> pSetup;
  //    iSetup.get<SetupRecord>().get(pSetup);
  // #endif

  // for (const auto& track : iEvent.get(tracksToken_)) {
  //   // do something with track parameters, e.g, plot the charge.
  //   int charge = track.charge();
  //   histo->Fill( charge );
  //   hPt->Fill(track.pt());
  // }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
CellHitSum::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
CellHitSum::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CellHitSum::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CellHitSum);
