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
  
  struct hitsinfo {
    hitsinfo() {
      x = y = z = phi = eta = 0.0;
      cell = cell2 = sector = sector2 = type = layer = 0;
      hitid = nhits = 0;
      isMu = false;
    }
    double x, y, z, phi, eta;
    int cell, cell2, sector, sector2, type, layer;
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
  
  TH1D *histo, *hPt; 
  TH1D *hELossEE;
  TH1D *hELossEEF ;
  TH1D *hELossEECN ;
  TH1D *hELossEECK ;
  TH1D *hELossHEF ;
  TH1D *hELossHEFF ;
  TH1D *hELossHEFCN ;
  TH1D *hELossHEFCK ;
  TH1D *hELossHEB ;

  TH1D *hELossCellSummedEE;
  TH1D *hELossCellSummedEEF ;
  TH1D *hELossCellSummedEECN ;
  TH1D *hELossCellSummedEECK ;
  TH1D *hELossCellSummedHEF ;
  TH1D *hELossCellSummedHEFF ;
  TH1D *hELossCellSummedHEFCN ;
  TH1D *hELossCellSummedHEFCK ;
  
  TH1D **hELossDQMEqV ;
  TH1D **hELossLayer ;
  
  // TH2D *hYZhits;
  TH2D **hXYhits;

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
  histo = fs->make<TH1D>("charge" , "Charges" , 200 , -20 , 20 );
  hPt = fs->make<TH1D>("hPt" , "hPt" , 1000 , 0. , 1000. );

  hELossEE = fs->make<TH1D>("hELossEE","hELossEE", 1000, 0., 1000.);
  hELossEEF = fs->make<TH1D>("hELossEEF","hELossEEF", 1000, 0., 1000.);
  hELossEECN = fs->make<TH1D>("hELossEECN","hELossEECN", 1000, 0., 1000.);
  hELossEECK = fs->make<TH1D>("hELossEECK","hELossEECK", 1000, 0., 1000.);

  hELossHEF = fs->make<TH1D>("hELossHEF","hELossHEF", 1000, 0., 1000.);
  hELossHEFF = fs->make<TH1D>("hELossHEFF","hELossHEFF", 1000, 0., 1000.);
  hELossHEFCN = fs->make<TH1D>("hELossHEFCN","hELossHEFCN", 1000, 0., 1000.);
  hELossHEFCK = fs->make<TH1D>("hELossHEFCK","hELossHEFCK", 1000, 0., 1000.);
  
  hELossHEB = fs->make<TH1D>("hELossHEB","hELossHEB", 1000, 0., 1000.);

  hELossCellSummedEE = fs->make<TH1D>("hELossCellSummedEE","hELossCellSummedEE", 1000, 0., 1000.);
  hELossCellSummedEEF = fs->make<TH1D>("hELossCellSummedEEF","hELossCellSummedEEF", 1000, 0., 1000.);
  hELossCellSummedEECN = fs->make<TH1D>("hELossCellSummedEECN","hELossCellSummedEECN", 1000, 0., 1000.);
  hELossCellSummedEECK = fs->make<TH1D>("hELossCellSummedEECK","hELossCellSummedEECK", 1000, 0., 1000.);

  hELossCellSummedHEF = fs->make<TH1D>("hELossCellSummedHEF","hELossCellSummedHEF", 1000, 0., 1000.);
  hELossCellSummedHEFF = fs->make<TH1D>("hELossCellSummedHEFF","hELossCellSummedHEFF", 1000, 0., 1000.);
  hELossCellSummedHEFCN = fs->make<TH1D>("hELossCellSummedHEFCN","hELossCellSummedHEFCN", 1000, 0., 1000.);
  hELossCellSummedHEFCK = fs->make<TH1D>("hELossCellSummedHEFCK","hELossCellSummedHEFCK", 1000, 0., 1000.);
  
  hELossDQMEqV = new TH1D*[50]; // for 50 layers in earch +/- z-direction
  hELossLayer = new TH1D*[50]; // for 50 layers in earch +/- z-direction
  hXYhits =  new TH2D*[50]; // for 50 layers in earch +/- z-direction
  for(int i=1;i<=50;i++)
    hELossDQMEqV[i] = fs->make<TH1D>(Form("hELossDQMEqV_layer_%02d",i),Form("hELossDQMEqV_layer_%02d",i), 100, 0, 0.1);
  for(int i=1;i<=50;i++)
    hELossLayer[i] = fs->make<TH1D>(Form("hELossLayer_%02d",i),Form("hELossLayer_%02d",i), 1000, 0., 1000.);
  for(int i=1;i<=50;i++)
    hXYhits[i] = fs->make<TH2D>(Form("hXYhits_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  
  
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

    
  Handle<SimTrackContainer> simtrack;
  iEvent.getByToken(tSimTrackContainer, simtrack);
  int itrk = 0;
  for(SimTrackContainer::const_iterator itTrack = simtrack->begin(); itTrack != simtrack->end(); ++itTrack) {
    // int charge = itTrack->charge();
    int charge = itTrack->charge();  
    histo->Fill( charge );
    if(itTrack->noGenpart())
      hPt->Fill(itTrack->momentum().pt());

    // if(abs(itTrack->type())==13){
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
  DetId::Detector det;
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
    if (name == "HGCalHESiliconSensitive")
      det = DetId::HGCalHSi;
    else if (name == "HGCalHEScintillatorSensitive")
      det = DetId::HGCalHSc;
    else
      det = DetId::HGCalEE;
    std::cout << "b) Perform test for " << name << " Detector " << det << " Mode "
              << geom->topology().dddConstants().geomMode() << std::endl;
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
    if(rhtools_.isSilicon(id1)){
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
	hinfo.cell = rhtools_.getCell(id1).first ;
	hinfo.cell2 = rhtools_.getCell(id1).second ;
	// hinfo.type = type;
	hinfo.layer = rhtools_.getLayerWithOffset(id1);
	hinfo.phi = rhtools_.getPhi(id1);
	hinfo.eta = rhtools_.getEta(id1);
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
	else
	  esum.eTime[k+2] += itHit->energy();
      }
      
      // if (verbosity_ > 1)
      //   edm::LogVerbatim("HGCalValidation") << " -----------------------   gx = " << hinfo.x << " gy = " << hinfo.y
      //                                       << " gz = " << hinfo.z << " phi = " << hinfo.phi << " eta = " << hinfo.eta;
      // printf("Det : %s, trackid : %d, ihit : %d, layer : %d, id : %u, time : %lf, tof : %lf, Eloss : %5.2lf (keV), (x,y,z) : (%lf,%lf,%lf)\n", 
      // 	     name.c_str(), itHit->geantTrackId(), nofSiHits, hinfo.layer, id_, itHit->time(), tof, itHit->energy()*1.e6, hinfo.x, hinfo.y, hinfo.z);
      
      map_hits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);
      nofSiHits++;
    }
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

      hXYhits[rhtools_.getLayerWithOffset(id1)]->Fill(global1.x(),global1.y());
			 
    }
    
    
    hDiffX->Fill(global1.x()-global2.x());
    hDiffY->Fill(global1.y()-global2.y());
    hDiffZ->Fill(global1.z()-global2.z());
    
  }
  //std::cout << "simhit size : " << simhit->size() << ", nof hits in Si : " << nofSiHits << ", map size : " << map_hits.size() << std::endl;
  
  
  std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr;
  for (itr = map_hits.begin(); itr != map_hits.end(); ++itr) {
    //uint32_t id_ = (*itr).first;
    hitsinfo hinfo = (*itr).second.first;
    energysum esum = (*itr).second.second;
    hELossDQMEqV[hinfo.layer]->Fill(esum.eTime[0]);
  }
  
  std::vector<uint32_t> cellMaxEdep;
  cellMaxEdep.clear();
  for(int il=1 ; il<=50 ; il++){
    double energy =  0.;
    uint32_t maxid = 0;
    double maxEsum = 0.0 ;
    for (itr = map_hits.begin(); itr != map_hits.end(); ++itr) {
      //uint32_t id_ = (*itr).first;
      hitsinfo hinfo = (*itr).second.first;
      energysum esum = (*itr).second.second;
      // printf("\tDet : %s, first hit : %d, nhits : %u, id : %u, Edep : %5.2lf (keV), (x,y,z) : (%lf,%lf,%lf)\n", 
      // 	   name.c_str(), hinfo.hitid, hinfo.nhits, (*itr).first, esum.etotal*1.e6, hinfo.x, hinfo.y, hinfo.z);
    
      if(hinfo.layer==il and hinfo.z > 0.){
	energy += esum.eTime[0];

	if(esum.eTime[0] > maxEsum){
	  maxEsum = esum.eTime[0] ; 
	  maxid = (*itr).first ; 
	}

      }//match layer and z-direction

    }//map loop
    
    if(maxEsum > 0.)
      cellMaxEdep.push_back(maxid);
    if(energy > 0.)
      hELossLayer[il]->Fill(energy*1.e6); //in keV
  }

  for ( unsigned int ic = 0 ; ic < cellMaxEdep.size() ; ic++ ){
    uint32_t id_ = cellMaxEdep[ic];
    energysum esum = map_hits[id_].second;
    HGCSiliconDetId id(id_);
    
    if(name == "HGCalEESensitive"){
      hELossCellSummedEE->Fill(esum.eTime[0]*1.e6);
      if(id.type()==HGCSiliconDetId::HGCalFine)
	hELossCellSummedEEF->Fill(esum.eTime[0]*1.e6); //in keV
      if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	hELossCellSummedEECN->Fill(esum.eTime[0]*1.e6); //in keV
      if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	hELossCellSummedEECK->Fill(esum.eTime[0]*1.e6); //in keV
    }
    
    if(name == "HGCalHESiliconSensitive"){
      hELossCellSummedHEF->Fill(esum.eTime[0]*1.e6);
      if(id.type()==HGCSiliconDetId::HGCalFine)
	hELossCellSummedHEFF->Fill(esum.eTime[0]*1.e6); //in keV
      if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	hELossCellSummedHEFCN->Fill(esum.eTime[0]*1.e6); //in keV
      if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	hELossCellSummedHEFCK->Fill(esum.eTime[0]*1.e6); //in keV
    }
    
  }

  cellMaxEdep.clear();
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
