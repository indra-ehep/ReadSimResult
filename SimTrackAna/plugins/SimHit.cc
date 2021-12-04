// -*- C++ -*-
//
// Package:    Demo/SimHit
// Class:      SimHit
//
/**\class SimHit TrackAnalyzer.cc Track/TrackAnalyzer/plugins/TrackAnalyzer.cc

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

#include "Geometry/HGCalCommonData/interface/HGCalGeometryMode.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

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

class SimHit : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
  
  explicit SimHit(const edm::ParameterSet&);
  ~SimHit();

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
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> tok_hgcal_;
  
  TH1D **hELCSFWF ;
  TH1D **hELCSFWCN ;
  TH1D **hELCSFWCK ;

  TH1D **hELCSPWF ;
  TH1D **hELCSPWCN ;
  TH1D **hELCSPWCK ;

  //FW = Full Wafer
  TH2D **hXYhitsFWF;
  TH2D **hXYhitsFWCN;
  TH2D **hXYhitsFWCK;

  TH2D **hXYhitsB;

  TH2D **hXYhitsPWF;
  TH2D **hXYhitsPWCN;
  TH2D **hXYhitsPWCK;

  // For rechittool z positions. The 0 and 1 are for -ve and +ve, respectively.
  TGraph **grXYhitsFWF0;
  TGraph **grXYhitsFWCN0;
  TGraph **grXYhitsFWCK0;
  TGraph **grXYhitsB0;
  int ixyFWF0[50], ixyFWCN0[50], ixyFWCK0[50], ixyB0[50];

  TGraph **grXYhitsFWF1;
  TGraph **grXYhitsFWCN1;
  TGraph **grXYhitsFWCK1;
  TGraph **grXYhitsB1;
  int ixyFWF1[50], ixyFWCN1[50], ixyFWCK1[50], ixyB1[50];

  TGraph **grXYhitsPWF0;
  TGraph **grXYhitsPWCN0;
  TGraph **grXYhitsPWCK0;
  int ixyPWF0[50], ixyPWCN0[50], ixyPWCK0[50];

  TGraph **grXYhitsPWF1;
  TGraph **grXYhitsPWCN1;
  TGraph **grXYhitsPWCK1;
  int ixyPWF1[50], ixyPWCN1[50], ixyPWCK1[50];
  /////////////////////////////////
  
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
SimHit::SimHit(const edm::ParameterSet& iConfig)
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
  
  hELCSFWF =  new TH1D*[50]; 
  hELCSFWCN =  new TH1D*[50]; 
  hELCSFWCK =  new TH1D*[50]; 
  for(int i=1;i<=50;i++)
    hELCSFWF[i] = fs->make<TH1D>(Form("hELCSFWF_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSFWCN[i] = fs->make<TH1D>(Form("hELCSFWCN_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSFWCK[i] = fs->make<TH1D>(Form("hELCSFWCK_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  
  hELCSPWF =  new TH1D*[50]; 
  hELCSPWCN =  new TH1D*[50]; 
  hELCSPWCK =  new TH1D*[50]; 
  for(int i=1;i<=50;i++)
    hELCSPWF[i] = fs->make<TH1D>(Form("hELCSPWF_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSPWCN[i] = fs->make<TH1D>(Form("hELCSPWCN_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  for(int i=1;i<=50;i++)
    hELCSPWCK[i] = fs->make<TH1D>(Form("hELCSPWCK_layer_%02d",i),Form("Energy for layer %d",i), 500, 0., 500.);
  
  // for 50 layers in earch +/- z-direction
  hXYhitsFWF =  new TH2D*[50]; 
  hXYhitsFWCN =  new TH2D*[50]; 
  hXYhitsFWCK =  new TH2D*[50]; 
  hXYhitsB =  new TH2D*[50]; 
  for(int i=1;i<=50;i++)
    hXYhitsFWF[i] = fs->make<TH2D>(Form("hXYhitsFWF_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsFWCN[i] = fs->make<TH2D>(Form("hXYhitsFWCN_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsFWCK[i] = fs->make<TH2D>(Form("hXYhitsFWCK_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsB[i] = fs->make<TH2D>(Form("hXYhitsB_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);

  hXYhitsPWF =  new TH2D*[50]; 
  hXYhitsPWCN =  new TH2D*[50]; 
  hXYhitsPWCK =  new TH2D*[50]; 
  for(int i=1;i<=50;i++)
    hXYhitsPWF[i] = fs->make<TH2D>(Form("hXYhitsPWF_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsPWCN[i] = fs->make<TH2D>(Form("hXYhitsPWCN_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  for(int i=1;i<=50;i++)
    hXYhitsPWCK[i] = fs->make<TH2D>(Form("hXYhitsPWCK_layer_%02d",i),Form("Hits in XY for layer %d",i), 600, -300., 300., 600, -300., 300.);
  
  grXYhitsFWF0 =  new TGraph*[50]; 
  grXYhitsFWCN0 =  new TGraph*[50]; 
  grXYhitsFWCK0 =  new TGraph*[50]; 
  grXYhitsB0 =  new TGraph*[50]; 
  for(int i=1;i<=50;i++){
    grXYhitsFWF0[i] = fs->make<TGraph>(0);
    grXYhitsFWF0[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grXYhitsFWCN0[i] = fs->make<TGraph>(0);
    grXYhitsFWCN0[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grXYhitsFWCK0[i] = fs->make<TGraph>(0);
    grXYhitsFWCK0[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    grXYhitsB0[i] = fs->make<TGraph>(0);
    grXYhitsB0[i]->SetNameTitle(Form("grXYhitsB0_layer_%02d",i),Form("HitsB0 in XY for layer %d",i));
    ixyFWF0[i-1] = 0; ixyFWCN0[i-1] = 0; ixyFWCK0[i-1] = 0; ixyB0[i-1] = 0; 
  }

  grXYhitsFWF1 =  new TGraph*[50]; 
  grXYhitsFWCN1 =  new TGraph*[50]; 
  grXYhitsFWCK1 =  new TGraph*[50]; 
  grXYhitsB1 =  new TGraph*[50]; 
  for(int i=1;i<=50;i++){
    grXYhitsFWF1[i] = fs->make<TGraph>(0);
    grXYhitsFWF1[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grXYhitsFWCN1[i] = fs->make<TGraph>(0);
    grXYhitsFWCN1[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grXYhitsFWCK1[i] = fs->make<TGraph>(0);
    grXYhitsFWCK1[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    grXYhitsB1[i] = fs->make<TGraph>(0);
    grXYhitsB1[i]->SetNameTitle(Form("grXYhitsB0_layer_%02d",i),Form("HitsB0 in XY for layer %d",i));
    ixyFWF1[i-1] = 0; ixyFWCN1[i-1] = 0; ixyFWCK1[i-1] = 0; ixyB1[i-1] = 0; 
  }

  grXYhitsPWF0 =  new TGraph*[50]; 
  grXYhitsPWCN0 =  new TGraph*[50]; 
  grXYhitsPWCK0 =  new TGraph*[50]; 
  for(int i=1;i<=50;i++){
    grXYhitsPWF0[i] = fs->make<TGraph>(0);
    grXYhitsPWF0[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grXYhitsPWCN0[i] = fs->make<TGraph>(0);
    grXYhitsPWCN0[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grXYhitsPWCK0[i] = fs->make<TGraph>(0);
    grXYhitsPWCK0[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    ixyPWF0[i-1] = 0; ixyPWCN0[i-1] = 0; ixyPWCK0[i-1] = 0; 
  }

  grXYhitsPWF1 =  new TGraph*[50]; 
  grXYhitsPWCN1 =  new TGraph*[50]; 
  grXYhitsPWCK1 =  new TGraph*[50]; 
  for(int i=1;i<=50;i++){
    grXYhitsPWF1[i] = fs->make<TGraph>(0);
    grXYhitsPWF1[i]->SetNameTitle(Form("grXYhitsF0_layer_%02d",i),Form("HitsF0 in XY for layer %d",i));
    grXYhitsPWCN1[i] = fs->make<TGraph>(0);
    grXYhitsPWCN1[i]->SetNameTitle(Form("grXYhitsCN0_layer_%02d",i),Form("HitsCN0 in XY for layer %d",i));
    grXYhitsPWCK1[i] = fs->make<TGraph>(0);
    grXYhitsPWCK1[i]->SetNameTitle(Form("grXYhitsCK0_layer_%02d",i),Form("HitsCK0 in XY for layer %d",i));
    ixyPWF1[i-1] = 0; ixyPWCN1[i-1] = 0; ixyPWCK1[i-1] = 0; 
  }
  
  
  evt = 0;
  winfo.clear();

  name = iConfig.getParameter<std::string>("Detector");
  //name = iConfig.getUntrackedParameter<string>("Detector", "");
  geomToken_ = esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", name});
    
  caloGeomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

  //tok_hgcal_ = esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(edm::ESInputTag{"", name});
  tok_hgcal_ = esConsumes<HGCalDDDConstants, IdealGeometryRecord>(edm::ESInputTag{"", name});
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
}


SimHit::~SimHit()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimHit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  
  // Handle<SimTrackContainer> simtrack;
  // iEvent.getByToken(tSimTrackContainer, simtrack);

  const HGCalDDDConstants& hgcons = iSetup.getData(tok_hgcal_);
  
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
    
    DetId id1 = static_cast<DetId>(itHit->id());
    GlobalPoint global = rhtools_.getPosition(id1);
    
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
    	hinfo.x = global.x();
    	hinfo.y = global.y();
    	hinfo.z = global.z();
    	hinfo.layer = rhtools_.getLayerWithOffset(id1);
    	hinfo.phi = rhtools_.getPhi(id1);
    	hinfo.eta = rhtools_.getEta(id1);
      }
      esum.etotal += itHit->energy();
      hinfo.nhits++;

      HepGeom::Point3D<float> gcoord = HepGeom::Point3D<float>(global.x(), global.y(), global.z());
      double tof = (gcoord.mag() * CLHEP::cm) / CLHEP::c_light;
      double time = itHit->time() ;
      time -= tof ; 
      
      for (unsigned int k = 0; k < 2; ++k) {
        if (time > 0 && time < 25.)
          esum.eTime[k] += itHit->energy();
    	else
    	  esum.eTime[k+2] += itHit->energy();
      }
      
      // printf("Det : %s, ihit:%d, layer:%d, id:%u, time : %5.3lf, tof : %5.3lf, ti-to : %5.3lf, Eloss : %5.2lf (keV), (x,y,z) : (%5.2lf,%5.2lf,%5.2lf)\n", 
      // 	     name.c_str(), nofSiHits, hinfo.layer, id_, itHit->time(), tof, time, itHit->energy()*1.e6, hinfo.x, hinfo.y, hinfo.z);
      
      map_hits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);
      nofSiHits++;
    }
    
    //if (geom->topology().valid(id1))
          
  }// hit loop
  //std::cout << "simhit size : " << simhit->size() << ", nof hits in Si : " << nofSiHits << ", map size : " << map_hits.size() << std::endl;
  
  bool isPWafer = false;
  bool isFWafer = false;
  int partialType = -1;
  
  std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr;
  for (itr = map_hits.begin(); itr != map_hits.end(); ++itr) {

    hitsinfo hinfo = (*itr).second.first;
    int ilayer = hinfo.layer;
    energysum esum = (*itr).second.second;
    double edep = esum.eTime[0] * CLHEP::GeV / CLHEP::keV; //index 0 and 1 corresponds to 25 ns and 1000 ns, respectively. In addititon, chaging energy loss unit to keV.

    if(TMath::AreEqualAbs(edep,0.0,1.e-5)) continue;

    DetId id1 = static_cast<DetId>((*itr).first);
    
    isPWafer = false;
    isFWafer = false;
    partialType = -1;
    if(rhtools_.isSilicon(id1)){
      for(unsigned int iw = 0 ; iw < winfo.size() ; iw++){
    	if(hinfo.layer == winfo[iw].layer and rhtools_.getWafer(id1).first == winfo[iw].u and rhtools_.getWafer(id1).second == winfo[iw].v){
    	  if(winfo[iw].type == 0)
    	    isPWafer = true;
    	  if(winfo[iw].type == 1)
    	    isFWafer = true;
    	}
      }
      
      int type, part, orient;
      HGCSiliconDetId detId = HGCSiliconDetId((*itr).first);
      std::tie(type, part, orient) = hgcons.waferType(detId) ;
      partialType = part ;
    }
    
    if(name == "HGCalEESensitive" or name == "HGCalHESiliconSensitive"){
      HGCSiliconDetId id((*itr).first);
      if(partialType==0){
	if(id.type()==HGCSiliconDetId::HGCalFine){
	  hELCSFWF[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsFWF0[ilayer]->SetPoint(ixyFWF0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsFWF1[ilayer]->SetPoint(ixyFWF1[ilayer-1]++,hinfo.x,hinfo.y);
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	  hELCSFWCN[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsFWCN0[ilayer]->SetPoint(ixyFWCN0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsFWCN1[ilayer]->SetPoint(ixyFWCN1[ilayer-1]++,hinfo.x,hinfo.y);
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
	  hELCSFWCK[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsFWCK0[ilayer]->SetPoint(ixyFWCK0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsFWCK1[ilayer]->SetPoint(ixyFWCK1[ilayer-1]++,hinfo.x,hinfo.y);
	}
      }//partialType == 0
    
      if(isFWafer){
	if(hinfo.z>0.0){
	  if(id.type()==HGCSiliconDetId::HGCalFine) hXYhitsFWF[ilayer]->Fill(hinfo.x,hinfo.y);
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThin) hXYhitsFWCN[ilayer]->Fill(hinfo.x,hinfo.y);
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThick) hXYhitsFWCK[ilayer]->Fill(hinfo.x,hinfo.y);
	}
      }//isFWafer
      
      if(partialType > 0){
	if(id.type()==HGCSiliconDetId::HGCalFine){
	  hELCSPWF[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsPWF0[ilayer]->SetPoint(ixyPWF0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsPWF1[ilayer]->SetPoint(ixyPWF1[ilayer-1]++,hinfo.x,hinfo.y);
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	  hELCSPWCN[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsPWCN0[ilayer]->SetPoint(ixyPWCN0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsPWCN1[ilayer]->SetPoint(ixyPWCN1[ilayer-1]++,hinfo.x,hinfo.y);
	}
	if(id.type()==HGCSiliconDetId::HGCalCoarseThick){
	  hELCSPWCK[ilayer]->Fill(edep);
	  if(hinfo.z<0.0)
	    grXYhitsPWCK0[ilayer]->SetPoint(ixyPWCK0[ilayer-1]++,hinfo.x,hinfo.y);
	  else
	    grXYhitsPWCK1[ilayer]->SetPoint(ixyPWCK1[ilayer-1]++,hinfo.x,hinfo.y);
	}
      }//partialType > 0

      if(isPWafer){
	if(hinfo.z>0.0){
	  if(id.type()==HGCSiliconDetId::HGCalFine) hXYhitsPWF[ilayer]->Fill(hinfo.x,hinfo.y);
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThin) hXYhitsPWCN[ilayer]->Fill(hinfo.x,hinfo.y);
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThick) hXYhitsPWCK[ilayer]->Fill(hinfo.x,hinfo.y);
	}      
      }
    }else{//scintillator
      
      if(hinfo.z<0.0){
	grXYhitsB0[ilayer]->SetPoint(ixyB0[ilayer-1]++,hinfo.x,hinfo.y);
	hXYhitsB[ilayer]->Fill(hinfo.x,hinfo.y);
      }else
	grXYhitsB1[ilayer]->SetPoint(ixyB1[ilayer-1]++,hinfo.x,hinfo.y);
      
    }
    
  }
  
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
SimHit::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SimHit::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimHit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(SimHit);
