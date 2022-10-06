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
      //double energy = it->energy();
      
      // float globalx = global.x();
      // float globaly = global.y();
      //float globalz = global.z();
      //h_RZ_->Fill(std::abs(globalz), global.perp());
      //hXYhits->Fill(global.x(),global.y());
      
      if (geom0->topology().valid(detId)) {
	//hXYhits->Fill(global1.x(),global1.y());
      
	if(rhtools_.isSilicon(detId)){
	  HGCSiliconDetId id(it.id());
	  HGCalDetId hid(it.id());
	  int il = rhtools_.getLayerWithOffset(detId);
	  if(id.type()==HGCSiliconDetId::HGCalFine){
	    //hXYhitsF[il]->Fill(global2.x(),global2.y());

	    if(global1.z()<0.0)
	      grXYhitsF0[il]->SetPoint(ixyF0[il-1]++,global1.x(),global1.y());
	    else
	      grXYhitsF1[il]->SetPoint(ixyF1[il-1]++,global1.x(),global1.y());
	    //if(id.zside() == -1)
	    //  gXYhitsF0[il]->SetPoint(ixydF0[il-1]++,global2.x(),global2.y());
	    //else
	    //  gXYhitsF1[il]->SetPoint(ixydF1[il-1]++,global2.x(),global2.y());
	    if(global1.z()<0.0)
	      grEtaPhihitsF0[il]->SetPoint(iepF0[il-1]++,global1.eta(), global1.phi());
	    else
	      grEtaPhihitsF1[il]->SetPoint(iepF1[il-1]++,global1.eta(), global1.phi());

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThin){
	    //hXYhitsCN[il]->Fill(global2.x(),global2.y());

	    if(global1.z()<0.0)
	      grXYhitsCN0[il]->SetPoint(ixyCN0[il-1]++,global1.x(),global1.y());
	    else
	      grXYhitsCN1[il]->SetPoint(ixyCN1[il-1]++,global1.x(),global1.y());
	    //if(id.zside() == -1)
	    //  gXYhitsCN0[il]->SetPoint(ixydCN0[il-1]++,global2.x(),global2.y());
	    //else
	    //  gXYhitsCN1[il]->SetPoint(ixydCN1[il-1]++,global2.x(),global2.y());
	    if(global1.z()<0.0)
	      grEtaPhihitsCN0[il]->SetPoint(iepCN0[il-1]++,global1.eta(), global1.phi());
	    else
	      grEtaPhihitsCN1[il]->SetPoint(iepCN1[il-1]++,global1.eta(), global1.phi());

	  }
	  if(id.type()==HGCSiliconDetId::HGCalCoarseThick){ //case 2 : 
	    //hXYhitsCK[il]->Fill(global2.x(),global2.y());

	    if(global1.z()<0.0)
	      grXYhitsCK0[il]->SetPoint(ixyCK0[il-1]++,global1.x(),global1.y());
	    else
	      grXYhitsCK1[il]->SetPoint(ixyCK1[il-1]++,global1.x(),global1.y());
	    //if(id.zside() == -1)
	    //  gXYhitsCK0[il]->SetPoint(ixydCK0[il-1]++,global2.x(),global2.y());
	    //else
	    //  gXYhitsCK1[il]->SetPoint(ixydCK1[il-1]++,global2.x(),global2.y());
	    if(global1.z()<0.0)
	      grEtaPhihitsCK0[il]->SetPoint(iepCK0[il-1]++,global1.eta(), global1.phi());
	    else
	      grEtaPhihitsCK1[il]->SetPoint(iepCK1[il-1]++,global1.eta(), global1.phi());

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
	
	  // if (global1.z() < 0.0)
	  //   hXYhitsB[il]->Fill(global1.x(),global1.y());

	  if (global1.z() < 0.0)
	    grXYhitsB0[il]->SetPoint(ixyB0[il]++, global1.x(), global1.y());
	  else
	    grXYhitsB1[il]->SetPoint(ixyB1[il]++, global1.x(), global1.y());

	  // if (id.zside() == -1)
	  //   gXYhitsB0[il]->SetPoint(ixydB0[il]++, global2.x(), global2.y());
	  // else
	  //   gXYhitsB1[il]->SetPoint(ixydB1[il]++, global2.x(), global2.y());
	  if (global1.z() < 0.0)
	    grEtaPhihitsB0[il]->SetPoint(iepB0[il]++, global1.eta(), global1.phi());
	  else
	    grEtaPhihitsB1[il]->SetPoint(iepB1[il]++, global1.eta(), global1.phi());

	  
	}//Silicon or scintillator

      }//valid topology

    }//loop over iterator
  }//is Valid container
  

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
