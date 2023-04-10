/**********************************************************************
 Created on : 06/10/2022
 Purpose    : Read the eta,phi plots from geantoutput
 Author     : Indranil Das, Visiting Fellow, TIFR
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

#include <string>
#include <iostream>

#include <TFile.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

int HistEtaPhi(int refLayer=41)
{
  //int AddGraphs(TGraph*& gr1, TGraph*& gr2);

  bool isEE = (refLayer<=26) ? true : false;
  bool isSCE = (refLayer>33) ? true : false;
  //int index[4] = {1,2,3,9}; //D92
  int maxfile = 20;
  //string geom = "D94" ;
  string geom = "D88" ;
  //string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_4_0_pre4/Extended2026"+geom+"/";
  string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_5_0_pre5/Extended2026"+geom+"/";
  //string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_6_X_2022-09-27-2300/Extended2026"+geom+"/";
  //string inputPath = "/home/indra/temp/Geant"+geom+"_L"+Form("%d",refLayer)+"/";

  //maxfile = 29;
  
  //TGraph *hEPtotF0,*hEPtotCN0, *hEPtotCK0, *hEPtotB0;
  TH2D *hEPGentotF0 = 0x0,*hEPGentotCN0 = 0x0, *hEPGentotCK0 = 0x0, *hEPGentotB0 = 0x0;
  TH2D *hEPGenFailtotF0 = 0x0,*hEPGenFailtotCN0 = 0x0, *hEPGenFailtotCK0 = 0x0, *hEPGenFailtotB0 = 0x0;
  TH2D *hEPRectotF0 = 0x0,*hEPRectotCN0 = 0x0, *hEPRectotCK0 = 0x0, *hEPRectotB0 = 0x0;
  TH2D *hEPRecFailtotF0 = 0x0,*hEPRecFailtotCN0 = 0x0, *hEPRecFailtotCK0 = 0x0, *hEPRecFailtotB0 = 0x0;

  for(int i=0;i<maxfile;i++){
    //int idx = index[i];
    int idx = i;
    string fname = inputPath+"geantoutput"+geom+"_"+to_string(idx)+".root";
    TFile *fGen = TFile::Open(fname.c_str());
    cout << "Openning file " << fGen <<", name : " << fGen->GetName() << endl;
    if(isEE){
      cout << "EE : Openning file " << fGen <<", name : " << fGen->GetName() << endl;
      TH2D *hEtaPhihitsF0 = (TH2D *)fGen->Get(Form("prodEE/hEPhitsF0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsCN0 = (TH2D *)fGen->Get(Form("prodEE/hEPhitsCN0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsCK0 = (TH2D *)fGen->Get(Form("prodEE/hEPhitsCK0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsB0 = (TH2D *)fGen->Get(Form("prodEE/hEPhitsB0_layer_%02d",refLayer));

      TH2D *hEtaPhiFailhitsF0 = (TH2D *)fGen->Get(Form("prodEE/hEPFailhitsF0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fGen->Get(Form("prodEE/hEPFailhitsCN0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fGen->Get(Form("prodEE/hEPFailhitsCK0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsB0 = (TH2D *)fGen->Get(Form("prodEE/hEPFailhitsB0_layer_%02d",refLayer));

      // TH2D *hEtaPhihitsF0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPhitsF0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsCN0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPhitsCN0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsCK0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPhitsCK0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsB0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPhitsB0_layer_%02d",refLayer));

      // TH2D *hEtaPhiFailhitsF0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPFailhitsF0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPFailhitsCN0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPFailhitsCK0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsB0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumEE/hEPFailhitsB0_layer_%02d",refLayer));

      if(i==0){
  	hEPGentotF0 = (TH2D *)hEtaPhihitsF0->Clone("hEPGentotF0") ;
  	hEPGentotCN0 = (TH2D *)hEtaPhihitsCN0->Clone("hEPGentotCN0") ;
  	hEPGentotCK0 = (TH2D *)hEtaPhihitsCK0->Clone("hEPGentotCK0") ;
  	hEPGentotB0 = (TH2D *)hEtaPhihitsB0->Clone("hEPGentotB0") ;

  	hEPGenFailtotF0 = (TH2D *)hEtaPhiFailhitsF0->Clone("hEPGenFailtotF0") ;
  	hEPGenFailtotCN0 = (TH2D *)hEtaPhiFailhitsCN0->Clone("hEPGenFailtotCN0") ;
  	hEPGenFailtotCK0 = (TH2D *)hEtaPhiFailhitsCK0->Clone("hEPGenFailtotCK0") ;
  	hEPGenFailtotB0 = (TH2D *)hEtaPhiFailhitsB0->Clone("hEPGenFailtotB0") ;
      }else {
  	hEPGentotF0->Add(hEtaPhihitsF0);
  	hEPGentotCN0->Add(hEtaPhihitsCN0);
  	hEPGentotCK0->Add(hEtaPhihitsCK0);
  	hEPGentotB0->Add(hEtaPhihitsB0);

  	hEPGenFailtotF0->Add(hEtaPhiFailhitsF0);
  	hEPGenFailtotCN0->Add(hEtaPhiFailhitsCN0);
  	hEPGenFailtotCK0->Add(hEtaPhiFailhitsCK0);
  	hEPGenFailtotB0->Add(hEtaPhiFailhitsB0);
      }
    }else{

      TH2D *hEtaPhihitsF0 = (TH2D *)fGen->Get(Form("prodHEF/hEPhitsF0_layer_%d",refLayer));
      TH2D *hEtaPhihitsCN0 = (TH2D *)fGen->Get(Form("prodHEF/hEPhitsCN0_layer_%d",refLayer));
      TH2D *hEtaPhihitsCK0 = (TH2D *)fGen->Get(Form("prodHEF/hEPhitsCK0_layer_%d",refLayer));
      TH2D *hEtaPhihitsB0 = (TH2D *)fGen->Get(Form("prodHEB/hEPhitsB0_layer_%d",refLayer));
      cout << "HE1 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;
      
      TH2D *hEtaPhiFailhitsF0 = (TH2D *)fGen->Get(Form("prodHEF/hEPFailhitsF0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fGen->Get(Form("prodHEF/hEPFailhitsCN0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fGen->Get(Form("prodHEF/hEPFailhitsCK0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsB0 = (TH2D *)fGen->Get(Form("prodHEB/hEPFailhitsB0_layer_%d",refLayer));
      cout << "HE2 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;

      // TH2D *hEtaPhihitsF0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPhitsF0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsCN0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPhitsCN0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsCK0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPhitsCK0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsB0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEB/hEPhitsB0_layer_%d",refLayer));
      // cout << "HE1 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;
      
      // TH2D *hEtaPhiFailhitsF0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPFailhitsF0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPFailhitsCN0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEF/hEPFailhitsCK0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsB0 = (TH2D *)fGen->Get(Form("hgcalCellHitSumHEB/hEPFailhitsB0_layer_%d",refLayer));
      // cout << "HE2 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;

      if(i==0){
	cout << "HE3.1 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;
  	hEPGentotF0 = (TH2D *)hEtaPhihitsF0->Clone("hEPGentotF0") ;
  	hEPGentotCN0 = (TH2D *)hEtaPhihitsCN0->Clone("hEPGentotCN0") ;
  	hEPGentotCK0 = (TH2D *)hEtaPhihitsCK0->Clone("hEPGentotCK0") ;
  	hEPGentotB0 = (TH2D *)hEtaPhihitsB0->Clone("hEPGentotB0") ;

  	hEPGenFailtotF0 = (TH2D *)hEtaPhiFailhitsF0->Clone("hEPGenFailtotF0") ;
  	hEPGenFailtotCN0 = (TH2D *)hEtaPhiFailhitsCN0->Clone("hEPGenFailtotCN0") ;
  	hEPGenFailtotCK0 = (TH2D *)hEtaPhiFailhitsCK0->Clone("hEPGenFailtotCK0") ;
  	hEPGenFailtotB0 = (TH2D *)hEtaPhiFailhitsB0->Clone("hEPGenFailtotB0") ;
      }else {
	cout << "HE3.2 : Openning file " << fGen <<", name : " << fGen->GetName() << endl;
  	hEPGentotF0->Add(hEtaPhihitsF0);
  	hEPGentotCN0->Add(hEtaPhihitsCN0);
  	hEPGentotCK0->Add(hEtaPhihitsCK0);
  	hEPGentotB0->Add(hEtaPhihitsB0);

  	hEPGenFailtotF0->Add(hEtaPhiFailhitsF0);
  	hEPGenFailtotCN0->Add(hEtaPhiFailhitsCN0);
  	hEPGenFailtotCK0->Add(hEtaPhiFailhitsCK0);
  	hEPGenFailtotB0->Add(hEtaPhiFailhitsB0);	
      }

    }//is HE
    
    hEPGentotCN0->SetTitle(Form("GenHits in layer %d (z < 0.0 cm)",refLayer));
  }//file loop

  
  for(int i=0;i<maxfile;i++){
    //int idx = index[i];
    int idx = i;
    string fname = inputPath+"hgcRecHit"+geom+"_"+to_string(idx)+".root";
    //string fname = inputPath+"recoutput"+geom+".root_"+to_string(idx)+".root";
    TFile *fRec = TFile::Open(fname.c_str());
    cout << "Openning file " << fRec <<", name : " << fRec->GetName() << endl;
    if(isEE){
      
      TH2D *hEtaPhihitsF0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPhitsF0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsCN0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPhitsCN0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsCK0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPhitsCK0_layer_%02d",refLayer));
      TH2D *hEtaPhihitsB0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPhitsB0_layer_%02d",refLayer));

      TH2D *hEtaPhiFailhitsF0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPFailhitsF0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPFailhitsCN0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPFailhitsCK0_layer_%02d",refLayer));
      TH2D *hEtaPhiFailhitsB0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyEE/hEPFailhitsB0_layer_%02d",refLayer));

      // TH2D *hEtaPhihitsF0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPhitsF0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsCN0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPhitsCN0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsCK0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPhitsCK0_layer_%02d",refLayer));
      // TH2D *hEtaPhihitsB0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPhitsB0_layer_%02d",refLayer));

      // TH2D *hEtaPhiFailhitsF0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPFailhitsF0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPFailhitsCN0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPFailhitsCK0_layer_%02d",refLayer));
      // TH2D *hEtaPhiFailhitsB0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyEE/hEPFailhitsB0_layer_%02d",refLayer));

      if(i==0){
  	hEPRectotF0 = (TH2D *)hEtaPhihitsF0->Clone("hEPRectotF0") ;
  	hEPRectotCN0 = (TH2D *)hEtaPhihitsCN0->Clone("hEPRectotCN0") ;
  	hEPRectotCK0 = (TH2D *)hEtaPhihitsCK0->Clone("hEPRectotCK0") ;
  	hEPRectotB0 = (TH2D *)hEtaPhihitsB0->Clone("hEPRectotB0") ;

  	hEPRecFailtotF0 = (TH2D *)hEtaPhiFailhitsF0->Clone("hEPRecFailtotF0") ;
  	hEPRecFailtotCN0 = (TH2D *)hEtaPhiFailhitsCN0->Clone("hEPRecFailtotCN0") ;
  	hEPRecFailtotCK0 = (TH2D *)hEtaPhiFailhitsCK0->Clone("hEPRecFailtotCK0") ;
  	hEPRecFailtotB0 = (TH2D *)hEtaPhiFailhitsB0->Clone("hEPRecFailtotB0") ;
      }else {
  	hEPRectotF0->Add(hEtaPhihitsF0);
  	hEPRectotCN0->Add(hEtaPhihitsCN0);
  	hEPRectotCK0->Add(hEtaPhihitsCK0);
  	hEPRectotB0->Add(hEtaPhihitsB0);

  	hEPRecFailtotF0->Add(hEtaPhiFailhitsF0);
  	hEPRecFailtotCN0->Add(hEtaPhiFailhitsCN0);
  	hEPRecFailtotCK0->Add(hEtaPhiFailhitsCK0);
  	hEPRecFailtotB0->Add(hEtaPhiFailhitsB0);
      }
    }else{
      TH2D *hEtaPhihitsF0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPhitsF0_layer_%d",refLayer));
      TH2D *hEtaPhihitsCN0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPhitsCN0_layer_%d",refLayer));
      TH2D *hEtaPhihitsCK0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPhitsCK0_layer_%d",refLayer));
      TH2D *hEtaPhihitsB0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyBH/hEPhitsB0_layer_%d",refLayer));

      TH2D *hEtaPhiFailhitsF0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPFailhitsF0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPFailhitsCN0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyFH/hEPFailhitsCK0_layer_%d",refLayer));
      TH2D *hEtaPhiFailhitsB0 = (TH2D *)fRec->Get(Form("hgcalRecHitStudyBH/hEPFailhitsB0_layer_%d",refLayer));      

      // TH2D *hEtaPhihitsF0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPhitsF0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsCN0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPhitsCN0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsCK0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPhitsCK0_layer_%d",refLayer));
      // TH2D *hEtaPhihitsB0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyBH/hEPhitsB0_layer_%d",refLayer));

      // TH2D *hEtaPhiFailhitsF0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPFailhitsF0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsCN0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPFailhitsCN0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsCK0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyFH/hEPFailhitsCK0_layer_%d",refLayer));
      // TH2D *hEtaPhiFailhitsB0 = (TH2D *)fRec->Get(Form("hgcalMTRecoStudyBH/hEPFailhitsB0_layer_%d",refLayer));      

      if(i==0){
  	hEPRectotF0 = (TH2D *)hEtaPhihitsF0->Clone("hEPRectotF0") ;
  	hEPRectotCN0 = (TH2D *)hEtaPhihitsCN0->Clone("hEPRectotCN0") ;
  	hEPRectotCK0 = (TH2D *)hEtaPhihitsCK0->Clone("hEPRectotCK0") ;
  	hEPRectotB0 = (TH2D *)hEtaPhihitsB0->Clone("hEPRectotB0") ;

  	hEPRecFailtotF0 = (TH2D *)hEtaPhiFailhitsF0->Clone("hEPRecFailtotF0") ;
  	hEPRecFailtotCN0 = (TH2D *)hEtaPhiFailhitsCN0->Clone("hEPRecFailtotCN0") ;
  	hEPRecFailtotCK0 = (TH2D *)hEtaPhiFailhitsCK0->Clone("hEPRecFailtotCK0") ;
  	hEPRecFailtotB0 = (TH2D *)hEtaPhiFailhitsB0->Clone("hEPRecFailtotB0") ;
      }else {
  	hEPRectotF0->Add(hEtaPhihitsF0);
  	hEPRectotCN0->Add(hEtaPhihitsCN0);
  	hEPRectotCK0->Add(hEtaPhihitsCK0);
  	hEPRectotB0->Add(hEtaPhihitsB0);

  	hEPRecFailtotF0->Add(hEtaPhiFailhitsF0);
  	hEPRecFailtotCN0->Add(hEtaPhiFailhitsCN0);
  	hEPRecFailtotCK0->Add(hEtaPhiFailhitsCK0);
  	hEPRecFailtotB0->Add(hEtaPhiFailhitsB0);	
      }

    }//is HE
    
    hEPRectotCN0->SetTitle(Form("RecHits in layer %d (z < 0.0 cm)",refLayer));
  }//file loop

  hEPGentotF0->SetMarkerColor(kRed);
  hEPGentotCN0->SetMarkerColor(kGreen+1);
  hEPGentotCK0->SetMarkerColor(kMagenta);
  hEPGentotB0->SetMarkerColor(kBlue);
  hEPGenFailtotF0->SetMarkerColor(kBlack);
  hEPGenFailtotCN0->SetMarkerColor(kRed+2);
  hEPGenFailtotCK0->SetMarkerColor(kBlack);
  hEPGenFailtotB0->SetMarkerColor(kBlack);
  
  hEPRectotF0->SetMarkerColor(kRed);
  hEPRectotCN0->SetMarkerColor(kGreen+1);
  hEPRectotCK0->SetMarkerColor(kMagenta);
  hEPRectotB0->SetMarkerColor(kBlue);
  hEPRecFailtotF0->SetMarkerColor(kBlack);
  hEPRecFailtotCN0->SetMarkerColor(kRed+2);
  hEPRecFailtotCK0->SetMarkerColor(kBlack);
  hEPRecFailtotB0->SetMarkerColor(kBlack);
  
  hEPGentotF0->SetMarkerStyle(kFullCircle);
  hEPGentotCN0->SetMarkerStyle(kFullCircle);
  hEPGentotCK0->SetMarkerStyle(kFullCircle);
  hEPGentotB0->SetMarkerStyle(kFullCircle);
  hEPGenFailtotF0->SetMarkerStyle(kFullCircle);
  hEPGenFailtotCN0->SetMarkerStyle(kFullCircle);
  hEPGenFailtotCK0->SetMarkerStyle(kFullCircle);
  hEPGenFailtotB0->SetMarkerStyle(kFullCircle);
  
  hEPRectotF0->SetMarkerStyle(kFullCircle);
  hEPRectotCN0->SetMarkerStyle(kFullCircle);
  hEPRectotCK0->SetMarkerStyle(kFullCircle);
  hEPRectotB0->SetMarkerStyle(kFullCircle);
  hEPRecFailtotF0->SetMarkerStyle(kFullCircle);
  hEPRecFailtotCN0->SetMarkerStyle(kFullCircle);
  hEPRecFailtotCK0->SetMarkerStyle(kFullCircle);
  hEPRecFailtotB0->SetMarkerStyle(kFullCircle);
  
  int rebin = 5; //1 nominal or 5 for DQMlike plots
  hEPGentotF0->Rebin2D(rebin);
  hEPGentotCN0->Rebin2D(rebin);
  hEPGentotCK0->Rebin2D(rebin);
  hEPGentotB0->Rebin2D(rebin);
  hEPGenFailtotF0->Rebin2D(rebin);
  hEPGenFailtotCN0->Rebin2D(rebin);
  hEPGenFailtotCK0->Rebin2D(rebin);
  hEPGenFailtotB0->Rebin2D(rebin);
  
  hEPRectotF0->Rebin2D(rebin);
  hEPRectotCN0->Rebin2D(rebin);
  hEPRectotCK0->Rebin2D(rebin);
  hEPRectotB0->Rebin2D(rebin);
  hEPRecFailtotF0->Rebin2D(rebin);
  hEPRecFailtotCN0->Rebin2D(rebin);
  hEPRecFailtotCK0->Rebin2D(rebin);
  hEPRecFailtotB0->Rebin2D(rebin);
  
  float markersize = .3;//0.1 nominal and 0.3 for DQM like
  hEPGentotF0->SetMarkerSize(markersize);
  hEPGentotCN0->SetMarkerSize(markersize);
  hEPGentotCK0->SetMarkerSize(markersize);
  hEPGentotB0->SetMarkerSize(markersize);
  hEPGenFailtotF0->SetMarkerSize(2*markersize);
  hEPGenFailtotCN0->SetMarkerSize(2*markersize);
  hEPGenFailtotCK0->SetMarkerSize(2*markersize);
  hEPGenFailtotB0->SetMarkerSize(2*markersize);

  hEPRectotF0->SetMarkerSize(markersize);
  hEPRectotCN0->SetMarkerSize(markersize);
  hEPRectotCK0->SetMarkerSize(markersize);
  hEPRectotB0->SetMarkerSize(markersize);
  hEPRecFailtotF0->SetMarkerSize(2*markersize);
  hEPRecFailtotCN0->SetMarkerSize(2*markersize);
  hEPRecFailtotCK0->SetMarkerSize(2*markersize);
  hEPRecFailtotB0->SetMarkerSize(2*markersize);
  
  hEPGentotCN0->SetFillColor(kGreen);
  hEPGentotCK0->SetFillColor(kMagenta);
  hEPGenFailtotF0->SetFillColor(kRed+2);
  hEPGenFailtotCK0->SetFillColor(kBlack);
  TLegend *leg1 = new TLegend(0.56,0.57,0.99,0.92);
  leg1->AddEntry(hEPGentotCN0,"Valid detIDs: 200 #mum", "f");
  leg1->AddEntry(hEPGentotCK0,"Valid detIDs: 300 #mum", "f");
  leg1->AddEntry(hEPGenFailtotF0,"Failed detIDs: 200 #mum", "f");
  leg1->AddEntry(hEPGenFailtotCK0,"Failed detIDs: 300 #mum", "f");

  hEPRectotCN0->SetFillColor(kGreen);
  hEPRectotCK0->SetFillColor(kMagenta);
  hEPRecFailtotF0->SetFillColor(kRed+2);
  hEPRecFailtotCK0->SetFillColor(kBlack);
  TLegend *leg2 = new TLegend(0.56,0.57,0.99,0.92);
  leg2->AddEntry(hEPRectotCN0,"Valid detIDs: 200 #mum", "f");
  leg2->AddEntry(hEPRectotCK0,"Valid detIDs: 300 #mum", "f");
  leg2->AddEntry(hEPRecFailtotF0,"Failed detIDs: 200 #mum", "f");
  leg2->AddEntry(hEPRecFailtotCK0,"Failed detIDs: 300 #mum", "f");

  gStyle->SetOptStat(00111);
  hEPGentotCN0->SetName("Valid 200 #mum");
  hEPGentotCK0->SetName("Valid 300 #mum");
  hEPGenFailtotCK0->SetName("Failed 300 #mum");
  hEPGenFailtotCN0->SetName("Failed 200 #mum");
  hEPRectotCN0->SetName("Valid 200 #mum");
  hEPRectotCK0->SetName("Valid 300 #mum");
  hEPRecFailtotCK0->SetName("Failed 300 #mum");
  hEPRecFailtotCN0->SetName("Failed 200 #mum");
  
  //hEPGentotCN0->GetXaxis()->SetLimits(-3.2, -1.3);
  hEPGentotCN0->GetXaxis()->SetTitle("#eta");
  hEPGentotCN0->GetYaxis()->SetTitle("#phi");
  hEPRectotCN0->GetXaxis()->SetTitle("#eta");
  hEPRectotCN0->GetYaxis()->SetTitle("#phi");

  hEPGentotCN0->GetXaxis()->SetTitleSize(0.06);
  hEPGentotCN0->GetXaxis()->SetTitleOffset(0.8);
  hEPGentotCN0->GetYaxis()->SetTitleSize(0.06);
  hEPGentotCN0->GetYaxis()->SetTitleOffset(0.8);
  hEPRectotCN0->GetXaxis()->SetTitleSize(0.06);
  hEPRectotCN0->GetXaxis()->SetTitleOffset(0.8);
  hEPRectotCN0->GetYaxis()->SetTitleSize(0.06);
  hEPRectotCN0->GetYaxis()->SetTitleOffset(0.8);

  
  TCanvas *c1 = new TCanvas("c1","c1",1800,800);
  c1->Divide(2,1);
  c1->cd(1);
  hEPGentotCN0->Draw();
  //hEPGentotF0->Draw("sames");
  hEPGentotCK0->Draw("sames");
  // hEPGentotB0->Draw("sames");
  hEPGenFailtotCN0->Draw("sames");
  //hEPGenFailtotF0->Draw("sames");
  hEPGenFailtotCK0->Draw("sames");
  // hEPGenFailtotB0->Draw("sames");
  leg1->Draw();
  //c1->cd(1)->Update();
  
  c1->cd(2);
  hEPRectotCN0->Draw();
  //hEPRectotF0->Draw("sames");
  hEPRectotCK0->Draw("sames");
  // hEPRectotB0->Draw("sames");
  hEPRecFailtotCN0->Draw("sames");
  //hEPRecFailtotF0->Draw("sames");
  hEPRecFailtotCK0->Draw("sames");
  // hEPRecFailtotB0->Draw("sames");
  leg2->Draw();
  
  hEPGentotCN0->GetXaxis()->SetRangeUser(-3.0, -1.3);
  hEPRectotCN0->GetXaxis()->SetRangeUser(-3.0, -1.3);
  c1->Update();

  TPaveStats *stats0 = (TPaveStats *)hEPRectotCN0->FindObject("stats");
  TPaveStats *stats1 = (TPaveStats *)hEPRectotCK0->FindObject("stats");
  TPaveStats *stats2 = (TPaveStats *)hEPRecFailtotCN0->FindObject("stats");
  TPaveStats *stats3 = (TPaveStats *)hEPRecFailtotCK0->FindObject("stats");
  stats0->SetX1NDC(.56);          stats0->SetX2NDC(.76);
  stats0->SetY1NDC(.4);           stats0->SetY2NDC(.55);
  stats1->SetX1NDC(.78);           stats1->SetX2NDC(.99);
  stats1->SetY1NDC(.4);           stats1->SetY2NDC(.55);
  stats2->SetX1NDC(.56);          stats2->SetX2NDC(.76);
  stats2->SetY1NDC(.2);           stats2->SetY2NDC(.35);
  stats3->SetX1NDC(.78);           stats3->SetX2NDC(.99);
  stats3->SetY1NDC(.2);           stats3->SetY2NDC(.35);
  stats3->Draw();
  
  c1->cd(1);
  TPaveStats *stats4 = (TPaveStats *)hEPGentotCN0->FindObject("stats");
  TPaveStats *stats5 = (TPaveStats *)hEPGentotCK0->FindObject("stats");
  TPaveStats *stats6 = (TPaveStats *)hEPGenFailtotCN0->FindObject("stats");
  TPaveStats *stats7 = (TPaveStats *)hEPGenFailtotCK0->FindObject("stats");
  stats4->SetX1NDC(.56);          stats4->SetX2NDC(.76);
  stats4->SetY1NDC(.4);           stats4->SetY2NDC(.55);
  stats5->SetX1NDC(.78);           stats5->SetX2NDC(.99);
  stats5->SetY1NDC(.4);           stats5->SetY2NDC(.55);
  stats6->SetX1NDC(.56);          stats6->SetX2NDC(.76);
  stats6->SetY1NDC(.2);           stats6->SetY2NDC(.35);
  stats7->SetX1NDC(.78);           stats7->SetX2NDC(.99);
  stats7->SetY1NDC(.2);           stats7->SetY2NDC(.35);
  stats7->Draw();
  
  
  return true;
}

// int AddGraphs(TGraph*& gr1, TGraph*& gr2)
// {
//   int npoints = gr1->GetN();

//   for(int i=0;i<gr2->GetN();i++)
//     gr1->SetPoint(npoints++, gr2->GetPointX(i), gr2->GetPointY(i));

//   return true;
// }
