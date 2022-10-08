/**********************************************************************
 Created on : 06/10/2022
 Purpose    : Read the eta,phi plots from geantoutput
 Author     : Indranil Das, Visiting Fellow, TIFR
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

#include <string>
#include <iostream>

#include <TFile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>

using namespace std;

int ReadEtaPhi(int refLayer=41)
{
  int AddGraphs(TGraph*& gr1, TGraph*& gr2);

  bool isEE = (refLayer<=26) ? true : false;
  bool isSCE = (refLayer>33) ? true : false;
  //int index[4] = {1,2,3,9}; //D92
  //string inputReco = "/home/indra/temp/Results/D92_Reco/";
  //string inputReco = "/home/indra/temp/Results/D92_Gen/";
  
  // int index[4] = {0,3,4,9}; //D88
  // string inputReco = "/home/indra/temp/Results/D88_Reco/";
  //string inputReco = "/home/indra/temp/Results/D88_Gen/";
  string inputReco = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/geomval/etaphi_debug_reeval/CMSSW_12_5_0_pre5/Extended2026D88/";

  TGraph *hEPtotF0,*hEPtotCN0, *hEPtotCK0, *hEPtotB0;


  // for(int i=0;i<4;i++){
  //   //int idx = index[i];
  //   int idx = i;
  //   //string fname = inputReco+"geantoutputD92_"+to_string(idx)+".root";
  //   string fname = inputReco+"geantoutputD88_"+to_string(idx)+".root";
  //   cout << "fname : " << fname << endl;
  //   TFile *fReco = TFile::Open(fname.c_str());
  //   cout << "Openning file " << fReco <<", name : " << fReco->GetName() << endl;
  //   if(!fReco) continue;

  //   if(isEE){
  //     TGraph *hEtaPhihitsF0 = (TGraph *)fReco->Get(Form("prodEE/grEtaPhihitsF0_layer_%02d",refLayer));
  //     TGraph *hEtaPhihitsCN0 = (TGraph *)fReco->Get(Form("prodEE/grEtaPhihitsCN0_layer_%02d",refLayer));
  //     TGraph *hEtaPhihitsCK0 = (TGraph *)fReco->Get(Form("prodEE/grEtaPhihitsCK0_layer_%02d",refLayer));
  //     TGraph *hEtaPhihitsB0 = (TGraph *)fReco->Get(Form("prodEE/grEtaPhihitsB0_layer_%02d",refLayer));
  //     if(i==0){
  // 	hEPtotF0 = (TGraph *)hEtaPhihitsF0->Clone("hEPtotF0") ;
  // 	hEPtotCN0 = (TGraph *)hEtaPhihitsCN0->Clone("hEPtotCN0") ;
  // 	hEPtotCK0 = (TGraph *)hEtaPhihitsCK0->Clone("hEPtotCK0") ;
  // 	hEPtotB0 = (TGraph *)hEtaPhihitsB0->Clone("hEPtotB0") ;
  //     }else {
  // 	AddGraphs(hEPtotF0, hEtaPhihitsF0);
  // 	AddGraphs(hEPtotCN0, hEtaPhihitsCN0);
  // 	AddGraphs(hEPtotCK0, hEtaPhihitsCK0);
  // 	AddGraphs(hEPtotB0, hEtaPhihitsB0);
  //     }
  //   }else{
  //     TGraph *hEtaPhihitsF0 = (TGraph *)fReco->Get(Form("prodHEF/grEtaPhihitsF0_layer_%d",refLayer));
  //     TGraph *hEtaPhihitsCN0 = (TGraph *)fReco->Get(Form("prodHEF/grEtaPhihitsCN0_layer_%d",refLayer));
  //     TGraph *hEtaPhihitsCK0 = (TGraph *)fReco->Get(Form("prodHEF/grEtaPhihitsCK0_layer_%d",refLayer));
  //     TGraph *hEtaPhihitsB0 = (TGraph *)fReco->Get(Form("prodHEB/grEtaPhihitsB0_layer_%d",refLayer));

  //     if(i==0){
  // 	hEPtotF0 = (TGraph *)hEtaPhihitsF0->Clone("hEPtotF0") ;
  // 	hEPtotCN0 = (TGraph *)hEtaPhihitsCN0->Clone("hEPtotCN0") ;
  // 	hEPtotCK0 = (TGraph *)hEtaPhihitsCK0->Clone("hEPtotCK0") ;
  // 	hEPtotB0 = (TGraph *)hEtaPhihitsB0->Clone("hEPtotB0") ;
  //     }else {
  // 	AddGraphs(hEPtotF0, hEtaPhihitsF0);
  // 	AddGraphs(hEPtotCN0, hEtaPhihitsCN0);
  // 	AddGraphs(hEPtotCK0, hEtaPhihitsCK0);
  // 	AddGraphs(hEPtotB0, hEtaPhihitsB0);
  //     }

  //   }//is HE
    
  //   hEPtotCN0->SetTitle(Form("GenHits in layer %d (z < 0.0 cm)",refLayer));
  // }//file loop

  for(int i=0;i<3;i++){
    //int idx = index[i];
    int idx = i;
    //string fname = inputReco+"hgcRecHitD92_"+to_string(idx)+".root";
    string fname = inputReco+"hgcRecHitD88_"+to_string(idx)+".root";
    cout << "fname : " << fname << endl;
    TFile *fReco = TFile::Open(fname.c_str());
    cout << "Openning file " << fReco <<", name : " << fReco->GetName() << endl;
    if(!fReco) continue;
    if(isEE){
      TGraph *hEtaPhihitsF0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyEE/grEtaPhihitsF0_layer_%02d",refLayer));
      TGraph *hEtaPhihitsCN0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyEE/grEtaPhihitsCN0_layer_%02d",refLayer));
      TGraph *hEtaPhihitsCK0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyEE/grEtaPhihitsCK0_layer_%02d",refLayer));
      TGraph *hEtaPhihitsB0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyEE/grEtaPhihitsB0_layer_%02d",refLayer));
      if(i==0){
  	hEPtotF0 = (TGraph *)hEtaPhihitsF0->Clone("hEPtotF0") ;
  	hEPtotCN0 = (TGraph *)hEtaPhihitsCN0->Clone("hEPtotCN0") ;
  	hEPtotCK0 = (TGraph *)hEtaPhihitsCK0->Clone("hEPtotCK0") ;
  	hEPtotB0 = (TGraph *)hEtaPhihitsB0->Clone("hEPtotB0") ;
      }else {
  	AddGraphs(hEPtotF0, hEtaPhihitsF0);
  	AddGraphs(hEPtotCN0, hEtaPhihitsCN0);
  	AddGraphs(hEPtotCK0, hEtaPhihitsCK0);
  	AddGraphs(hEPtotB0, hEtaPhihitsB0);
      }
    }else{
      TGraph *hEtaPhihitsF0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyFH/grEtaPhihitsF0_layer_%d",refLayer));
      TGraph *hEtaPhihitsCN0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyFH/grEtaPhihitsCN0_layer_%d",refLayer));
      TGraph *hEtaPhihitsCK0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyFH/grEtaPhihitsCK0_layer_%d",refLayer));
      TGraph *hEtaPhihitsB0 = (TGraph *)fReco->Get(Form("hgcalRecHitStudyBH/grEtaPhihitsB0_layer_%d",refLayer));

      if(i==0){
  	hEPtotF0 = (TGraph *)hEtaPhihitsF0->Clone("hEPtotF0") ;
  	hEPtotCN0 = (TGraph *)hEtaPhihitsCN0->Clone("hEPtotCN0") ;
  	hEPtotCK0 = (TGraph *)hEtaPhihitsCK0->Clone("hEPtotCK0") ;
  	hEPtotB0 = (TGraph *)hEtaPhihitsB0->Clone("hEPtotB0") ;
      }else {
  	AddGraphs(hEPtotF0, hEtaPhihitsF0);
  	AddGraphs(hEPtotCN0, hEtaPhihitsCN0);
  	AddGraphs(hEPtotCK0, hEtaPhihitsCK0);
  	AddGraphs(hEPtotB0, hEtaPhihitsB0);
      }

    }//is HE
    
    hEPtotCN0->SetTitle(Form("RecHits in layer %d (z < 0.0 cm)",refLayer));
  }//file loop
  


  hEPtotF0->SetMarkerColor(kRed);
  hEPtotCN0->SetMarkerColor(kGreen+1);
  hEPtotCK0->SetMarkerColor(kMagenta);
  hEPtotB0->SetMarkerColor(kBlue);

  hEPtotF0->SetMarkerStyle(kFullCircle);
  hEPtotCN0->SetMarkerStyle(kFullCircle);
  hEPtotCK0->SetMarkerStyle(kFullCircle);
  hEPtotB0->SetMarkerStyle(kFullCircle);

  hEPtotF0->SetMarkerSize(0.1);
  hEPtotCN0->SetMarkerSize(0.1);
  hEPtotCK0->SetMarkerSize(0.1);
  hEPtotB0->SetMarkerSize(0.1);


  hEPtotCN0->GetXaxis()->SetLimits(-3.2, -1.3);
  hEPtotCN0->GetXaxis()->SetTitle("#eta");
  hEPtotCN0->GetYaxis()->SetTitle("#phi");
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hEPtotCN0->Draw("AP");
  hEPtotF0->Draw("P");
  hEPtotCK0->Draw("P");
  hEPtotB0->Draw("P");
  c1->SaveAs("c1.png");
  c1->SaveAs("c1.pdf");

  return true;
}

int AddGraphs(TGraph*& gr1, TGraph*& gr2)
{
  int npoints = gr1->GetN();

  for(int i=0;i<gr2->GetN();i++)
    gr1->SetPoint(npoints++, gr2->GetPointX(i), gr2->GetPointY(i));

  return true;
}
