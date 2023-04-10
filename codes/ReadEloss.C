/**********************************************************************
 Created on : 09/04/2023
 Purpose    : Read the Eloss
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

int ReadEloss(int refType=1)
{
  int maxfile = 20;
  string geom = "D94" ;
  //string geom = "D88" ;
  //string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_4_0_pre4/Extended2026"+geom+"/";
  //string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_5_0_pre5/Extended2026"+geom+"/";
  //string inputPath = "/home/indra/temp/etaphi_debug_reeval/CMSSW_12_6_X_2022-09-27-2300/Extended2026"+geom+"/";
  //string inputPath = "/home/indra/temp/Geant"+geom+"_L"+Form("%d",refLayer)+"/";
  string inputPath = "/home/indra/temp/Geant"+geom+"_All/";

  maxfile = 49;
  int maxlayer = 47;

  
  string siType = "";
  if(refType==1)
    siType = "F";
  else if(refType==2)
    siType = "CN";
  else
    siType = "CK";

  
  TH1D *hElossT = 0x0,*hElossF = 0x0, *hElossP = 0x0;
  
  for(int i=0;i<maxfile;i++){
    //int idx = index[i];
    int idx = i;
    string fname = inputPath+"geantoutput"+geom+"_"+to_string(idx)+".root";
    TFile *fGen = TFile::Open(fname.c_str());
    cout << "Openning file " << fGen <<", name : " << fGen->GetName() << endl;
    for(int il=1;il<maxlayer;il++){
      
      TH1D *hElossT120EE = (TH1D *)fGen->Get(Form("hgcalCellHitSumEE/hELCSMax%s_layer_%02d",siType.c_str(),il));
      TH1D *hElossF120EE = (TH1D *)fGen->Get(Form("hgcalCellHitSumEE/hHxELCSMax%s_layer_%02d",siType.c_str(),il));
      TH1D *hElossP120EE = (TH1D *)fGen->Get(Form("hgcalCellHitSumEE/hNHxELCSMax%s_layer_%02d",siType.c_str(),il));
      
      TH1D *hElossT120HEF = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEF/hELCSMax%s_layer_%02d",siType.c_str(),il));
      TH1D *hElossF120HEF = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEF/hHxELCSMax%s_layer_%02d",siType.c_str(),il));
      TH1D *hElossP120HEF = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEF/hNHxELCSMax%s_layer_%02d",siType.c_str(),il));

      // TH1D *hElossT120HEB = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEB/hELCSMax%s_layer_%02d",siType.c_str(),il));
      // TH1D *hElossF120HEB = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEB/hHxELCSMax%s_layer_%02d",siType.c_str(),il));
      // TH1D *hElossP120HEB = (TH1D *)fGen->Get(Form("hgcalCellHitSumHEB/hNHxELCSMax%s_layer_%02d",siType.c_str(),il));

      if(i==0 and il==1){
	
	hElossT = (TH1D *)hElossT120EE->Clone("hElossT") ;
	hElossF = (TH1D *)hElossF120EE->Clone("hElossF") ;
	hElossP = (TH1D *)hElossP120EE->Clone("hElossP") ;
	
	hElossT->Add(hElossT120HEF);
	hElossF->Add(hElossF120HEF);
	hElossP->Add(hElossP120HEF);
	
	// hElossT->Add(hElossT120HEB);
	// hElossF->Add(hElossF120HEB);
	// hElossP->Add(hElossP120HEB);

      }else {

	hElossT->Add(hElossT120EE);
	hElossF->Add(hElossF120EE);
	hElossP->Add(hElossP120EE);

	hElossT->Add(hElossT120HEF);
	hElossF->Add(hElossF120HEF);
	hElossP->Add(hElossP120HEF);
	
	// hElossT->Add(hElossT120HEB);
	// hElossF->Add(hElossF120HEB);
	// hElossP->Add(hElossP120HEB);

      }

    }//is HE
    
  }//file loop

  hElossT->SetLineColor(kBlack);
  hElossF->SetLineColor(kGreen+1);
  hElossP->SetLineColor(kRed);

  hElossT->SetLineWidth(3);
  hElossF->SetLineWidth(3);
  hElossP->SetLineWidth(3);

  if(refType==1)
    hElossT->SetTitle("Energy loss for 120 #mum Si");
  else if(refType==2)
    hElossT->SetTitle("Energy loss for 200 #mum Si");
  else
    hElossT->SetTitle("Energy loss for 300 #mum Si");
  
  hElossF->SetName("Full wafers");
  hElossP->SetName("Partial wafers");

  TLegend *leg1 = new TLegend(0.56,0.57,0.99,0.92);
  leg1->AddEntry(hElossT,"Inclusive", "L");
  leg1->AddEntry(hElossF,"Full wafers", "L");
  leg1->AddEntry(hElossP,"Partial wafers", "L");

  TF1 *fnF = new TF1("fnF","landau",0,500.0);
  fnF->SetNpx(1000);
  fnF->SetLineColor(kGreen+1);
  fnF->SetLineWidth(3);
  fnF->SetParameter(1,hElossF->GetMean());
  fnF->SetParameter(2,hElossF->GetRMS());

  TF1 *fnP = (TF1 *)fnF->Clone("fnP");
  fnP->SetLineColor(kRed);
  fnP->SetParameter(1,hElossP->GetMean());
  fnP->SetParameter(2,hElossP->GetRMS());

  hElossF->Fit(fnF,"");
  hElossP->Fit(fnP,"");
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  // c1->Divide(2,1);
  // c1->cd(1);
  hElossT->Draw();
  hElossF->Draw("sames");

  hElossP->Draw("sames");
  
  leg1->Draw();

  
  hElossT->GetXaxis()->SetRangeUser(0.0, 200.0);
  hElossT->GetXaxis()->SetTitle("Eloss (keV)");
  hElossT->GetYaxis()->SetTitle("Entries");
  
  // c1->Update();

  // TPaveStats *stats0 = (TPaveStats *)hEPRectotCN0->FindObject("stats");
  // TPaveStats *stats1 = (TPaveStats *)hEPRectotCK0->FindObject("stats");
  // TPaveStats *stats2 = (TPaveStats *)hEPRecFailtotCN0->FindObject("stats");
  // TPaveStats *stats3 = (TPaveStats *)hEPRecFailtotCK0->FindObject("stats");
  // stats0->SetX1NDC(.56);          stats0->SetX2NDC(.76);
  // stats0->SetY1NDC(.4);           stats0->SetY2NDC(.55);
  // stats1->SetX1NDC(.78);           stats1->SetX2NDC(.99);
  // stats1->SetY1NDC(.4);           stats1->SetY2NDC(.55);
  // stats2->SetX1NDC(.56);          stats2->SetX2NDC(.76);
  // stats2->SetY1NDC(.2);           stats2->SetY2NDC(.35);
  // stats3->SetX1NDC(.78);           stats3->SetX2NDC(.99);
  // stats3->SetY1NDC(.2);           stats3->SetY2NDC(.35);
  // stats3->Draw();
  
  // c1->cd(1);
  // TPaveStats *stats4 = (TPaveStats *)hEPGentotCN0->FindObject("stats");
  // TPaveStats *stats5 = (TPaveStats *)hEPGentotCK0->FindObject("stats");
  // TPaveStats *stats6 = (TPaveStats *)hEPGenFailtotCN0->FindObject("stats");
  // TPaveStats *stats7 = (TPaveStats *)hEPGenFailtotCK0->FindObject("stats");
  // stats4->SetX1NDC(.56);          stats4->SetX2NDC(.76);
  // stats4->SetY1NDC(.4);           stats4->SetY2NDC(.55);
  // stats5->SetX1NDC(.78);           stats5->SetX2NDC(.99);
  // stats5->SetY1NDC(.4);           stats5->SetY2NDC(.55);
  // stats6->SetX1NDC(.56);          stats6->SetX2NDC(.76);
  // stats6->SetY1NDC(.2);           stats6->SetY2NDC(.35);
  // stats7->SetX1NDC(.78);           stats7->SetX2NDC(.99);
  // stats7->SetY1NDC(.2);           stats7->SetY2NDC(.35);
  // stats7->Draw();
  
  
  return true;
}

// int AddGraphs(TGraph*& gr1, TGraph*& gr2)
// {
//   int npoints = gr1->GetN();

//   for(int i=0;i<gr2->GetN();i++)
//     gr1->SetPoint(npoints++, gr2->GetPointX(i), gr2->GetPointY(i));

//   return true;
// }
