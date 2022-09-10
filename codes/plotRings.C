/**********************************************************************
 Created on : 05/09/2022
 Purpose    : Plot the rings of scintillator tiles
 Author     : Indranil Das, Visiting Fellow
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <bitset>
#include <string>

int GetCartesian(double r, double theta, double& x, double& y){
  x = r*TMath::Cos(theta*TMath::Pi()/180.);
  y = r*TMath::Sin(theta*TMath::Pi()/180.);
  return true;
}

int plotRings(int refLayer=34)
{
  
  ifstream fin("Validation/HGCalValidation/data/scintillatorV0.txt");
  string s;
  vector<string> SciLayouts;
  while(getline(fin,s)){
    //cout << s << ",  SIZE : "<<s.size()<< endl;
    if(s.size()==55){
      //cout << s << endl;
      SciLayouts.push_back(s);
    }
  }
  // TFile *fGeant = TFile::Open("geantoutput_merged_D49.root");
  // TH2D *hXYhits = (TH2D *)fGeant->Get(Form("prodHEB/hXYhits_layer_%d",refLayer));

  TFile *fGeant = TFile::Open("geantoutputD88_gr-merged_old_cmssw.root");
  TH2D *hXYhits = (TH2D *)fGeant->Get(Form("hgcalCellHitSumHEB/hXYhits_layer_%d",refLayer));
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gStyle->SetOptStat(100011);
  hXYhits->SetMarkerColor(kBlue);
  hXYhits->Draw("sames");

  // TFile *fGeant = TFile::Open("geantoutputD88_gr-merged_old_cmssw.root");
  // TGraph *hXYhits = (TGraph *)fGeant->Get(Form("hgcalCellHitSumHEB/gXYhitsB1_layer_%d",refLayer));
  // TCanvas *c1 = new TCanvas("c1","c1",600,600);
  // hXYhits->SetTitle(Form("Hits in layer %d (z > 0.0 cm)",refLayer));
  // hXYhits->SetMarkerColor(kBlue);
  // hXYhits->SetMarkerStyle(kFullCircle);
  // hXYhits->SetMarkerSize(0.1);
  // hXYhits->Draw("AP");

  for(int ilir=0;ilir<SciLayouts.size();ilir++){
    int layer, ring;
    float rmin, rmax, sipmsize;
    long int hex1, hex2, hex3, hex4;
    string hexcomb;
    char scitype;
    sscanf((const char *)SciLayouts.at(ilir).c_str(), "%d %d %f %f %f %lX %lX %lX %lX %c\n", &layer, &ring, &rmin, &rmax, &sipmsize, &hex1, &hex2, &hex3, &hex4, &scitype);
    if(layer==refLayer){
      printf("%d %d %6.2f %6.2f %3.1f %lX %lX %lX %lX %c\n",  layer,  ring,  rmin,  rmax,  sipmsize,  hex1,  hex2,  hex3,  hex4,  scitype);
      //hexcomb = Form("%lX%lX%lX%lX",hex1,  hex2,  hex3,  hex4);
      hexcomb = std::bitset<24>(hex1).to_string() + std::bitset<24>(hex2).to_string() + std::bitset<24>(hex3).to_string() + std::bitset<24>(hex4).to_string() ;
      //cout<<bitset<24>(hex1)<<endl;
      cout<<hexcomb<<endl;
      rmax /= 10.0;
      rmin /= 10.0;
      //if(ring==29){
      for(int isect=0;isect<3;isect++){
	for(int iphi=0;iphi<hexcomb.size();iphi++){
	  float phimin = isect*120. + iphi*1.25;
	  float phimax = isect*120. + (iphi+1)*1.25;
	  bool isDraw = (hexcomb[iphi]=='1') ? true : false; 
	  if(isDraw){
	    
	    TArc *arc1 = new TArc(0.0, 0.0, rmax, phimin,phimax);
	    arc1->SetFillStyle(0);
	    arc1->SetNoEdges();
	    arc1->SetLineWidth(1);
	    arc1->Draw("sames");
	    
	    double x1,y1,x2,y2;
	    GetCartesian(rmax, phimin, x1, y1);
	    // printf("x1,y1 : %f,%f\n",x1,y1);
	    GetCartesian(rmin, phimin, x2, y2);
	    // printf("x2,y2 : %f,%f\n",x2,y2);
	    TLine *l1 = new TLine(x1,y1,x2,y2);
	    //l1->SetLineColor(kBlue);
	    //l1->SetLineWidth(2);
	    l1->SetLineWidth(1);
	    l1->Draw("sames");

	    GetCartesian(rmax, phimax, x1, y1);
	    // printf("x1,y1 : %f,%f\n",x1,y1);
	    GetCartesian(rmin, phimax, x2, y2);
	    // printf("x2,y2 : %f,%f\n",x2,y2);
	    TLine *l2 = new TLine(x1,y1,x2,y2);
	    //l1->SetLineColor(kBlue);
	    //l1->SetLineWidth(2);
	    l2->SetLineWidth(1);
	    l2->Draw("sames");

	    TArc *arc2 = new TArc(0.0, 0.0, rmin, phimin,phimax);
	    arc2->SetFillStyle(0);
	    arc2->SetNoEdges();
	    arc2->SetLineWidth(1);
	    arc2->Draw("sames");
	    
	  }
	}
      }
    }
  }
  fin.close();
  SciLayouts.clear();
  c1->SetTickx();
  c1->SetTicky();
  hXYhits->GetXaxis()->SetTitle("x (cm)");
  hXYhits->GetYaxis()->SetTitleOffset(1.4);
  hXYhits->GetYaxis()->SetTitle("y (cm)");
  c1->SaveAs(Form("Scintillator_merged_Layer_%d.png",refLayer));
  c1->SaveAs(Form("Scintillator_merged_Layer_%d.pdf",refLayer));
  //c1->SaveAs(Form("Scintillator_merged_Layer_%d.svg",refLayer));
  delete fGeant;
  
  return true;
}
