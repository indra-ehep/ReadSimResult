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

int GetPolar(double x, double y, double& r, double& theta){
  r = TMath::Sqrt(x*x + y*y);
  theta =  180.0*TMath::ATan2(y, x)/TMath::Pi();
  if(theta<0.0)
    theta = (180. - TMath::Abs(theta))  + 180.;
  return true;
}

int plotRings(int refLayer=47)
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

  // TFile *fGeant = TFile::Open("geantoutputD88_gr-merged_ReEval.root");
  // TH2D *hXYhits = (TH2D *)fGeant->Get(Form("hgcalCellHitSumHEB/hXYhits_layer_%d",refLayer));
  
  // TCanvas *c1 = new TCanvas("c1","c1",600,600);
  // gStyle->SetOptStat(100011);
  // hXYhits->SetMarkerColor(kBlue);
  // hXYhits->Draw("sames");

  
  //TFile *fGeant = TFile::Open("geantoutputD92_val_cassette_shift.root");
  TFile *fGeant = TFile::Open("geantoutputD92_gr-merged.root");
  TGraph *hXYhits = (TGraph *)fGeant->Get(Form("hgcalCellHitSumHEB/gXYhitsB0_layer_%d",refLayer));
  //TGraph *hXYhits = (TGraph *)fGeant->Get(Form("hgcalCellHitSumHEF/gXYhitsCK1_layer_%d",refLayer));
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hXYhits->SetTitle(Form("Hits in layer %d (z > 0.0 cm)",refLayer));
  hXYhits->SetMarkerColor(kBlue);
  //hXYhits->SetMarkerColor(kMagenta);
  hXYhits->SetMarkerStyle(kFullCircle);
  hXYhits->SetMarkerSize(0.1);
  hXYhits->Draw("AP");

  float xshift = -60.0, yshift = -20.0;
  int tcassette = 12;
  for(int ilir=0;ilir<SciLayouts.size();ilir++){
    int layer, ring;
    float rmin, rmax, sipmsize;
    long int hex1, hex2, hex3, hex4;
    string hexcomb;
    char scitype;
    int icassette = 0 ;
    sscanf((const char *)SciLayouts.at(ilir).c_str(), "%d %d %f %f %f %lX %lX %lX %lX %c\n", &layer, &ring, &rmin, &rmax, &sipmsize, &hex1, &hex2, &hex3, &hex4, &scitype);
    if(layer==refLayer){
      printf("%d %d %6.2f %6.2f %3.1f %lX %lX %lX %lX %c\n",  layer,  ring,  rmin,  rmax,  sipmsize,  hex1,  hex2,  hex3,  hex4,  scitype);
      //hexcomb = Form("%lX%lX%lX%lX",hex1,  hex2,  hex3,  hex4);
      hexcomb = std::bitset<24>(hex1).to_string() + std::bitset<24>(hex2).to_string() + std::bitset<24>(hex3).to_string() + std::bitset<24>(hex4).to_string() ;
      //cout<<bitset<24>(hex1)<<endl;
      cout<<hexcomb<<endl;
      rmax /= 10.0;
      rmin /= 10.0;
      // xshift /= 10.0;
      // yshift /= 10.0;
      //if(ring==29){
      int iseg = 0 ;
      for(int isect=0;isect<3;isect++){
	cout<<" hexcomb.size() : " << hexcomb.size() << endl;
	for(int iphi=0;iphi<hexcomb.size();iphi++){

	  double phimin = isect*120. + iphi*1.25;
	  double phimax = isect*120. + (iphi+1)*1.25;
	  double nrmin = rmin ;
	  double nrmax = rmax ;
	  
	  if(iseg%24==0)
	    icassette++;
	  //cout << "LProcesing cassette : " <<icassette<<endl;
	  double offsetX = 0.0;
	  double offsetY = 0.0;
	  double offsetR = 0.0;
	  double offsetT = 0.0;
	  if(icassette==tcassette){
	    offsetX = xshift;
	    offsetY = yshift;
	    double x1,y1,x2,y2,x3,y3,x4,y4;
	    GetCartesian(rmin, phimin, x1, y1); x1 += xshift; y1 += yshift; 
	    GetCartesian(rmin, phimax, x2, y2); x2 += xshift; y2 += yshift; 
	    GetCartesian(rmax, phimin, x3, y3); x3 += xshift; y3 += yshift; 
	    GetCartesian(rmax, phimax, x4, y4); x4 += xshift; y4 += yshift; 
	    
	    GetPolar(x1, y1, nrmin, phimin);
	    GetPolar(x2, y2, nrmin, phimax);
	    GetPolar(x3, y3, nrmax, phimin);
	    GetPolar(x4, y4, nrmax, phimax);
	    
	    GetPolar(offsetX, offsetY, offsetR, offsetT);
	    cout << "offsetX " << offsetX << ", offsetY " << offsetY << ", offsetR " << offsetR << ", offsetT " << offsetT << endl;
	  }else{
	    offsetX = 0.0;
	    offsetY = 0.0;
	  }
	  
	  
	  bool isDraw = (hexcomb[iphi]=='1') ? true : false;
	  if(isDraw){
	    
	    TArc *arc1 = new TArc(offsetX, offsetY, nrmax, phimin,phimax);
	    arc1->SetFillStyle(0);
	    arc1->SetNoEdges();
	    arc1->SetLineWidth(1);
	    arc1->Draw("sames");
	    
	    double x1,y1,x2,y2;
	    GetCartesian(nrmax, phimin, x1, y1);
	    // printf("x1,y1 : %f,%f\n",x1,y1);
	    GetCartesian(nrmin, phimin, x2, y2);
	    // printf("x2,y2 : %f,%f\n",x2,y2);
	    x1 += offsetX; y1 += offsetY;
	    x2 += offsetX; y2 += offsetY; 
	    TLine *l1 = new TLine(x1,y1,x2,y2);
	    //l1->SetLineColor(kBlue);
	    //l1->SetLineWidth(2);
	    l1->SetLineWidth(1);
	    l1->Draw("sames");

	    GetCartesian(nrmax, phimax, x1, y1);
	    // printf("x1,y1 : %f,%f\n",x1,y1);
	    GetCartesian(nrmin, phimax, x2, y2);
	    // printf("x2,y2 : %f,%f\n",x2,y2);
	    x1 += offsetX; y1 += offsetY;
	    x2 += offsetX; y2 += offsetY; 
	    TLine *l2 = new TLine(x1,y1,x2,y2);
	    //l1->SetLineColor(kBlue);
	    //l1->SetLineWidth(2);
	    l2->SetLineWidth(1);
	    l2->Draw("sames");

	    TArc *arc2 = new TArc(offsetX, offsetY, nrmin, phimin,phimax);
	    arc2->SetFillStyle(0);
	    arc2->SetNoEdges();
	    arc2->SetLineWidth(1);
	    arc2->Draw("sames");
	    
	  }//isDraw
	  iseg++;
	}
      }
    }
  }
  fin.close();
  SciLayouts.clear();
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx();
  c1->SetGridy();
  hXYhits->GetXaxis()->SetTitle("x (cm)");
  hXYhits->GetYaxis()->SetTitleOffset(1.4);
  hXYhits->GetYaxis()->SetTitle("y (cm)");
  c1->SaveAs(Form("Scintillator_merged_Layer_%d.png",refLayer));
  c1->SaveAs(Form("Scintillator_merged_Layer_%d.pdf",refLayer));
  //c1->SaveAs(Form("Silicon_merged_Layer_%d.png",refLayer));
  //c1->SaveAs(Form("Silicon_merged_Layer_%d.pdf",refLayer));
  //c1->SaveAs(Form("Scintillator_merged_Layer_%d.svg",refLayer));
  delete fGeant;
  
  return true;
}
