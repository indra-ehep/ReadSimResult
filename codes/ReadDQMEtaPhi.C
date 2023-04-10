/**********************************************************************
 Created on : 10/04/2023
 Purpose    : Read the eta-phi plots prepared by DQM
 Author     : Indranil Das, Visiting Fellow
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

int ReadDQMEtaPhi(const char *infile = "DQM_plots.root", const char *hName = "h2_HEScint_layer21_simHits_12_5_0_pre1_default", const char* lTitle = "12_5_0_pre1_Default")
{
  TFile *fin = TFile::Open(infile);
  TH2F *h2 = (TH2F *)fin->Get(hName);
  //h2->SetTitle("Release validation Simhits : Layer 41");
  h2->SetTitle("Release validation Rechits : Layer 41");
  
  TLegend *leg1 = new TLegend(0.56,0.75,0.99,0.88);
  leg1->AddEntry(h2,Form("%s",lTitle), "L");

  Color_t color = kRed;
  h2->SetLineColor(color);
  h2->SetLineWidth(3);
  h2->SetMarkerColor(color);

  //h2->Rebin2D(2);
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  h2->Draw();
  leg1->Draw();
  h2->GetXaxis()->SetRangeUser(-3.0, -1.5);
  h2->GetXaxis()->SetTitle("#eta");
  h2->GetYaxis()->SetTitle("#phi");
  h2->GetXaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleSize(0.06);
  h2->GetYaxis()->SetTitleOffset(0.8);

  return true;
}
