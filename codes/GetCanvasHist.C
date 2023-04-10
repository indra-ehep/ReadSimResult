/**********************************************************************
 Created on : 09/04/2023
 Purpose    : Get the histogram from saved canvas
 Author     : Indranil Das, Visiting Fellow
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

int GetCanvasHist(const char *infile = "CombinedCSMax_Si300_D83_1M.root", const char *cName = "c3_3")
{
  TFile *fin = TFile::Open(infile);
  TCanvas *c1 = (TCanvas *)fin->Get(cName);
  TH1D *hTot = (TH1D *)c1->GetListOfPrimitives()->At(1);
  gStyle->SetOptStat(0);
  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  hTot->Draw();
  gStyle->SetOptStat(0);
  hTot->Draw();
  
  return true;
}
