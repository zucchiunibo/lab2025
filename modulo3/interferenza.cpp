#include "Fit_Diffrazione.h"

void grafico(TString filename) {
  TCanvas *c3 = new TCanvas("c3", "Interferenza", 800, 600);
  c3->SetFillColor(0);
  c3->SetLeftMargin(0.12);
  c3->SetBottomMargin(0.12);
  c3->SetRightMargin(0.05);
  c3->SetTopMargin(0.08);

  c3->SetGrid();
  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(3);
  TGraph *interferenza = new TGraph(filename);
  interferenza->SetTitle("Grafico di Interferenza");
  interferenza->SetMarkerStyle(20);
  interferenza->SetMarkerSize(0.8);
  interferenza->SetMarkerColor(kBlue + 2);
  interferenza->SetLineColor(kBlue + 2);
  interferenza->SetLineWidth(2);
  interferenza->GetXaxis()->SetTitle("Posizione [m]");
  interferenza->GetYaxis()->SetTitle("Voltaggio (V)");
  interferenza->GetXaxis()->CenterTitle();
  interferenza->GetYaxis()->CenterTitle();
  interferenza->GetXaxis()->SetTitleSize(0.05);
  interferenza->GetYaxis()->SetTitleSize(0.05);
  interferenza->GetXaxis()->SetLabelSize(0.04);
  interferenza->GetYaxis()->SetLabelSize(0.04);
  interferenza->Draw("APL");
  c3->Update();
  interferenza->GetXaxis()->SetRangeUser(0.031, 0.045); //Commentare
  c3->Modified();
  c3->Update();
  c3->SaveAs("Interferenza.pdf");
}