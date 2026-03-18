#include "Fit_Diffrazione.h"

/*****************************************************************************************
 *****************************************************************************************
 ** **
 ** Questo programma esegue il fit della funzione che descrive la figura di
diffrazione **
 ** prodotta da una fenditura lineare. Due funzioni accessorie permettono la
visualiz-  **
 ** zazione dei dati sperimentali (utile per stimare i valori iniziali delle
grandezze  **
 ** coinvolte) e di produrre il grafico della funzione da fittare. **
 ** **
 *****************************************************************************************
  C.M. 8 Giugno 2015
  C.M. 4 Febbraio 2016, aggiunta del fondo par[5]
  C.M. 14 Marzo 2024, aggiornamento linea 48, char --> Tstring
 *****************************************************************************************
Si possono utilizzare questi comandi da terminale
> root -l --> Avvia ROOT
root[0] .L Fit_Diffrazione.cpp  --> Carica il file Fit_Diffrazione.cpp
root[1] mydata("nomefile.txt") --> costruisce il grafico dei dati sperimentali
root[2] myfunc() --> visualizza la funzione con valori di default
root[3] myfunc(background, normalizzazione, lambda, larghezza fenditura, shift
lungo x, distanza fenditura - schermo) --> visualizza la funzione con valori
forniti dall'utente myfit("nomefile.txt", background, normalizzazione, lambda,
larghezza fenditura, shift lungo x, distanza fenditura - schermo) --> esegue il
fit Tutte le grandezze sono in metri.
 ****************************************************************************************/

// Double_t Diffrazione(double *x, double *par) {
//   double Arg =
//       TMath::Pi() * par[0] * (x[0] - par[1]) /
//       TMath::Sqrt((x[0] - par[1]) * (x[0] - par[1]) + par[2] * par[2]) / par[3];
//   double Diffr = par[5] + par[4] * sin(Arg) * sin(Arg) / Arg / Arg;
//   return Diffr;
// }

Double_t Diffrazione(double *x, double *par) {
    double dx = x[0] - par[1]; // x - x0
    double Arg = TMath::Pi() * par[0] * dx / (par[2] * par[3]); // pi*a*(x-x0)/(L*lambda)

    double sinc2;
    if (fabs(Arg) < 1e-6)
        sinc2 = 1.0; // limite centrale
    else
        sinc2 = TMath::Power(TMath::Sin(Arg)/Arg, 2);

    double Diffr = par[5] + par[4] * sinc2;
    return Diffr;
}


void myfunc(const FitParams &p) {
  TF1 *f1 = new TF1("myfunc", Diffrazione, p.x0 - 0.03, p.x0 + 0.03, 6);
  f1->SetParameter(0, p.d);
  f1->SetParameter(1, p.x0);
  f1->SetParameter(2, p.L);
  f1->SetParameter(3, p.lambda);
  f1->SetParameter(4, p.I0);
  f1->SetParameter(5, p.bkg);
  f1->SetNpx(1000);
  f1->SetParName(0, "Larghezza fenditura");
  f1->SetParName(1, "shift lungo x");
  f1->SetParName(2, "Distanza fenditura - schermo");
  f1->SetParName(3, "Lambda");
  f1->SetParName(4, "Normalizzazione");
  f1->SetParName(5, "fondo");
  f1->Draw("same");
}

void mydata(const FitParams& p) {
  TCanvas *c1 = new TCanvas("c1", "Diffrazione", 800, 600);
  c1->SetFillColor(0);
  c1->SetGrid();
  TGraphErrors *data = new TGraphErrors(p.fname, "%lg %lg %lg");
  data->SetTitle("Figura di diffrazione");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.5);
  data->SetMarkerColor(kBlue + 2);
  data->SetLineColor(kBlue + 2);
  data->SetLineWidth(2);
  data->GetXaxis()->SetTitle("Posizione [m]");
  data->GetYaxis()->SetTitle("Voltaggio (V)");
  data->GetXaxis()->CenterTitle();
  data->GetYaxis()->CenterTitle();
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetYaxis()->SetTitleSize(0.05);
  data->GetXaxis()->SetLabelSize(0.04);
  data->GetYaxis()->SetLabelSize(0.04);
  data->Draw("AP");
  c1->SaveAs("Diffrazione.pdf");
}

void myfit(const FitParams &p) {
  TGraphErrors *data = new TGraphErrors(p.fname, "%lg %lg %lg");
  TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
  f1->SetParameter(0, p.d);
  f1->SetParameter(1, p.x0);
  f1->SetParameter(2, p.L);
  f1->SetParameter(3, p.lambda);
  f1->SetParameter(4, p.I0);
  f1->SetParameter(5, p.bkg);
  f1->FixParameter(0, p.d);
  f1->FixParameter(2, p.L);

  // f1->SetParLimits(1, x0 - 0.001, x0 + 0.001);
  // f1->SetParLimits(4, I0 - 10., I0 + 10.);
  // f1->SetParLimits(5, bkg - 5., bkg + 5.);

  TCanvas *c2 = new TCanvas("c2", "Diffrazione_fit", 800, 600);
  data->Fit("myfunc", "R");
  data->Draw("AP");
  data->SetLineColor(4);
  data->SetMarkerColor(4);
  f1->Draw("same");
  data->SetTitle("Figura di diffrazione");
  data->GetXaxis()->SetTitle("Posizione, m");
  data->GetYaxis()->SetTitle("Voltaggio (V)");
  data->GetXaxis()->CenterTitle(true);
  data->GetXaxis()->CenterTitle(true);
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0); 
  leg->SetFillColor(0);  
  leg->AddEntry(data, "L= ... m, d=... mm", "p");
  leg->AddEntry(f1, "fit", "l");
  leg->Draw();
  c2->SaveAs("Diffrazione_fit.pdf");
}
