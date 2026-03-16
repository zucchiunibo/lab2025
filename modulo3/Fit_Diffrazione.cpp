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

Double_t Diffrazione(double *x, double *par) {
  double Arg =
      TMath::Pi() * par[0] * (x[0] - par[1]) /
      TMath::Sqrt((x[0] - par[1]) * (x[0] - par[1]) + par[2] * par[2]) / par[3];
  double Diffr = par[5] + par[4] * sin(Arg) * sin(Arg) / Arg / Arg;
  return Diffr;
}

void myfunc(double bkg = 0., double I0 = 1., double lambda = 600.E-9,
            double d = 1.E-4, double x0 = 0.05, double L = 1.0) {
  TF1 *f1 = new TF1("myfunc", Diffrazione, x0 - 0.03, x0 + 0.03, 6);
  f1->SetParameter(0, d);
  f1->SetParameter(1, x0);
  f1->SetParameter(2, L);
  f1->SetParameter(3, lambda);
  f1->SetParameter(4, I0);
  f1->SetParameter(5, bkg);
  f1->SetParName(0, "Larghezza fenditura");
  f1->SetParName(1, "shift lungo x");
  f1->SetParName(2, "Distanza fenditura - schermo");
  f1->SetParName(3, "Lambda");
  f1->SetParName(4, "Normalizzazione");
  f1->SetParName(5, "fondo");
  f1->Draw("same");
}

void mydata(TString fname = "diff1(step100)") {
  TGraph *data = new TGraphErrors(fname, "%lg %lg");
  data->Draw("AP");
  data->SetLineColor(4);
  data->SetMarkerColor(4);
  data->SetTitle("Figura di diffrazione");
  data->GetXaxis()->SetTitle("Posizione, m");
  data->GetYaxis()->SetTitle("Int. luminosa, unit. arb.");
  data->GetXaxis()->CenterTitle(true);
  data->GetXaxis()->CenterTitle(true);
  data->SaveAs("Grafico.jpg");
}

void myfit(TString fname = "diff1(step100)", double bkg = 0., double I0 = 1.,
           double lambda = 600.E-9, double d = 1.E-4, double x0 = 0.05,
           double L = 1.0) {
  TGraphErrors *data = new TGraphErrors(fname, "%lg %lg");
  TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
  f1->SetParameter(0, d);
  f1->SetParameter(1, x0);
  f1->SetParameter(2, L);
  f1->SetParameter(3, lambda);
  f1->SetParameter(4, I0);
  f1->SetParameter(5, bkg);

  f1->FixParameter(0, d);
  f1->FixParameter(2, L);

       f1->SetParLimits(1,x0-0.001,x0+0.001);
       f1->SetParLimits(4,I0-10., I0+10.);
  //     f1->SetParLimits(5,bkg-5., bkg+5.);

  data->Fit("myfunc", "R");
  data->Draw("AP");
  data->SetLineColor(4);
  data->SetMarkerColor(4);
  f1->Draw("same");
  data->SetTitle("Figura di diffrazione");
  data->GetXaxis()->SetTitle("Posizione, m");
  data->GetYaxis()->SetTitle(" ");
  data->GetXaxis()->CenterTitle(true);
  data->GetXaxis()->CenterTitle(true);
  TLegend *leg = new TLegend(.6, .7, .9, .9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0); // no border for legend
  leg->SetFillColor(0);  // fill color is white
  leg->AddEntry(data, "L= ... m, d=... mm", "p");
  leg->AddEntry(f1, "fit", "l");
  leg->Draw();
}
