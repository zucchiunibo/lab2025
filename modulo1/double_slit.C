#include "TGraph.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <vector>

class double_slit {
  double B_;
  double I_0_;
  double k_;
  double a_;
  double d_;

 public:
  double_slit()
      : B_{0.0}
      , I_0_{0.0}
      , k_{0.0}
      , a_{0.0}
      , d_{0.0} {
  }

  double_slit(double B, double I_0, double k, double a, double d)
      : B_{B}
      , I_0_{I_0}
      , k_{k}
      , a_{a}
      , d_{d} {
  }

  void fit() {
    TF1* f     = new TF1("Funzione da simulare", "[0] + [1] * ((sin([2]*[3]*x))/([2]*[3]*x)) * (cos([2]*[4]*x))^2");
    TH1F* hist = new TH1F("hist", "Istogramma delle occorrenze", 400, -0.05f, 0.05f);

    f->SetParameters(30, 150, 4, 10, 157/4);  //TOP: (30, 150, 4, 10, 157 / 4) (0, 1, 2, 3, 4)
    std::cout << f->GetParameter(1);
    std::ifstream inputFile("data_double_slit.txt");
    double value;

    if (inputFile.is_open()) {
      while (inputFile >> value) {
        hist->Fill(value);
      }
      inputFile.close();
    } else {
      std::cerr << "Errore: Impossibile aprire il file dati.txt" << '\n';
      return;
    }

    hist->Fit(f);

    TCanvas* DS = new TCanvas("DS", "Hiso_DB_Slit", 800, 600);
    hist->Draw();
    f->Draw("same");
    DS->SaveAs("Double_slit.png");
  }
};